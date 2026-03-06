import streamlit as st
import pandas as pd
import numpy as np
import io
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from pybiomart import Dataset
from typing import Optional
from collections import Counter

# =========================================================
# Streamlit setup
# =========================================================
st.set_page_config(
    page_title="Chromosome Amino Acid Composition and E/Q Ratio Explorer",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 Chromosome Amino Acid Composition and E/Q Ratio Explorer (Streamlit + Ensembl)")

st.write("""
Enter a **human chromosome (1–22, X, Y)**.

Button-driven workflow (fast UI + heavy computations only when requested):
1) Generate the **raw Ensembl protein table** (Swiss-Prot mapping + peptide) — show preview + downloads
2) Generate the **E/Q analysis** (AA counts + E/Q + statistics + plots) — show preview + downloads
3) Generate the **Walking (Surfing) analysis** (AA selected + E + Q only) + quick chromosome-wide plot
""")

chromosome = st.text_input("Chromosome (1–22, X, Y):", "").strip()
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

# =========================================================
# Session state init
# =========================================================
for k in [
    "chrom_ok",
    "ensembl_df",
    "df_all",
    "df_valid",
    "eq_summary_df",
    "avg_comp_all",
    "walk_df",
    "fig_corr",
    "fig_hist",
    "fig_walk",
]:
    if k not in st.session_state:
        st.session_state[k] = None


# =========================================================
# Amino Acid Reference Table (collapsed by default)
# =========================================================
with st.expander("📘 Amino Acid Reference Table (click to view)", expanded=False):
    aa_data = [
        ["Alanine", "Ala", "A"],
        ["Arginine", "Arg", "R"],
        ["Asparagine", "Asn", "N"],
        ["Aspartic acid", "Asp", "D"],
        ["Cysteine", "Cys", "C"],
        ["Glutamine", "Gln", "Q"],
        ["Glutamic acid", "Glu", "E"],
        ["Glycine", "Gly", "G"],
        ["Histidine", "His", "H"],
        ["Isoleucine", "Ile", "I"],
        ["Leucine", "Leu", "L"],
        ["Lysine", "Lys", "K"],
        ["Methionine", "Met", "M"],
        ["Phenylalanine", "Phe", "F"],
        ["Proline", "Pro", "P"],
        ["Serine", "Ser", "S"],
        ["Threonine", "Thr", "T"],
        ["Tryptophan", "Trp", "W"],
        ["Tyrosine", "Tyr", "Y"],
        ["Valine", "Val", "V"],
    ]
    aa_df_ref = pd.DataFrame(aa_data, columns=["Full Name", "3-letter", "1-letter"])
    st.dataframe(aa_df_ref, hide_index=True, use_container_width=True)


# =========================================================
# Data fetch (Ensembl / BioMart) -- keep http
# =========================================================
@st.cache_data(show_spinner=False)
def fetch_ensembl_chromosome_data(chrom: str) -> pd.DataFrame:
    dataset = Dataset(name="hsapiens_gene_ensembl", host="http://www.ensembl.org")
    df = dataset.query(
        attributes=[
            "ensembl_gene_id",
            "external_gene_name",
            "start_position",
            "end_position",
            "chromosome_name",
            "uniprotswissprot",
            "peptide",
        ],
        filters={"chromosome_name": [chrom]},
    )

    df = df.rename(
        columns={
            "ensembl_gene_id": "Gene stable ID",
            "external_gene_name": "Gene name",
            "start_position": "Gene start (bp)",
            "end_position": "Gene end (bp)",
            "chromosome_name": "Chromosome/scaffold name",
            "uniprotswissprot": "Entry",
            "peptide": "Sequence",
            # friendly alt names (if returned)
            "UniProtKB/Swiss-Prot ID": "Entry",
            "Peptide": "Sequence",
        }
    )

    df["Entry"] = df["Entry"].replace("", np.nan)
    df["Sequence"] = df["Sequence"].replace("", np.nan)
    df = df.dropna(subset=["Entry", "Sequence"]).copy()

    df["Sequence"] = df["Sequence"].astype(str).str.rstrip("*")
    df["Gene start (bp)"] = pd.to_numeric(df["Gene start (bp)"], errors="coerce")
    df["Gene end (bp)"] = pd.to_numeric(df["Gene end (bp)"], errors="coerce")

    df = df.sort_values(by=["Gene start (bp)"], na_position="last").reset_index(drop=True)
    return df


# =========================================================
# FAST AA counting + E/Q (no BioPython)
# =========================================================
def fast_count_aa(seq: str) -> dict:
    # Keep only uppercase letters; remove terminal stop marker if any
    seq = str(seq).rstrip("*").upper()
    c = Counter(seq)
    return {aa: int(c.get(aa, 0)) for aa in AA_ORDER}


def add_aa_counts_and_eq_fast(df: pd.DataFrame, batch_size: int = 800) -> pd.DataFrame:
    """
    Much faster than BioPython ProteinAnalysis for thousands of sequences.
    """
    n = len(df)
    progress = st.progress(0)
    status = st.empty()

    rows = []
    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        batch = df.iloc[start:end]

        batch_out = []
        for seq in batch["Sequence"].astype(str):
            counts = fast_count_aa(seq)
            e = counts["E"]
            q = counts["Q"]
            eq = (e / q) if q != 0 else np.nan
            batch_out.append(
                {
                    **counts,
                    "Length": len(str(seq).rstrip("*")),
                    "E/Q ratio": eq,
                    "E/Q valid": (q != 0),
                }
            )

        rows.append(pd.DataFrame(batch_out))
        progress.progress(end / n)
        status.write(f"Computed {end}/{n} proteins...")

    progress.empty()
    status.empty()

    aa_df = pd.concat(rows, ignore_index=True)
    out = pd.concat([df.reset_index(drop=True), aa_df], axis=1)
    return out


# =========================================================
# Walking / Surfing analysis (AA selected + E + Q ONLY)
# =========================================================
def make_walk_table(
    df_all_sorted: pd.DataFrame,
    window: int,
    threshold: float,
    aa_selected: str,
    accession_col: str = "Entry",
    pos_col: str = "Gene start (bp)",
) -> pd.DataFrame:
    n = len(df_all_sorted)
    if n < window:
        return pd.DataFrame()

    needed = [aa_selected, "E", "Q", accession_col, pos_col]
    missing = [c for c in needed if c not in df_all_sorted.columns]
    if missing:
        raise ValueError(f"Missing required columns for walking analysis: {missing}")

    aa_vals = pd.to_numeric(df_all_sorted[aa_selected], errors="coerce").to_numpy()
    e_vals = pd.to_numeric(df_all_sorted["E"], errors="coerce").to_numpy()
    q_vals = pd.to_numeric(df_all_sorted["Q"], errors="coerce").to_numpy()
    acc = df_all_sorted[accession_col].astype(str).to_numpy()
    pos = pd.to_numeric(df_all_sorted[pos_col], errors="coerce").to_numpy()

    rows = []
    for i in range(0, n - window + 1):
        j = i + window

        frame_aa = aa_vals[i:j]
        frame_e = e_vals[i:j]
        frame_q = q_vals[i:j]
        frame_acc = acc[i:j]
        frame_pos = pos[i:j]

        surfing_aa = int(np.sum(frame_aa > threshold))
        surfing_e = int(np.sum(frame_e > threshold))
        surfing_q = int(np.sum(frame_q > threshold))

        frame_start_bp = float(frame_pos[0]) if len(frame_pos) else np.nan
        frame_end_bp = float(frame_pos[-1]) if len(frame_pos) else np.nan

        row = {
            "Frame #": i + 1,
            "Frame start bp": frame_start_bp,
            "Frame end bp": frame_end_bp,
            f"Surfing_{aa_selected}": surfing_aa,
            "Surfing_E": surfing_e,
            "Surfing_Q": surfing_q,
        }
        for k in range(window):
            row[f"Protein_{k+1}"] = frame_acc[k]
        rows.append(row)

    return pd.DataFrame(rows)


# =========================================================
# Excel export (final)
# =========================================================
def to_excel_bytes(
    df_all: pd.DataFrame,
    df_valid: pd.DataFrame,
    avg_comp_all: pd.Series,
    eq_summary_df: pd.DataFrame,
    walk_df: Optional[pd.DataFrame] = None,
    fig_hist=None,
    fig_corr=None,
    fig_walk=None,
):
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
        df_all.to_excel(writer, index=False, sheet_name="AllProteins")
        df_valid.to_excel(writer, index=False, sheet_name="Valid_EQ")
        avg_comp_all.to_frame("Average (all proteins)").to_excel(writer, sheet_name="AverageAA_ALL")
        eq_summary_df.to_excel(writer, sheet_name="EQ_Summary_VALID")

        if walk_df is not None and not walk_df.empty:
            walk_df.to_excel(writer, index=False, sheet_name="Walking")

        if fig_hist is not None:
            sheet = writer.sheets["EQ_Summary_VALID"]
            imgdata = io.BytesIO()
            fig_hist.savefig(imgdata, format="png", bbox_inches="tight")
            imgdata.seek(0)
            sheet.insert_image("D2", "eq_hist.png", {"image_data": imgdata})

        if fig_corr is not None:
            sheet = writer.sheets["EQ_Summary_VALID"]
            imgdata2 = io.BytesIO()
            fig_corr.savefig(imgdata2, format="png", bbox_inches="tight")
            imgdata2.seek(0)
            sheet.insert_image("D22", "corr.png", {"image_data": imgdata2})

        if fig_walk is not None and walk_df is not None and not walk_df.empty:
            sheetw = writer.sheets.get("Walking")
            if sheetw is None:
                sheetw = writer.book.add_worksheet("Walking")
                writer.sheets["Walking"] = sheetw

            imgdata3 = io.BytesIO()
            fig_walk.savefig(imgdata3, format="png", bbox_inches="tight")
            imgdata3.seek(0)
            sheetw.insert_image("A2", "walk_plot.png", {"image_data": imgdata3})

    return output.getvalue()


# =========================================================
# Utility: preview + downloads (CSV + Excel) without rendering full table
# =========================================================
def preview_and_download(df: pd.DataFrame, preview_cols: list[str], title: str, base_filename: str):
    st.subheader(title)
    preview_n = st.number_input(
        f"Preview rows ({title})",
        min_value=1,
        max_value=100,
        value=5,
        step=1,
        key=f"preview_{title}",
    )
    st.caption("Preview only to keep the app responsive. Download the full table below.")
    st.dataframe(df[preview_cols].head(preview_n), use_container_width=True)

    csv_bytes = df[preview_cols].to_csv(index=False).encode("utf-8")
    st.download_button(
        label=f"⬇️ Download full table (CSV): {title}",
        data=csv_bytes,
        file_name=f"{base_filename}.csv",
        mime="text/csv",
        key=f"csv_{title}",
    )

    out = io.BytesIO()
    with pd.ExcelWriter(out, engine="xlsxwriter") as w:
        df[preview_cols].to_excel(w, index=False, sheet_name="Table")
    st.download_button(
        label=f"⬇️ Download full table (Excel): {title}",
        data=out.getvalue(),
        file_name=f"{base_filename}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key=f"xlsx_{title}",
    )


# =========================================================
# STEP 1
# =========================================================
st.header("Step 1 — Raw Ensembl protein table")
run_step1 = st.button("1) Generate raw Ensembl table")

if run_step1:
    if not chromosome:
        st.error("Please enter a chromosome first.")
        st.stop()

    chrom_ok = chromosome.upper()
    if chrom_ok not in ["X", "Y"]:
        chrom_ok = chromosome
    st.session_state["chrom_ok"] = chrom_ok

    with st.spinner(f"Fetching Ensembl data for chromosome {chrom_ok}..."):
        try:
            ensembl_df = fetch_ensembl_chromosome_data(chrom_ok)
        except Exception as e:
            st.error(f"Ensembl/pybiomart error: {e}")
            st.stop()

    if ensembl_df.empty:
        st.error("❌ No Ensembl rows found for this chromosome (or no Swiss-Prot mapping / peptide available).")
        st.stop()

    st.session_state["ensembl_df"] = ensembl_df

    # Reset downstream results
    for kk in ["df_all", "df_valid", "walk_df", "fig_corr", "fig_hist", "fig_walk", "eq_summary_df", "avg_comp_all"]:
        st.session_state[kk] = None

    st.success(f"✅ Found {len(ensembl_df)} rows (gene/protein) with Swiss-Prot mapping + peptide on chromosome {chrom_ok}.")

if st.session_state["ensembl_df"] is not None:
    ensembl_df = st.session_state["ensembl_df"]
    chrom_ok = st.session_state["chrom_ok"]

    cols_step1 = ["Gene stable ID", "Gene name", "Gene start (bp)", "Gene end (bp)", "Entry", "Sequence"]
    cols_step1 = [c for c in cols_step1 if c in ensembl_df.columns and c != "Chromosome/scaffold name"]

    preview_and_download(
        df=ensembl_df,
        preview_cols=cols_step1,
        title="Raw Ensembl table (Swiss-Prot mapping + peptide)",
        base_filename=f"chromosome_{chrom_ok}_ensembl_raw",
    )


# =========================================================
# STEP 2
# =========================================================
if st.session_state["ensembl_df"] is not None:
    st.header("Step 2 — E/Q analysis (fast AA counting)")
    batch_size = st.slider("Batch size (speed vs responsiveness)", 200, 2000, 800, 100)

    run_step2 = st.button("2) Generate E/Q analysis (fast)")

    if run_step2:
        ensembl_df = st.session_state["ensembl_df"]
        chrom_ok = st.session_state["chrom_ok"]

        with st.spinner("Computing amino acid counts + E/Q ratio (fast)..."):
            df_all = add_aa_counts_and_eq_fast(ensembl_df, batch_size=int(batch_size))

        st.session_state["df_all"] = df_all

        cols_all = [
            "Gene stable ID",
            "Gene name",
            "Gene start (bp)",
            "Gene end (bp)",
            "Entry",
            "Length",
            "E/Q ratio",
            "E/Q valid",
        ] + AA_ORDER
        cols_all = [c for c in cols_all if c in df_all.columns and c != "Chromosome/scaffold name"]

        preview_and_download(
            df=df_all,
            preview_cols=cols_all,
            title="AA counts + E/Q table (full table available via download)",
            base_filename=f"chromosome_{chrom_ok}_AAcounts_EQ",
        )

        invalid_mask = ~df_all["E/Q valid"]
        n_invalid = int(invalid_mask.sum())
        n_total = int(len(df_all))
        n_valid = n_total - n_invalid

        st.subheader("✅/❌ E/Q usability")
        st.write(f"Total proteins found: **{n_total}**")
        st.write(f"Proteins with **non-usable** E/Q (Q = 0): **{n_invalid}**")
        st.write(f"Proteins with **valid** E/Q: **{n_valid}**")

        if n_invalid > 0:
            with st.expander("Show rows with non-usable E/Q", expanded=False):
                show_bad_cols = [c for c in ["Gene stable ID", "Gene name", "Entry", "Length", "E", "Q", "E/Q ratio"] if c in df_all.columns]
                st.dataframe(df_all.loc[invalid_mask, show_bad_cols], use_container_width=True)

        st.subheader("📈 Average amino acid composition (all proteins)")
        avg_comp_all = df_all[AA_ORDER].mean(numeric_only=True)
        st.session_state["avg_comp_all"] = avg_comp_all
        st.bar_chart(avg_comp_all)

        st.subheader("📋 E/Q ratio summary statistics (valid E/Q only)")
        df_valid = df_all[df_all["E/Q valid"]].copy().reset_index(drop=True)
        st.session_state["df_valid"] = df_valid

        if df_valid.empty:
            st.error("No proteins with valid E/Q on this chromosome.")
            st.stop()

        eq_summary = {
            "Mean": df_valid["E/Q ratio"].mean(),
            "Median": df_valid["E/Q ratio"].median(),
            "Standard deviation": df_valid["E/Q ratio"].std(),
            "Minimum": df_valid["E/Q ratio"].min(),
            "Maximum": df_valid["E/Q ratio"].max(),
            "Count (valid EQ)": len(df_valid),
            "Total proteins": n_total,
            "Invalid EQ (Q=0)": n_invalid,
        }
        eq_summary_df = pd.DataFrame(eq_summary, index=["Value"]).T
        st.session_state["eq_summary_df"] = eq_summary_df
        st.table(eq_summary_df)

        st.subheader("🔗 Correlation between protein length and E/Q ratio (valid proteins)")
        corr_df = df_valid.dropna(subset=["Length", "E/Q ratio"])
        fig_corr = None
        if len(corr_df) > 2:
            r, p = pearsonr(corr_df["Length"], corr_df["E/Q ratio"])
            st.write(f"**Pearson r** = `{r:.3f}` (p = `{p:.3e}`)")
            fig_corr, ax = plt.subplots(figsize=(8, 5))
            sns.regplot(data=corr_df, x="Length", y="E/Q ratio", scatter_kws={"alpha": 0.6})
            ax.set_title(f"Length vs E/Q (r={r:.2f})")
            st.pyplot(fig_corr)
        else:
            st.warning("Not enough valid data points to compute correlation.")
        st.session_state["fig_corr"] = fig_corr

        st.subheader("📊 E/Q ratio distribution (valid proteins)")
        fig_hist, axh = plt.subplots(figsize=(10, 5))
        sns.histplot(df_valid["E/Q ratio"], bins=30, kde=True, ax=axh)
        axh.set_title(f"E/Q ratio distribution (Chromosome {chrom_ok})")
        axh.set_xlabel("E/Q ratio")
        axh.set_ylabel("Protein count")
        st.pyplot(fig_hist)
        st.session_state["fig_hist"] = fig_hist

        st.success("✅ Step 2 completed.")


# =========================================================
# STEP 3
# =========================================================
if st.session_state["df_all"] is not None:
    st.header("Step 3 — Walking (Surfing) analysis + chromosome-wide plot")

    df_all = st.session_state["df_all"]
    chrom_ok = st.session_state["chrom_ok"]

    st.write("""
Walking definition (implemented):
A sliding window of fixed size (number of proteins) is moved along the chromosome-ordered protein list (sorted by Gene start bp).
For each window, we count how many proteins have an absolute count > T for:
- selected amino acid (AA)
- E
- Q
""")

    aa_selected = st.selectbox("Amino acid for walking (AA):", options=AA_ORDER, index=AA_ORDER.index("E"))
    window = st.slider("Frame size (window)", min_value=5, max_value=200, value=25, step=1)
    threshold = st.number_input("Threshold T (AA count > T)", min_value=0.0, value=49.0, step=1.0)

    x_axis_choice = st.radio("Walking plot X-axis:", options=["Frame start bp", "Frame #"], horizontal=True)
    smooth = st.checkbox("Smooth curves (rolling mean)", value=True)
    smooth_window = st.slider("Smoothing window (frames)", 1, 200, 25, 1) if smooth else 1

    run_step3 = st.button("3) Generate walking analysis")

    if run_step3:
        df_for_walk = df_all.sort_values(by=["Gene start (bp)"], na_position="last").reset_index(drop=True)

        try:
            walk_df = make_walk_table(
                df_all_sorted=df_for_walk,
                window=window,
                threshold=threshold,
                aa_selected=aa_selected,
                accession_col="Entry",
                pos_col="Gene start (bp)",
            )
        except Exception as e:
            st.error(f"Walking analysis error: {e}")
            st.stop()

        if walk_df.empty:
            st.warning("Not enough proteins to compute the walking table with the selected frame size.")
            st.stop()

        st.session_state["walk_df"] = walk_df

        st.caption(f"Walking table rows: {len(walk_df)} (computed as N - window + 1)")
        st.dataframe(walk_df.head(50), use_container_width=True)
        st.caption("Showing first 50 rows only. Download the full walking table below.")

        csv_bytes = walk_df.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="⬇️ Download full Walking table (CSV)",
            data=csv_bytes,
            file_name=f"chromosome_{chrom_ok}_walking_{aa_selected}_w{window}_T{int(threshold)}.csv",
            mime="text/csv",
        )

        st.subheader("📈 Walking quick plot (chromosome overview)")
        x_col = "Frame start bp" if x_axis_choice == "Frame start bp" else "Frame #"
        y_cols = [f"Surfing_{aa_selected}", "Surfing_E", "Surfing_Q"]

        plot_df = walk_df[[x_col] + y_cols].copy().sort_values(by=x_col)

        if smooth and smooth_window > 1:
            for yc in y_cols:
                plot_df[yc] = plot_df[yc].rolling(window=smooth_window, min_periods=1).mean()

        fig_walk, axw = plt.subplots(figsize=(12, 5))
        for yc in y_cols:
            axw.plot(plot_df[x_col], plot_df[yc], label=yc)
        axw.set_title(f"Walking (window={window}, T={threshold}) — Chromosome {chrom_ok}")
        axw.set_xlabel(x_col)
        axw.set_ylabel("Count of proteins in frame with (AA count > T)")
        axw.legend()
        st.pyplot(fig_walk)

        st.session_state["fig_walk"] = fig_walk
        st.success("✅ Step 3 completed.")


# =========================================================
# FINAL Excel export (generated only when requested)
# =========================================================
if (
    st.session_state["df_all"] is not None
    and st.session_state["df_valid"] is not None
    and st.session_state["avg_comp_all"] is not None
    and st.session_state["eq_summary_df"] is not None
):
    st.header("Final export")

    if st.button("Generate final Excel file (All + Valid EQ + AverageAA + Stats + Walking + Plots)"):
        with st.spinner("Creating Excel file..."):
            df_all = st.session_state["df_all"]
            df_valid = st.session_state["df_valid"]
            avg_comp_all = st.session_state["avg_comp_all"]
            eq_summary_df = st.session_state["eq_summary_df"]

            walk_df = st.session_state.get("walk_df")
            fig_hist = st.session_state.get("fig_hist")
            fig_corr = st.session_state.get("fig_corr")
            fig_walk = st.session_state.get("fig_walk")

            excel_bytes = to_excel_bytes(
                df_all=df_all,
                df_valid=df_valid,
                avg_comp_all=avg_comp_all,
                eq_summary_df=eq_summary_df,
                walk_df=walk_df if isinstance(walk_df, pd.DataFrame) else None,
                fig_hist=fig_hist,
                fig_corr=fig_corr,
                fig_walk=fig_walk,
            )

        st.download_button(
            label="⬇️ Download Excel",
            data=excel_bytes,
            file_name=f"chromosome_{st.session_state['chrom_ok']}_EQ_summary.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )
