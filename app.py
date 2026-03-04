import streamlit as st
import pandas as pd
import numpy as np
import io
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.stats import pearsonr
from pybiomart import Dataset

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

Workflow:
1) Show **all** proteins found (Swiss-Prot mapping + peptide)  
2) Compute **E/Q ratio** and mark proteins where E/Q is **not usable** (typically **Q = 0**)  
3) Compute **average amino acid composition** using **all** proteins found  
4) Compute **E/Q summary statistics** using **only proteins with valid E/Q**  
5) Correlation + distribution on valid E/Q proteins  
6) **Walking (Surfing) analysis** on the chromosome-ordered protein list (AA selected + E + Q only) + quick chromosome-wide plot
""")

chromosome = st.text_input("Chromosome (1–22, X, Y):", "").strip()
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


# =========================================================
# Data fetch (Ensembl / BioMart)
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
            # Friendly labels (if returned)
            "Gene stable ID": "Gene stable ID",
            "Gene name": "Gene name",
            "Gene start (bp)": "Gene start (bp)",
            "Gene end (bp)": "Gene end (bp)",
            "Chromosome/scaffold name": "Chromosome/scaffold name",
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
# Core analysis: AA counts + E/Q validity
# =========================================================
def add_aa_counts_and_eq(df: pd.DataFrame) -> pd.DataFrame:
    aa_counts_list = []
    eq_list = []
    eq_valid_list = []
    length_list = []

    for seq in df["Sequence"].astype(str):
        seq = seq.rstrip("*")
        analysis = ProteinAnalysis(seq)
        counts = analysis.count_amino_acids()

        e = counts.get("E", 0)
        q = counts.get("Q", 0)
        eq = e / q if q != 0 else np.nan
        eq_valid = (q != 0)

        aa_counts_list.append({aa: int(counts.get(aa, 0)) for aa in AA_ORDER})
        eq_list.append(eq)
        eq_valid_list.append(eq_valid)
        length_list.append(len(seq))

    aa_df = pd.DataFrame(aa_counts_list)
    out = pd.concat([df.reset_index(drop=True), aa_df], axis=1)
    out["Length"] = length_list
    out["E/Q ratio"] = eq_list
    out["E/Q valid"] = eq_valid_list
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

        # Frame start/end bp (current definition: start positions of first and last protein in the frame)
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
# Excel export
# =========================================================
def to_excel_bytes(
    df_all: pd.DataFrame,
    df_valid: pd.DataFrame,
    avg_comp_all: pd.Series,
    eq_summary_df: pd.DataFrame,
    walk_df: pd.DataFrame | None = None,
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
            # Put walk plot into Walking sheet if exists, else create it
            if "Walking" in writer.sheets:
                sheetw = writer.sheets["Walking"]
            else:
                sheetw = writer.book.add_worksheet("Walking")
                writer.sheets["Walking"] = sheetw

            imgdata3 = io.BytesIO()
            fig_walk.savefig(imgdata3, format="png", bbox_inches="tight")
            imgdata3.seek(0)
            sheetw.insert_image("A2", "walk_plot.png", {"image_data": imgdata3})

    return output.getvalue()


# =========================================================
# Run
# =========================================================
if chromosome:
    chrom_ok = chromosome.upper()
    if chrom_ok not in ["X", "Y"]:
        chrom_ok = chromosome

    with st.spinner(f"Fetching Ensembl data for chromosome {chrom_ok}..."):
        try:
            ensembl_df = fetch_ensembl_chromosome_data(chrom_ok)
        except Exception as e:
            st.error(f"Ensembl/pybiomart error: {e}")
            st.stop()

    if ensembl_df.empty:
        st.error("❌ No Ensembl rows found for this chromosome (or no Swiss-Prot mapping / peptide available).")
        st.stop()

    st.success(f"✅ Found {len(ensembl_df)} rows (gene/protein) with Swiss-Prot mapping + peptide on chromosome {chrom_ok}.")

    with st.spinner("Computing amino acid counts + E/Q ratio..."):
        df_all = add_aa_counts_and_eq(ensembl_df)

    # ---------------------------------------------------------
    # Table 1: ALL proteins (no Chromosome/scaffold column)
    # ---------------------------------------------------------
    st.subheader("📊 Full table (all proteins found)")
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
    st.dataframe(df_all[cols_all])

    # ---------------------------------------------------------
    # E/Q usability report
    # ---------------------------------------------------------
    invalid_mask = ~df_all["E/Q valid"]
    n_invalid = int(invalid_mask.sum())
    n_total = int(len(df_all))
    n_valid = n_total - n_invalid

    st.subheader("✅/❌ E/Q usability")
    st.write(f"Total proteins found: **{n_total}**")
    st.write(f"Proteins with **non-usable** E/Q (Q = 0): **{n_invalid}**")
    st.write(f"Proteins with **valid** E/Q: **{n_valid}**")

    if n_invalid > 0:
        with st.expander("Show rows with non-usable E/Q"):
            show_bad_cols = [c for c in ["Gene stable ID", "Gene name", "Entry", "Length", "E", "Q", "E/Q ratio"] if c in df_all.columns]
            st.dataframe(df_all.loc[invalid_mask, show_bad_cols])

    # ---------------------------------------------------------
    # Average AA composition on ALL proteins
    # ---------------------------------------------------------
    st.subheader("📈 Average amino acid composition (all proteins)")
    avg_comp_all = df_all[AA_ORDER].mean(numeric_only=True)
    st.bar_chart(avg_comp_all)

    # ---------------------------------------------------------
    # E/Q summary on VALID proteins only
    # ---------------------------------------------------------
    st.subheader("📋 E/Q ratio summary statistics (valid E/Q only)")
    df_valid = df_all[df_all["E/Q valid"]].copy().reset_index(drop=True)

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
    st.table(eq_summary_df)

    # ---------------------------------------------------------
    # Correlation (valid proteins)
    # ---------------------------------------------------------
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

    # ---------------------------------------------------------
    # Histogram (valid proteins)
    # ---------------------------------------------------------
    st.subheader("📊 E/Q ratio distribution (valid proteins)")
    fig_hist, axh = plt.subplots(figsize=(10, 5))
    sns.histplot(df_valid["E/Q ratio"], bins=30, kde=True, ax=axh)
    axh.set_title(f"E/Q ratio distribution (Chromosome {chrom_ok})")
    axh.set_xlabel("E/Q ratio")
    axh.set_ylabel("Protein count")
    st.pyplot(fig_hist)

    # ---------------------------------------------------------
    # Walking (Surfing) analysis + quick chromosome plot
    # ---------------------------------------------------------
    st.subheader("🏄 Walking (Surfing) analysis along chromosome-ordered proteins")

    st.write("""
**Walking definition (implemented):**
Given a selected frame size (window), the software slides a block along the protein list ordered by chromosome position.
For each frame, it counts how many proteins have an **absolute count** above a fixed threshold **T** for:
- the selected amino acid (AA),
- **E**,
- **Q**.

The results are stored in **Surfing_*** columns.
The table also includes:
- **Frame start bp** and **Frame end bp**,
- the **protein accessions** in the frame (**Protein_1 ... Protein_window**).
""")

    aa_selected = st.selectbox("Amino acid for walking (AA):", options=AA_ORDER, index=AA_ORDER.index("E"))
    window = st.slider("Frame size (window)", min_value=5, max_value=200, value=25, step=1)
    threshold = st.number_input("Threshold T (AA count > T)", min_value=0.0, value=49.0, step=1.0)

    # Plot controls
    x_axis_choice = st.radio("Walking plot X-axis:", options=["Frame start bp", "Frame #"], horizontal=True)
    smooth = st.checkbox("Smooth curves (rolling mean)", value=True)
    smooth_window = st.slider("Smoothing window (frames)", 1, 200, 25, 1) if smooth else 1

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
        walk_df = pd.DataFrame()

    fig_walk = None
    if walk_df.empty:
        st.warning("Not enough proteins to compute the walking table with the selected frame size.")
    else:
        st.caption(f"Walking table rows: {len(walk_df)} (computed as N - window + 1)")
        st.dataframe(walk_df)

        # -------- Quick chromosome-wide plot --------
        st.subheader("📈 Walking quick plot (chromosome overview)")

        x_col = "Frame start bp" if x_axis_choice == "Frame start bp" else "Frame #"
        y_cols = [f"Surfing_{aa_selected}", "Surfing_E", "Surfing_Q"]

        plot_df = walk_df[[x_col] + y_cols].copy()
        plot_df = plot_df.sort_values(by=x_col)

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

        # Downloads
        csv_bytes = walk_df.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="Download Walking table (CSV)",
            data=csv_bytes,
            file_name=f"chromosome_{chrom_ok}_walking_{aa_selected}_w{window}_T{int(threshold)}.csv",
            mime="text/csv",
        )

    # ---------------------------------------------------------
    # Excel export
    # ---------------------------------------------------------
    st.subheader("⬇️ Export Excel")
    excel_bytes = to_excel_bytes(
        df_all=df_all,
        df_valid=df_valid,
        avg_comp_all=avg_comp_all,
        eq_summary_df=eq_summary_df,
        walk_df=walk_df if not walk_df.empty else None,
        fig_hist=fig_hist,
        fig_corr=fig_corr,
        fig_walk=fig_walk,
    )

    st.download_button(
        label="Download Excel (All + Valid EQ + AverageAA + Stats + Walking + Plots)",
        data=excel_bytes,
        file_name=f"chromosome_{chrom_ok}_EQ_summary.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

    st.success("✅ Analysis completed.")
