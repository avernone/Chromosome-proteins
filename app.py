import streamlit as st
import pandas as pd
import numpy as np
import io
import hashlib
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

Workflow:
- **Step 1 (gene-level):** gene table with **Swiss-Prot AC list** per gene (ACs separated by `;`)
- **Step 2 (walking):** walking/surfing on the **protein list** ordered by chromosome position
- **Step 3 (E/Q last):** full AA composition + E/Q ratio statistics and plots
""")

chromosome = st.text_input("Chromosome (1–22, X, Y):", "").strip()
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

# =========================================================
# Session state init
# =========================================================
for k in [
    "chrom_ok",
    "df_gene",
    "df_protein_raw",
    "walk_df",
    "eq_window_df",
    "df_all",
    "df_valid",
    "eq_summary_df",
    "avg_comp_all",
    "fig_corr",
    "fig_hist",
    "fig_walk",
    "fig_eq_window",
]:
    if k not in st.session_state:
        st.session_state[k] = None

# =========================================================
# Amino Acid Reference Table
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
    st.dataframe(aa_df_ref, hide_index=True, width="stretch")

# =========================================================
# Helpers
# =========================================================
def normalize_biomart_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {
        "ensembl_gene_id": "Gene stable ID",
        "external_gene_name": "Gene name",
        "start_position": "Gene start (bp)",
        "end_position": "Gene end (bp)",
        "chromosome_name": "Chromosome/scaffold name",
        "uniprotswissprot": "Swiss-Prot ACs",
        "peptide": "Sequence",
        "UniProtKB/Swiss-Prot ID": "Swiss-Prot ACs",
        "Peptide": "Sequence",
    }
    cols = {c: rename_map[c] for c in df.columns if c in rename_map}
    return df.rename(columns=cols)


def sha1_10(s: str) -> str:
    return hashlib.sha1(s.encode("utf-8")).hexdigest()[:10]


# =========================================================
# Step 1 fetch (gene-level)
# =========================================================
@st.cache_data(show_spinner=False, ttl=24 * 3600)
def fetch_ensembl_gene_level(chrom: str) -> pd.DataFrame:
    dataset = Dataset(name="hsapiens_gene_ensembl", host="http://www.ensembl.org")
    df = dataset.query(
        attributes=[
            "ensembl_gene_id",
            "external_gene_name",
            "start_position",
            "end_position",
            "chromosome_name",
            "uniprotswissprot",
        ],
        filters={"chromosome_name": [chrom]},
    )
    df = normalize_biomart_columns(df)

    df["Swiss-Prot ACs"] = df["Swiss-Prot ACs"].replace("", np.nan)
    df = df.dropna(subset=["Swiss-Prot ACs"]).copy()

    df["Gene start (bp)"] = pd.to_numeric(df["Gene start (bp)"], errors="coerce")
    df["Gene end (bp)"] = pd.to_numeric(df["Gene end (bp)"], errors="coerce")

    df = df.drop_duplicates(
        subset=["Gene stable ID", "Gene start (bp)", "Gene end (bp)", "Swiss-Prot ACs"],
        keep="first",
    )

    df["Isoform detected"] = df["Swiss-Prot ACs"].astype(str).str.contains("-", regex=False)

    return df.sort_values(by=["Gene start (bp)"], na_position="last").reset_index(drop=True)


# =========================================================
# Protein-level raw fetch (with peptides)
# =========================================================
@st.cache_data(show_spinner=False, ttl=24 * 3600)
def fetch_ensembl_protein_level_raw(chrom: str) -> pd.DataFrame:
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
    df = normalize_biomart_columns(df)

    df["Swiss-Prot ACs"] = df["Swiss-Prot ACs"].replace("", np.nan)
    df["Sequence"] = df["Sequence"].replace("", np.nan)
    df = df.dropna(subset=["Swiss-Prot ACs", "Sequence"]).copy()

    df["Sequence"] = df["Sequence"].astype(str).str.rstrip("*")
    df["Gene start (bp)"] = pd.to_numeric(df["Gene start (bp)"], errors="coerce")
    df["Gene end (bp)"] = pd.to_numeric(df["Gene end (bp)"], errors="coerce")

    df["Entry"] = df["Swiss-Prot ACs"].astype(str).str.split(";")
    df = df.explode("Entry")
    df["Entry"] = df["Entry"].astype(str).str.strip()
    df = df[df["Entry"].ne("")].copy()

    df["SeqLen"] = df["Sequence"].astype(str).str.len()
    df["SeqHash"] = df["Sequence"].astype(str).apply(lambda s: sha1_10(s))

    df["Entry_has_multiple_sequences"] = df.groupby("Entry")["SeqHash"].transform("nunique") > 1
    df["Isoform detected"] = (
        df["Entry"].astype(str).str.contains("-", regex=False)
        | df["Entry_has_multiple_sequences"]
    )

    df = df.drop_duplicates(
        subset=["Entry", "Sequence", "Gene start (bp)", "Gene end (bp)"],
        keep="first",
    )

    return df.sort_values(by=["Gene start (bp)"], na_position="last").reset_index(drop=True)


# =========================================================
# FAST AA counting
# =========================================================
def fast_count_aa(seq: str) -> dict:
    seq = str(seq).rstrip("*").upper()
    c = Counter(seq)
    return {aa: int(c.get(aa, 0)) for aa in AA_ORDER}


def add_aa_counts_and_eq_fast(df: pd.DataFrame, batch_size: int = 1200) -> pd.DataFrame:
    n = len(df)
    progress = st.progress(0)
    status = st.empty()

    out_chunks = []
    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        batch = df.iloc[start:end]

        rows = []
        for seq in batch["Sequence"].astype(str):
            s = str(seq).rstrip("*")
            counts = fast_count_aa(s)
            e = counts["E"]
            q = counts["Q"]
            eq = (e / q) if q != 0 else np.nan
            rows.append(
                {
                    **counts,
                    "Length": len(s),
                    "E/Q ratio": eq,
                    "E/Q valid": (q != 0),
                }
            )

        out_chunks.append(pd.DataFrame(rows))
        progress.progress(end / n)
        status.write(f"Computed {end}/{n} proteins...")

    progress.empty()
    status.empty()

    aa_df = pd.concat(out_chunks, ignore_index=True)
    out = pd.concat([df.reset_index(drop=True), aa_df], axis=1)
    return out


# =========================================================
# Walking / Surfing
# =========================================================
def make_walk_table_fast_with_list(
    df_sorted: pd.DataFrame,
    window: int,
    threshold: float,
    aa_selected: str,
    accession_col: str = "Entry",
    pos_col: str = "Gene start (bp)",
) -> pd.DataFrame:
    n = len(df_sorted)
    if n < window:
        return pd.DataFrame()

    seqs = df_sorted["Sequence"].astype(str).tolist()
    pos = pd.to_numeric(df_sorted[pos_col], errors="coerce").to_numpy()
    acc = df_sorted[accession_col].astype(str).to_numpy()

    def count_letter(s: str, letter: str) -> int:
        return str(s).count(letter)

    aa_counts = np.array([count_letter(s, aa_selected) for s in seqs], dtype=np.int32)
    e_counts = np.array([count_letter(s, "E") for s in seqs], dtype=np.int32)
    q_counts = np.array([count_letter(s, "Q") for s in seqs], dtype=np.int32)

    b_aa = (aa_counts > threshold).astype(np.int32)
    b_e = (e_counts > threshold).astype(np.int32)
    b_q = (q_counts > threshold).astype(np.int32)

    def window_sum_bool(b: np.ndarray, w: int) -> np.ndarray:
        cs = np.concatenate([[0], np.cumsum(b)])
        return cs[w:] - cs[:-w]

    surf_aa = window_sum_bool(b_aa, window)
    surf_e = window_sum_bool(b_e, window)
    surf_q = window_sum_bool(b_q, window)

    frame_start_bp = pos[: n - window + 1]
    frame_end_bp = pos[window - 1 :]

    prot_list = [";".join(acc[i : i + window]) for i in range(0, n - window + 1)]

    out = pd.DataFrame(
        {
            "Frame #": np.arange(1, n - window + 2, dtype=np.int32),
            "Frame start bp": frame_start_bp,
            "Frame end bp": frame_end_bp,
            f"Surfing_{aa_selected}": surf_aa,
            "Surfing_E": surf_e,
            "Surfing_Q": surf_q,
            "Proteins in window": prot_list,
        }
    )
    return out


# =========================================================
# E/Q sliding-window table
# =========================================================
def make_eq_window_table(
    df_sorted: pd.DataFrame,
    window: int,
    accession_col: str = "Entry",
    pos_col: str = "Gene start (bp)",
) -> pd.DataFrame:
    n = len(df_sorted)
    if n < window:
        return pd.DataFrame()

    work = df_sorted.copy()

    work["E_count"] = work["Sequence"].astype(str).str.count("E")
    work["Q_count"] = work["Sequence"].astype(str).str.count("Q")
    work["E/Q ratio"] = work["E_count"] / work["Q_count"]
    work.loc[work["Q_count"] == 0, "E/Q ratio"] = np.nan

    pos = pd.to_numeric(work[pos_col], errors="coerce").to_numpy()
    acc = work[accession_col].astype(str).to_numpy()
    eq_vals = pd.to_numeric(work["E/Q ratio"], errors="coerce").to_numpy()

    rows = []
    for i in range(0, n - window + 1):
        j = i + window

        frame_eq = eq_vals[i:j]
        frame_pos = pos[i:j]
        frame_acc = acc[i:j]

        row = {
            "Frame #": i + 1,
            "Frame start bp": float(frame_pos[0]) if len(frame_pos) else np.nan,
            "Frame end bp": float(frame_pos[-1]) if len(frame_pos) else np.nan,
            "Mean E/Q ratio": np.nanmean(frame_eq),
            "Valid E/Q proteins in window": int(np.sum(~np.isnan(frame_eq))),
            "Proteins in window": ";".join(frame_acc),
        }
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
    walk_df: Optional[pd.DataFrame] = None,
    eq_window_df: Optional[pd.DataFrame] = None,
    fig_hist=None,
    fig_corr=None,
    fig_walk=None,
    fig_eq_window=None,
):
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
        df_all.to_excel(writer, index=False, sheet_name="AllProteins")
        df_valid.to_excel(writer, index=False, sheet_name="Valid_EQ")
        avg_comp_all.to_frame("Average (all proteins)").to_excel(writer, sheet_name="AverageAA_ALL")
        eq_summary_df.to_excel(writer, sheet_name="EQ_Summary_VALID")

        if walk_df is not None and not walk_df.empty:
            walk_df.to_excel(writer, index=False, sheet_name="Walking")

        if eq_window_df is not None and not eq_window_df.empty:
            eq_window_df.to_excel(writer, index=False, sheet_name="EQ_Window")

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

        if fig_eq_window is not None and eq_window_df is not None and not eq_window_df.empty:
            sheeteq = writer.sheets.get("EQ_Window")
            if sheeteq is None:
                sheeteq = writer.book.add_worksheet("EQ_Window")
                writer.sheets["EQ_Window"] = sheeteq
            imgdata4 = io.BytesIO()
            fig_eq_window.savefig(imgdata4, format="png", bbox_inches="tight")
            imgdata4.seek(0)
            sheeteq.insert_image("A2", "eq_window_plot.png", {"image_data": imgdata4})

    return output.getvalue()


# =========================================================
# Utility: preview + downloads
# =========================================================
def preview_and_download(
    df: pd.DataFrame,
    preview_cols: list[str],
    title: str,
    base_filename: str,
    preview_default: int = 5,
):
    st.subheader(title)
    preview_n = st.number_input(
        f"Preview rows ({title})",
        min_value=1,
        max_value=200,
        value=preview_default,
        step=1,
        key=f"preview_{base_filename}",
    )
    st.caption("Preview only to keep the app responsive. Download the full table below.")
    st.dataframe(df[preview_cols].head(int(preview_n)), width="stretch")

    csv_bytes = df[preview_cols].to_csv(index=False).encode("utf-8")
    st.download_button(
        label=f"⬇️ Download full table (CSV): {title}",
        data=csv_bytes,
        file_name=f"{base_filename}.csv",
        mime="text/csv",
        key=f"csv_{base_filename}",
    )

    out = io.BytesIO()
    with pd.ExcelWriter(out, engine="xlsxwriter") as w:
        df[preview_cols].to_excel(w, index=False, sheet_name="Table")
    st.download_button(
        label=f"⬇️ Download full table (Excel): {title}",
        data=out.getvalue(),
        file_name=f"{base_filename}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key=f"xlsx_{base_filename}",
    )


# =========================================================
# STEP 1 — Gene-level table
# =========================================================
st.header("Step 1 — Gene-level table (show everything)")
st.caption("Gene-centric table (NOT protein-centric). **Swiss-Prot ACs** may contain multiple accessions separated by ';'.")

run_step1 = st.button("1) Generate gene-level table")

if run_step1:
    if not chromosome:
        st.error("Please enter a chromosome first.")
        st.stop()

    chrom_ok = chromosome.upper()
    if chrom_ok not in ["X", "Y"]:
        chrom_ok = chromosome
    st.session_state["chrom_ok"] = chrom_ok

    with st.spinner(f"Fetching Ensembl gene-level data for chromosome {chrom_ok}..."):
        try:
            df_gene = fetch_ensembl_gene_level(chrom_ok)
        except Exception as e:
            st.error(f"Ensembl/pybiomart error (Step 1): {e}")
            st.stop()

    if df_gene.empty:
        st.error("❌ No Ensembl rows found for this chromosome (or no Swiss-Prot mapping available).")
        st.stop()

    st.session_state["df_gene"] = df_gene

    for kk in [
        "df_protein_raw",
        "walk_df",
        "eq_window_df",
        "df_all",
        "df_valid",
        "eq_summary_df",
        "avg_comp_all",
        "fig_corr",
        "fig_hist",
        "fig_walk",
        "fig_eq_window",
    ]:
        st.session_state[kk] = None

    st.success(f"✅ Found {len(df_gene)} gene-level rows with Swiss-Prot mapping on chromosome {chrom_ok}.")

if st.session_state["df_gene"] is not None:
    df_gene = st.session_state["df_gene"]
    chrom_ok = st.session_state["chrom_ok"]

    cols_step1 = ["Gene stable ID", "Gene name", "Gene start (bp)", "Gene end (bp)", "Swiss-Prot ACs", "Isoform detected"]
    cols_step1 = [c for c in cols_step1 if c in df_gene.columns and c != "Chromosome/scaffold name"]

    preview_and_download(
        df=df_gene,
        preview_cols=cols_step1,
        title="Table 1 — Gene-level mapping (Swiss-Prot AC list per gene)",
        base_filename=f"chromosome_{chrom_ok}_table1_gene_level",
        preview_default=10,
    )


# =========================================================
# STEP 2 — Walking
# =========================================================
if st.session_state["df_gene"] is not None:
    st.header("Step 2 — Walking (Surfing) analysis")

    st.checkbox(
        "Use Swiss-Prot reviewed only (recommended)",
        value=True,
        disabled=True,
        help="BioMart attribute 'uniprotswissprot' already returns Swiss-Prot (reviewed) accessions.",
    )

    aa_selected = st.selectbox("Amino acid for walking (AA):", options=AA_ORDER, index=AA_ORDER.index("D"))
    window = st.slider("Frame size (window)", min_value=5, max_value=300, value=25, step=1)
    threshold = st.number_input("Threshold T (AA count > T)", min_value=0.0, value=49.0, step=1.0)

    x_axis_choice = st.radio("Walking plot X-axis:", options=["Frame start bp", "Frame #"], horizontal=True)
    smooth = st.checkbox("Smooth curves (rolling mean)", value=True)
    smooth_window = st.slider("Smoothing window (frames)", 1, 300, 25, 1) if smooth else 1

    run_step2 = st.button("2) Generate walking analysis")

    if run_step2:
        chrom_ok = st.session_state["chrom_ok"]

        with st.spinner(f"Fetching Ensembl protein-level data (with peptides) for chromosome {chrom_ok}..."):
            try:
                df_prot_raw = fetch_ensembl_protein_level_raw(chrom_ok)
            except Exception as e:
                st.error(f"Ensembl/pybiomart error (Step 2 fetch): {e}")
                st.stop()

        if df_prot_raw.empty:
            st.error("❌ No protein-level rows found (or no peptide available).")
            st.stop()

        df_prot_raw["Is_reviewed"] = True
        st.session_state["df_protein_raw"] = df_prot_raw

        cols_raw = [
            "Gene stable ID", "Gene name", "Gene start (bp)", "Gene end (bp)",
            "Entry", "Is_reviewed", "SeqLen", "SeqHash", "Entry_has_multiple_sequences", "Isoform detected"
        ]
        preview_and_download(
            df=df_prot_raw.drop(columns=["Sequence"], errors="ignore"),
            preview_cols=[c for c in cols_raw if c in df_prot_raw.columns],
            title="Protein-level raw table (used for walking & E/Q) — Sequence not shown in preview",
            base_filename=f"chromosome_{chrom_ok}_protein_raw",
            preview_default=10,
        )

        df_for_walk = df_prot_raw.sort_values(by=["Gene start (bp)"], na_position="last").reset_index(drop=True)

        with st.spinner("Computing Walking table (fast)..."):
            walk_df = make_walk_table_fast_with_list(
                df_sorted=df_for_walk,
                window=int(window),
                threshold=float(threshold),
                aa_selected=str(aa_selected),
            )

        if walk_df.empty:
            st.warning("Not enough proteins to compute the walking table with the selected frame size.")
            st.stop()

        st.session_state["walk_df"] = walk_df

        st.subheader("Walking table preview")
        st.caption(f"Walking table rows: {len(walk_df)} (computed as N - window + 1)")
        st.dataframe(walk_df.head(50), width="stretch")
        st.caption("Showing first 50 rows only. Download the full walking table below.")

        csv_bytes = walk_df.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="⬇️ Download full Walking table (CSV)",
            data=csv_bytes,
            file_name=f"chromosome_{chrom_ok}_walking_{aa_selected}_w{int(window)}_T{int(threshold)}.csv",
            mime="text/csv",
        )

        st.subheader("📈 Walking quick plot (chromosome overview)")
        x_col = "Frame start bp" if x_axis_choice == "Frame start bp" else "Frame #"
        y_cols = [f"Surfing_{aa_selected}", "Surfing_E", "Surfing_Q"]

        plot_df = walk_df[[x_col] + y_cols].copy().sort_values(by=x_col)

        if smooth and smooth_window > 1:
            for yc in y_cols:
                plot_df[yc] = plot_df[yc].rolling(window=int(smooth_window), min_periods=1).mean()

        fig_walk, axw = plt.subplots(figsize=(12, 5))
        for yc in y_cols:
            axw.plot(plot_df[x_col], plot_df[yc], label=yc)
        axw.set_title(f"Walking (window={int(window)}, T={float(threshold)}) — Chromosome {chrom_ok}")
        axw.set_xlabel(x_col)
        axw.set_ylabel("Count of proteins in frame with (AA count > T)")
        axw.legend()
        st.pyplot(fig_walk)

        st.session_state["fig_walk"] = fig_walk

        # ---------------------------------------------------------
        # E/Q sliding-window table + plot
        # ---------------------------------------------------------
        st.subheader("📊 E/Q sliding-window table")

        eq_window_df = make_eq_window_table(
            df_sorted=df_for_walk,
            window=int(window),
            accession_col="Entry",
            pos_col="Gene start (bp)",
        )

        if eq_window_df.empty:
            st.warning("Not enough proteins to compute the E/Q sliding-window table.")
        else:
            st.session_state["eq_window_df"] = eq_window_df

            st.caption(f"E/Q window table rows: {len(eq_window_df)} (computed as N - window + 1)")
            st.dataframe(eq_window_df.head(50), width="stretch")
            st.caption("Showing first 50 rows only. Download the full E/Q window table below.")

            eq_csv_bytes = eq_window_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="⬇️ Download full E/Q window table (CSV)",
                data=eq_csv_bytes,
                file_name=f"chromosome_{chrom_ok}_EQ_window_w{int(window)}.csv",
                mime="text/csv",
            )

            st.subheader("📈 E/Q sliding-window plot (chromosome overview)")

            eq_plot_df = eq_window_df.copy()
            eq_x_col = "Frame start bp" if x_axis_choice == "Frame start bp" else "Frame #"

            if smooth and smooth_window > 1:
                eq_plot_df["Mean E/Q smoothed"] = eq_plot_df["Mean E/Q ratio"].rolling(
                    window=int(smooth_window),
                    min_periods=1
                ).mean()
                y_col = "Mean E/Q smoothed"
            else:
                y_col = "Mean E/Q ratio"

            fig_eq_window, axeq = plt.subplots(figsize=(12, 5))
            axeq.plot(eq_plot_df[eq_x_col], eq_plot_df[y_col])
            axeq.set_title(f"E/Q sliding-window mean along chromosome {chrom_ok}")
            axeq.set_xlabel(eq_x_col)
            axeq.set_ylabel("Mean E/Q ratio")
            st.pyplot(fig_eq_window)

            st.session_state["fig_eq_window"] = fig_eq_window

        st.success("✅ Step 2 (walking) completed.")


# =========================================================
# STEP 3 — E/Q last
# =========================================================
if st.session_state["df_protein_raw"] is not None:
    st.header("Step 3 — E/Q analysis (last step)")

    batch_size = st.slider("Batch size (speed vs responsiveness)", 200, 3000, 1200, 100)
    run_step3 = st.button("3) Generate E/Q analysis")

    if run_step3:
        chrom_ok = st.session_state["chrom_ok"]
        df_prot_raw = st.session_state["df_protein_raw"].copy()

        with st.spinner("Computing full AA counts + E/Q ratio (fast)..."):
            df_all = add_aa_counts_and_eq_fast(df_prot_raw, batch_size=int(batch_size))

        df_all = df_all.drop(columns=["Sequence"], errors="ignore")
        st.session_state["df_all"] = df_all

        cols_all = [
            "Gene stable ID", "Gene name", "Gene start (bp)", "Gene end (bp)",
            "Entry", "Is_reviewed", "SeqLen", "SeqHash", "Entry_has_multiple_sequences", "Isoform detected",
            "Length", "E/Q ratio", "E/Q valid",
        ] + AA_ORDER
        cols_all = [c for c in cols_all if c in df_all.columns and c != "Chromosome/scaffold name"]

        preview_and_download(
            df=df_all,
            preview_cols=cols_all,
            title="Table — Protein-level AA counts + E/Q (one row per accession)",
            base_filename=f"chromosome_{chrom_ok}_EQ_full",
            preview_default=10,
        )

        invalid_mask = ~df_all["E/Q valid"]
        n_invalid = int(invalid_mask.sum())
        n_total = int(len(df_all))
        n_valid = n_total - n_invalid

        st.subheader("✅/❌ E/Q usability")
        st.write(f"Total proteins: **{n_total}**")
        st.write(f"Proteins with **non-usable** E/Q (Q = 0): **{n_invalid}**")
        st.write(f"Proteins with **valid** E/Q: **{n_valid}**")

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
            r, p = pearsonr(corr_df["Length"], df_valid.loc[corr_df.index, "E/Q ratio"])
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

        st.success("✅ Step 3 (E/Q) completed.")


# =========================================================
# FINAL Excel export
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
            excel_bytes = to_excel_bytes(
                df_all=st.session_state["df_all"],
                df_valid=st.session_state["df_valid"],
                avg_comp_all=st.session_state["avg_comp_all"],
                eq_summary_df=st.session_state["eq_summary_df"],
                walk_df=st.session_state.get("walk_df"),
                eq_window_df=st.session_state.get("eq_window_df"),
                fig_hist=st.session_state.get("fig_hist"),
                fig_corr=st.session_state.get("fig_corr"),
                fig_walk=st.session_state.get("fig_walk"),
                fig_eq_window=st.session_state.get("fig_eq_window"),
            )

        st.download_button(
            label="⬇️ Download Excel",
            data=excel_bytes,
            file_name=f"chromosome_{st.session_state['chrom_ok']}_final.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )
