import io
from typing import Tuple

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import streamlit as st
from pybiomart import Dataset
from bioservices import UniProt
from Bio.SeqUtils.ProtParam import ProteinAnalysis


st.set_page_config(page_title="Chromosome Protein Analyzer", layout="wide")

AA_COLUMNS = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
]

# Stable Ensembl archive. Ensembl 115 was released in September 2025.
# pybiomart expects an HTTP-style BioMart host URL here.
ENSEMBL_HOST = "http://sep2025.archive.ensembl.org"


def init_session_state() -> None:
    defaults = {
        "analysis_done": False,
        "show_heatmap": False,
        "canonical_df": None,
        "df_ensembl": None,
        "raw_uniprot": "",
        "last_params": {},
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def clean_sequence(seq: str) -> str:
    if pd.isna(seq):
        return ""
    return str(seq).rstrip("*").strip()


def validate_chromosome(chromosome: str) -> str:
    """Validate and normalize a human chromosome name.

    Accepted values are 1-22, X, Y and MT.
    Inputs such as chr1, chrX and chrMT are normalized.
    """
    if chromosome is None:
        return ""

    chrom = str(chromosome).strip().upper()

    if chrom.startswith("CHR"):
        chrom = chrom[3:]

    valid_chromosomes = {str(i) for i in range(1, 23)}
    valid_chromosomes.update({"X", "Y", "MT"})

    if chrom in valid_chromosomes:
        return chrom

    return ""


@st.cache_data(show_spinner=False)
def fetch_ensembl_data(chromosome_number: str) -> pd.DataFrame:
    dataset = Dataset(
        name="hsapiens_gene_ensembl",
        host=ENSEMBL_HOST,
    )

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
        filters={"chromosome_name": [chromosome_number]},
    )

    if df.empty:
        return df

    df = df.copy()

    if "UniProtKB/Swiss-Prot ID" in df.columns:
        df["UniProtKB/Swiss-Prot ID"] = df["UniProtKB/Swiss-Prot ID"].replace("", np.nan)

    if "Peptide" in df.columns:
        df["Peptide"] = df["Peptide"].apply(clean_sequence)

    df = df.dropna(subset=["UniProtKB/Swiss-Prot ID", "Peptide"])
    df = df[df["Peptide"] != ""].copy()
    df["length"] = df["Peptide"].str.len()

    df.rename(
        columns={
            "UniProtKB/Swiss-Prot ID": "id",
            "Peptide": "sequence",
        },
        inplace=True,
    )
    return df


@st.cache_data(show_spinner=False)
def fetch_uniprot_data(chromosome_number: str) -> Tuple[pd.DataFrame, str]:
    service = UniProt()

    query = (
        f'proteomecomponent:"chromosome {chromosome_number}" '
        f'AND organism_id:"9606" '
        f'AND proteome:up000005640 '
        f'AND reviewed:true'
    )

    result = service.search(query, frmt="tsv", columns="accession,id,length,sequence")

    if not result or not result.strip():
        return pd.DataFrame(columns=["id", "entry name", "length", "sequence"]), result

    df = pd.read_table(io.StringIO(result))

    rename_map = {}
    if "Entry" in df.columns:
        rename_map["Entry"] = "id"
    if "Entry Name" in df.columns:
        rename_map["Entry Name"] = "entry name"
    if "Length" in df.columns:
        rename_map["Length"] = "length"
    if "Sequence" in df.columns:
        rename_map["Sequence"] = "sequence"

    df = df.rename(columns=rename_map)

    for col in ["id", "entry name", "length", "sequence"]:
        if col not in df.columns:
            df[col] = np.nan

    df = df[["id", "entry name", "length", "sequence"]].copy()
    df["sequence"] = df["sequence"].apply(clean_sequence)
    return df, result


def merge_ensembl_uniprot(df_ensembl: pd.DataFrame, df_uniprot: pd.DataFrame) -> pd.DataFrame:
    if df_ensembl.empty or df_uniprot.empty:
        return pd.DataFrame()

    merged = pd.merge(
        df_ensembl,
        df_uniprot,
        on=["id", "length", "sequence"],
        how="inner",
    ).copy()

    return merged


def compute_amino_acid_counts(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=AA_COLUMNS)

    rows = []
    for seq in df["sequence"]:
        counts = ProteinAnalysis(seq).count_amino_acids()
        rows.append(counts)

    aa_df = pd.DataFrame(rows)
    for col in AA_COLUMNS:
        if col not in aa_df.columns:
            aa_df[col] = 0

    return aa_df[AA_COLUMNS].copy()


def build_canonical_table(merged_df: pd.DataFrame, aa_df: pd.DataFrame) -> pd.DataFrame:
    if merged_df.empty:
        return pd.DataFrame()

    final_df = pd.concat(
        [merged_df.reset_index(drop=True), aa_df.reset_index(drop=True)],
        axis=1,
    )

    desired_columns = [
        "Gene stable ID",
        "Gene name",
        "Gene start (bp)",
        "Gene end (bp)",
        "id",
        "length",
        "entry name",
        *AA_COLUMNS,
        "Chromosome/scaffold name",
        "sequence",
    ]

    available = [c for c in desired_columns if c in final_df.columns]
    final_df = final_df[available].copy()

    if "Gene start (bp)" in final_df.columns:
        final_df = final_df.sort_values(by="Gene start (bp)").reset_index(drop=True)

    return final_df


def _safe_join(values: pd.Series) -> str:
    vals = [str(v) for v in values.tolist() if pd.notna(v) and str(v).strip() != ""]
    return "; ".join(vals)


def build_walking_table(
    canonical_df: pd.DataFrame,
    aa_selected: str,
    window_size: int,
    threshold: int,
) -> pd.DataFrame:
    if canonical_df.empty or aa_selected not in canonical_df.columns or len(canonical_df) < window_size:
        return pd.DataFrame()

    rows = []
    values = canonical_df[aa_selected].reset_index(drop=True)

    for start_idx in range(0, len(canonical_df) - window_size + 1):
        end_idx = start_idx + window_size - 1
        block = canonical_df.iloc[start_idx:end_idx + 1]
        block_values = values.iloc[start_idx:end_idx + 1]
        walking_count = int((block_values > threshold).sum())

        rows.append(
            {
                "AA": aa_selected,
                "window_size": int(window_size),
                "threshold": int(threshold),
                "frame_number": start_idx + 1,
                "start_block_element": start_idx + 1,
                "final_block_element": end_idx + 1,
                "block_start_bp": int(block["Gene start (bp)"].iloc[0]),
                "block_end_bp": int(block["Gene end (bp)"].iloc[-1]),
                "walking_count": walking_count,
                "first_gene_name": block["Gene name"].iloc[0] if "Gene name" in block.columns else "",
                "last_gene_name": block["Gene name"].iloc[-1] if "Gene name" in block.columns else "",
                "gene_names_in_window": _safe_join(block["Gene name"]) if "Gene name" in block.columns else "",
                "gene_ids_in_window": _safe_join(block["Gene stable ID"]) if "Gene stable ID" in block.columns else "",
            }
        )

    return pd.DataFrame(rows)


def build_eq_ratio_table(canonical_df: pd.DataFrame, window_size: int, threshold: int) -> pd.DataFrame:
    walking_e = build_walking_table(canonical_df, "E", window_size, threshold)
    walking_q = build_walking_table(canonical_df, "Q", window_size, threshold)

    if walking_e.empty or walking_q.empty:
        return pd.DataFrame()

    eq_df = walking_e[
        [
            "frame_number",
            "start_block_element",
            "final_block_element",
            "block_start_bp",
            "block_end_bp",
            "first_gene_name",
            "last_gene_name",
            "gene_names_in_window",
            "gene_ids_in_window",
        ]
    ].copy()
    eq_df["E_count"] = walking_e["walking_count"].values
    eq_df["Q_count"] = walking_q["walking_count"].values
    q_counts = eq_df["Q_count"].replace(0, np.nan)
    eq_df["E_Q_ratio"] = (eq_df["E_count"] / q_counts).replace([np.inf, -np.inf], np.nan)
    return eq_df


def build_walking_heatmap(canonical_df: pd.DataFrame, window_size: int, threshold: int):
    if canonical_df.empty or len(canonical_df) < window_size:
        return None

    heatmap_data = []
    valid_aa = []

    for aa in AA_COLUMNS:
        walking_df = build_walking_table(canonical_df, aa, window_size, threshold)
        if walking_df.empty:
            continue
        heatmap_data.append(walking_df["walking_count"].values / window_size)
        valid_aa.append(aa)

    if not heatmap_data:
        return None

    heatmap_array = np.array(heatmap_data)
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(heatmap_array, aspect="auto", vmin=0, vmax=1)
    ax.set_yticks(range(len(valid_aa)))
    ax.set_yticklabels(valid_aa)
    ax.set_xlabel("Frame number")
    ax.set_ylabel("Amino acid")
    ax.set_title("Walking heatmap (frequency of proteins with AA > T)")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Frequency (0-1)")
    return fig


def dataframe_to_csv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8")


def dataframe_to_xlsx_bytes(df: pd.DataFrame, sheet_name: str = "Sheet1") -> bytes:
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name[:31])
    buffer.seek(0)
    return buffer.getvalue()


def build_figure_b(walking_df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(
        walking_df["frame_number"],
        walking_df["walking_count"],
        marker="o",
        linewidth=1.8,
    )
    ax.set_xlabel("Frame number")
    ax.set_ylabel("Number of proteins with X > T")
    ax.set_title("Walking analysis")
    ax.set_xticks(walking_df["frame_number"])
    ax.grid(True, axis="y")
    return fig


def build_figure_a(canonical_df: pd.DataFrame, aa_selected: str, threshold: int, window_size: int):
    if canonical_df.empty:
        return None

    n_show = len(canonical_df)
    values = canonical_df[aa_selected].reset_index(drop=True)
    enriched = values > threshold

    n_frames = max(0, n_show - window_size + 1)
    if n_frames == 0:
        fig, ax = plt.subplots(figsize=(10, 2.5))
        ax.text(0.5, 0.5, "Not enough proteins for schematic walking panel.", ha="center", va="center")
        ax.axis("off")
        return fig

    frame_counts = []
    for start_idx in range(n_frames):
        end_idx = start_idx + window_size - 1
        frame_counts.append(int((values.iloc[start_idx:end_idx + 1] > threshold).sum()))
    max_frame = int(np.argmax(frame_counts))

    fig_w = min(28, max(14, 0.22 * n_show + 6))
    fig_h = min(28, max(5, 2.5 + 0.18 * n_frames))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    y_top = n_frames + 2
    row_font = max(5, min(9, 10 - n_frames * 0.06))
    label_font = max(5, min(9, 10 - n_show * 0.02))

    for i in range(n_show):
        face = "red" if enriched.iloc[i] else "white"
        ax.add_patch(
            Rectangle(
                (i + 1 - 0.45, y_top),
                0.9,
                0.42,
                facecolor=face,
                edgecolor="black",
                linewidth=0.8,
            )
        )

    ax.text((n_show + 1) / 2, y_top + 0.68, "ORDERED LIST OF PROTEINS CODIFIED BY A CHROMOSOME",
            ha="center", va="bottom", fontsize=label_font)
    ax.text(-1.6, y_top + 0.2, "FRAME\nNUMBER", ha="center", va="center", fontsize=label_font)

    for frame_idx in range(n_frames):
        row_y = n_frames - frame_idx + 1
        start = frame_idx + 1
        end = frame_idx + window_size
        color = "red" if frame_idx == max_frame else "black"
        count = frame_counts[frame_idx]
        ax.text(-0.2, row_y, str(frame_idx + 1), ha="right", va="center", fontsize=row_font, color=color)
        ax.plot([start, start], [y_top, row_y], color="gray", linewidth=1)
        ax.plot([end, end], [y_top, row_y], color="gray", linewidth=1)
        ax.plot([start, end], [row_y, row_y], color="gray", linewidth=1)
        ax.text((start + end) / 2, row_y + 0.1, f"{count} proteins X > T",
                ha="center", va="bottom", fontsize=row_font, color=color)

    ax.set_xlim(-2, n_show + 2)
    ax.set_ylim(0.5, y_top + 1.2)
    ax.axis("off")
    ax.set_title("Figure A style schematic walking panel")
    return fig


def build_aa_profile_plot(canonical_df: pd.DataFrame, aa_selected: str, threshold: int):
    profile_df = canonical_df.reset_index(drop=True).copy()
    profile_df["protein_order"] = np.arange(1, len(profile_df) + 1)
    enriched = profile_df[aa_selected] > threshold

    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.plot(profile_df["protein_order"], profile_df[aa_selected], marker="o", linewidth=1.4)
    ax.axhline(y=threshold, linestyle="--", linewidth=1.2)

    if enriched.any():
        ax.scatter(profile_df.loc[enriched, "protein_order"], profile_df.loc[enriched, aa_selected], s=36, zorder=3)

    ax.set_xlabel("Proteins ordered along the chromosome")
    ax.set_ylabel(f"Absolute content of {aa_selected}")
    ax.set_title(f"{aa_selected} content profile with threshold T={threshold}")
    ax.grid(True, axis="y")
    return fig


def build_eq_comparison_plot(eq_df: pd.DataFrame):
    if eq_df.empty:
        return None

    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(eq_df["frame_number"], eq_df["E_count"], marker="o", linewidth=1.8, label="E")
    ax1.plot(eq_df["frame_number"], eq_df["Q_count"], marker="s", linewidth=1.8, label="Q")
    ax1.set_xlabel("Frame number")
    ax1.set_ylabel("Number of proteins with X > T")
    ax1.grid(True, axis="y")

    ax2 = ax1.twinx()
    ax2.plot(eq_df["frame_number"], eq_df["E_Q_ratio"], marker="^", linewidth=1.5, linestyle="--", label="E/Q")
    ax2.set_ylabel("E/Q ratio")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")
    ax1.set_title("Walking comparison: E, Q and E/Q ratio")
    return fig


def run_analysis(chromosome_number: str):
    df_ensembl = fetch_ensembl_data(chromosome_number)
    df_uniprot, raw_uniprot = fetch_uniprot_data(chromosome_number)
    merged_df = merge_ensembl_uniprot(df_ensembl, df_uniprot)
    aa_df = compute_amino_acid_counts(merged_df)
    canonical_df = build_canonical_table(merged_df, aa_df)
    return df_ensembl, canonical_df, raw_uniprot


def main():
    init_session_state()

    st.title("Chromosome Protein Analyzer")
    st.write(
        "Canonical table, downloadable WALKING tables, Figure B-style plot, "
        "Figure A-style schematic panel with all iterations, E/Q/EQ comparison, "
        "AA content profile with T line and normalized walking heatmap."
    )

    with st.sidebar:
        st.header("Parameters")
        chromosome_number = st.text_input("Chromosome number", value="1")
        aa_selected = st.selectbox("Amino acid", AA_COLUMNS, index=3)
        window_size = st.number_input("Window size (LB)", min_value=2, max_value=500, value=25, step=1)
        threshold = st.number_input("Threshold (T)", min_value=0, max_value=10000, value=49, step=1)
        run_button = st.button("Run analysis", use_container_width=True)

    validated_chromosome = validate_chromosome(chromosome_number)

    params = {
        "chromosome_number": validated_chromosome,
        "aa_selected": aa_selected,
        "window_size": int(window_size),
        "threshold": int(threshold),
    }

    if run_button:
        if not params["chromosome_number"]:
            st.error("Invalid chromosome. Please enter one of: 1-22, X, Y or MT. Inputs such as chr1 or chrX are accepted.")
            st.session_state.analysis_done = False
            return

        with st.spinner("Retrieving data from Ensembl and UniProt..."):
            try:
                df_ensembl, canonical_df, raw_uniprot = run_analysis(params["chromosome_number"])
            except Exception as e:
                st.error(f"Error during analysis: {e}")
                st.session_state.analysis_done = False
                return

        st.session_state.df_ensembl = df_ensembl
        st.session_state.canonical_df = canonical_df
        st.session_state.raw_uniprot = raw_uniprot
        st.session_state.analysis_done = True
        st.session_state.show_heatmap = False
        st.session_state.last_params = params.copy()

    if not st.session_state.analysis_done:
        st.info("Set the parameters in the sidebar and run the analysis.")
        return

    canonical_df = st.session_state.canonical_df
    df_ensembl = st.session_state.df_ensembl
    raw_uniprot = st.session_state.raw_uniprot
    stored_params = st.session_state.last_params or params

    if canonical_df is None or canonical_df.empty:
        st.warning("No data available after merging Ensembl and UniProt.")
        return

    st.subheader(f"Results for chromosome {stored_params['chromosome_number']}")

    walking_df = build_walking_table(
        canonical_df=canonical_df,
        aa_selected=stored_params["aa_selected"],
        window_size=stored_params["window_size"],
        threshold=stored_params["threshold"],
    )

    eq_df = build_eq_ratio_table(
        canonical_df=canonical_df,
        window_size=stored_params["window_size"],
        threshold=stored_params["threshold"],
    )

    profile_preview = canonical_df[
        ["Gene name", "Gene stable ID", "Gene start (bp)", "Gene end (bp)", stored_params["aa_selected"]]
    ].copy()
    profile_preview["above_threshold"] = profile_preview[stored_params["aa_selected"]] > stored_params["threshold"]
    profile_preview.insert(0, "protein_order", np.arange(1, len(profile_preview) + 1))

    source_export = df_ensembl.copy() if df_ensembl is not None else pd.DataFrame()

    tab1, tab2, tab3, tab4 = st.tabs([
        "Canonical table",
        "Walking and plots",
        "AA profile",
        "Source data",
    ])

    with tab1:
        st.dataframe(canonical_df, use_container_width=True)
        st.write(f"Number of rows: {len(canonical_df)}")
        c1, c2 = st.columns(2)
        with c1:
            st.download_button("Download canonical table CSV",
                data=dataframe_to_csv_bytes(canonical_df),
                file_name=f"CHR{stored_params['chromosome_number'].zfill(2)}_canonical.csv",
                mime="text/csv")
        with c2:
            st.download_button("Download canonical table Excel",
                data=dataframe_to_xlsx_bytes(canonical_df, sheet_name="canonical"),
                file_name=f"CHR{stored_params['chromosome_number'].zfill(2)}_canonical.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

    with tab2:
        if walking_df.empty:
            st.warning("Insufficient data for walking with the selected parameters.")
        else:
            st.write(
                f"Walking on AA = {stored_params['aa_selected']}, "
                f"LB = {stored_params['window_size']}, T = {stored_params['threshold']}"
            )
            fig_a = build_figure_a(canonical_df, stored_params["aa_selected"], stored_params["threshold"], stored_params["window_size"])
            if fig_a is not None:
                st.pyplot(fig_a)
            st.pyplot(build_figure_b(walking_df))
            fig_eq = build_eq_comparison_plot(eq_df)
            if fig_eq is not None:
                st.pyplot(fig_eq)

            c1, c2 = st.columns(2)
            with c1:
                st.download_button("Download WALKING table CSV",
                    data=dataframe_to_csv_bytes(walking_df),
                    file_name=f"CHR{stored_params['chromosome_number'].zfill(2)}_canonical_WALKING.csv",
                    mime="text/csv")
            with c2:
                st.download_button("Download WALKING table Excel",
                    data=dataframe_to_xlsx_bytes(walking_df, sheet_name="walking"),
                    file_name=f"CHR{stored_params['chromosome_number'].zfill(2)}_canonical_WALKING.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

            if not eq_df.empty:
                c3, c4 = st.columns(2)
                with c3:
                    st.download_button("Download E_Q table CSV",
                        data=dataframe_to_csv_bytes(eq_df),
                        file_name=f"CHR{stored_params['chromosome_number'].zfill(2)}_EQ_ratio.csv",
                        mime="text/csv")
                with c4:
                    st.download_button("Download E_Q table Excel",
                        data=dataframe_to_xlsx_bytes(eq_df, sheet_name="EQ_ratio"),
                        file_name=f"CHR{stored_params['chromosome_number'].zfill(2)}_EQ_ratio.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

            with st.expander("WALKING table preview"):
                st.dataframe(walking_df, use_container_width=True)
            if not eq_df.empty:
                with st.expander("E/Q table preview"):
                    st.dataframe(eq_df, use_container_width=True)

            st.markdown("---")
            if st.button("Generate normalized walking heatmap", key="generate_heatmap"):
                st.session_state.show_heatmap = True
            if st.session_state.show_heatmap:
                fig_heatmap = build_walking_heatmap(canonical_df, stored_params["window_size"], stored_params["threshold"])
                if fig_heatmap is None:
                    st.warning("Unable to generate the heatmap with the selected parameters.")
                else:
                    st.pyplot(fig_heatmap)

    with tab3:
        st.pyplot(build_aa_profile_plot(canonical_df, stored_params["aa_selected"], stored_params["threshold"]))
        c1, c2 = st.columns(2)
        with c1:
            st.download_button("Download AA profile CSV",
                data=dataframe_to_csv_bytes(profile_preview),
                file_name=f"CHR{stored_params['chromosome_number'].zfill(2)}_AA_profile.csv",
                mime="text/csv")
        with c2:
            st.download_button("Download AA profile Excel",
                data=dataframe_to_xlsx_bytes(profile_preview, sheet_name="AA_profile"),
                file_name=f"CHR{stored_params['chromosome_number'].zfill(2)}_AA_profile.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        with st.expander("Profile data preview"):
            st.dataframe(profile_preview, use_container_width=True)

    with tab4:
        st.write("### Ensembl data")
        if df_ensembl is None or df_ensembl.empty:
            st.warning("No Ensembl data found.")
        else:
            st.dataframe(df_ensembl, use_container_width=True)

        st.write("### UniProt output")
        if raw_uniprot and raw_uniprot.strip():
            st.text(raw_uniprot)
        else:
            st.warning("No UniProt output available.")

        c1, c2 = st.columns(2)
        with c1:
            st.download_button("Download source data CSV",
                data=dataframe_to_csv_bytes(source_export),
                file_name=f"CHR{stored_params['chromosome_number'].zfill(2)}_ensembl_source.csv",
                mime="text/csv")
        with c2:
            st.download_button("Download source data Excel",
                data=dataframe_to_xlsx_bytes(source_export, sheet_name="ensembl_source"),
                file_name=f"CHR{stored_params['chromosome_number'].zfill(2)}_ensembl_source.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")


if __name__ == "__main__":
    main()
