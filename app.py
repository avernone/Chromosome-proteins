import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import io
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.stats import pearsonr

# Streamlit setup
st.set_page_config(page_title="Chromosome E/Q Ratio Explorer", page_icon="🧬", layout="wide")
st.title("🧬 Chromosome Amino Acid Composition and E/Q Ratio Explorer (Enhanced Version)")

st.write("""
Enter a **human chromosome number (1–22, X, or Y)** to analyze all reviewed UniProt entries for that chromosome.

The app:
1. Retrieves all **reviewed (Swiss-Prot)** human proteins for the selected chromosome  
2. Computes amino acid composition for each sequence  
3. Calculates the **E/Q ratio (glutamate/glutamine)**  
4. Displays summary tables including genomic coordinates and amino acid counts  
5. Evaluates **correlation between protein length and E/Q ratio**  
6. Exports results, statistics and plots in Excel format
""")

# --- Input ---
chromosome = st.text_input("Chromosome number (1–22, X, Y):", "").strip()

# --- Cache API calls for speed ---
@st.cache_data(show_spinner=False)
def fetch_uniprot_data(chrom):
    """Fetch UniProt reviewed human proteins by chromosome."""
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    query = f'(proteomecomponent:"chromosome {chrom}") AND (organism_id:9606) AND (reviewed:true)'
    params = {"query": query, "format": "tsv", "fields": "accession,id,length,sequence,ft_gene,ft_gene_start,ft_gene_end"}
    r = requests.get(base_url, params=params, timeout=60)
    if r.status_code == 200 and len(r.text.splitlines()) > 1:
        return pd.read_csv(io.StringIO(r.text), sep="\t")
    else:
        return pd.DataFrame()

# --- Compute amino acid composition & E/Q ratio ---
def compute_composition(df):
    compositions, eq_ratios = [], []
    for seq in df["Sequence"]:
        try:
            analysis = ProteinAnalysis(seq)
            counts = analysis.count_amino_acids()
            e = counts.get("E", 0)
            q = counts.get("Q", 0)
            eq = e / q if q != 0 else np.nan
        except Exception:
            counts = {aa: np.nan for aa in list("ACDEFGHIKLMNPQRSTVWY")}
            eq = np.nan
        compositions.append(counts)
        eq_ratios.append(eq)
    comp_df = pd.DataFrame(compositions)
    df["E/Q ratio"] = eq_ratios
    return pd.concat([df, comp_df], axis=1)

# --- Run analysis ---
if chromosome:
    with st.spinner(f"Retrieving reviewed UniProt proteins for chromosome {chromosome}..."):
        df = fetch_uniprot_data(chromosome)

        if df.empty:
            st.error("❌ No reviewed UniProt entries found for this chromosome.")
        else:
            total_proteins = len(df)
            st.success(f"✅ Retrieved {total_proteins} protein entries.")

            # Compute amino acid composition and E/Q ratio
            df_full = compute_composition(df)
            df_full.dropna(subset=["E/Q ratio"], inplace=True)
            analyzed_proteins = len(df_full)

            st.caption(f"⚙️ Total retrieved proteins: {total_proteins} — proteins with valid E/Q ratio: {analyzed_proteins}")

            # Sort by E/Q ratio
            df_full.sort_values(by="E/Q ratio", inplace=True)

            # --- Summary table ---
            st.subheader(f"📊 Proteins and E/Q ratios (Chromosome {chromosome})")
            st.dataframe(df_full[
                [
                    "Gene Start", "Gene End", "Entry", "Entry Name", "Length", "E/Q ratio"
                ] + list("ACDEFGHIKLMNPQRSTVWY")
            ])

            # --- Average amino acid composition ---
            st.subheader(f"📈 Average amino acid composition")
            avg_comp = df_full[list("ACDEFGHIKLMNPQRSTVWY")].mean()
            st.bar_chart(avg_comp)

            # --- Summary statistics of E/Q ratio ---
            st.subheader("📋 E/Q ratio summary statistics")

            eq_summary = {
                "Mean": df_full["E/Q ratio"].mean(),
                "Median": df_full["E/Q ratio"].median(),
                "Standard deviation": df_full["E/Q ratio"].std(),
                "Minimum": df_full["E/Q ratio"].min(),
                "Maximum": df_full["E/Q ratio"].max(),
                "Count": analyzed_proteins
            }
            eq_summary_df = pd.DataFrame(eq_summary, index=["Value"]).T
            st.table(eq_summary_df)

            # --- Correlation analysis: E/Q ratio vs Length ---
            st.subheader("🔗 Correlation between protein length and E/Q ratio")
            corr_df = df_full.dropna(subset=["Length", "E/Q ratio"])
            if len(corr_df) > 2:
                r, p = pearsonr(corr_df["Length"], corr_df["E/Q ratio"])
                st.write(f"**Pearson correlation (r)** = `{r:.3f}` (p = {p:.3e})")

                # Scatter plot with regression line
                fig_corr, ax_corr = plt.subplots(figsize=(8, 5))
                sns.regplot(
                    data=corr_df,
                    x="Length",
                    y="E/Q ratio",
                    scatter_kws={"alpha": 0.6, "color": "skyblue"},
                    line_kws={"color": "red"}
                )
                ax_corr.set_title(f"Correlation between protein length and E/Q ratio (r = {r:.2f})")
                st.pyplot(fig_corr)
            else:
                st.warning("Not enough data points to compute correlation.")

            # --- Plot distribution of E/Q ratios ---
            fig, ax = plt.subplots(figsize=(10, 5))
            sns.histplot(df_full["E/Q ratio"], bins=30, kde=True, color="skyblue", ax=ax)
            ax.set_title(f"E/Q ratio distribution (Chromosome {chromosome})")
            ax.set_xlabel("E/Q ratio")
            ax.set_ylabel("Protein count")
            st.pyplot(fig)

            # --- Export Excel ---
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
                df_full.to_excel(writer, index=False, sheet_name="ChromosomeData")
                avg_comp.to_excel(writer, sheet_name="AverageComposition")
                eq_summary_df.to_excel(writer, sheet_name="EQ_Summary")

                if len(corr_df) > 2:
                    corr_summary = pd.DataFrame({
                        "Statistic": ["Pearson r", "p-value"],
                        "Value": [r, p]
                    })
                    corr_summary.to_excel(writer, sheet_name="Correlation", index=False)

                # Insert plots
                fig_buf = io.BytesIO()
                fig.savefig(fig_buf, format="png", bbox_inches="tight")
                fig_buf.seek(0)
                if len(corr_df) > 2:
                    corr_buf = io.BytesIO()
                    fig_corr.savefig(corr_buf, format="png", bbox_inches="tight")
                    corr_buf.seek(0)

                worksheet = writer.sheets["ChromosomeData"]
                worksheet.insert_image("I2", "eq_distribution.png", {"image_data": fig_buf})
                if len(corr_df) > 2:
                    writer.sheets["Correlation"].insert_image("D5", "correlation.png", {"image_data": corr_buf})

            st.download_button(
                label="⬇️ Download Excel (Results + Statistics + Correlation)",
                data=output.getvalue(),
                file_name=f"chromosome_{chromosome}_EQ_summary.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            st.success("✅ Analysis completed successfully.")

