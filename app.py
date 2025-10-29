import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from bioservices import UniProt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pybiomart import Server, Dataset
import io

# --- Streamlit setup ---
st.set_page_config(page_title="Chromosome Amino Acid & E/Q Analysis", page_icon="🧬", layout="wide")
st.title("🧬 Chromosome Amino Acid Composition and E/Q Ratio Explorer")

st.write("""
Enter a **human chromosome number (1–22, X, or Y)** to retrieve all protein-coding genes.  
The app will:
1. Retrieve protein sequences from Ensembl and UniProt  
2. Compute amino acid composition for each gene  
3. Calculate the **E/Q ratio** (glutamate/glutamine)  
4. Display summary tables and visualizations  
5. Export results and plots to Excel
""")

# --- Input ---
chromosome_number = st.text_input("Chromosome number:", "")

if chromosome_number:
    with st.spinner("Retrieving and analyzing protein data..."):
        try:
            # --- Step 1. Retrieve data from Ensembl ---
            server = Server(host='http://www.ensembl.org')
            dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')

            Ensembl_data = dataset.query(
                attributes=[
                    'ensembl_gene_id', 'external_gene_name', 'start_position', 'end_position',
                    'chromosome_name', 'uniprotswissprot', 'peptide'
                ],
                filters={'chromosome_name': [chromosome_number]}
            )

            # Clean Ensembl data
            Ensembl_data['UniProtKB/Swiss-Prot ID'].replace('', np.nan, inplace=True)
            Ensembl_data = Ensembl_data.dropna()
            Ensembl_data['Peptide'] = Ensembl_data['Peptide'].str.rstrip('*')
            Ensembl_data['length'] = Ensembl_data['Peptide'].apply(len)
            Ensembl_data.rename(columns={
                'UniProtKB/Swiss-Prot ID': 'id',
                'Peptide': 'sequence'
            }, inplace=True)

            # --- Step 2. Retrieve UniProt data for the same chromosome ---
            service = UniProt()
            query = f'proteomecomponent:"chromosome {chromosome_number}" AND organism_id:9606 AND proteome:up000005640 AND reviewed:true'
            cols = "accession,id,length,sequence"
            result_uniprot = service.search(query, frmt="tsv", columns=cols)

            df_uniprot = pd.read_table(io.StringIO(result_uniprot))
            df_uniprot.columns = ['id', 'entry name', 'length', 'sequence']

            # --- Step 3. Merge Ensembl + UniProt ---
            merged = pd.merge(Ensembl_data, df_uniprot, on=['id', 'length', 'sequence'], how='inner')

            # --- Step 4. Compute amino acid counts + E/Q ratios ---
            aa_counts = []
            eq_ratios = []

            for seq in merged['sequence']:
                try:
                    analysis = ProteinAnalysis(seq)
                    counts = analysis.count_amino_acids()
                    aa_counts.append(counts)
                    e = counts.get("E", 0)
                    q = counts.get("Q", 0)
                    eq_ratios.append(e / q if q != 0 else np.nan)
                except Exception:
                    aa_counts.append({aa: np.nan for aa in list("ACDEFGHIKLMNPQRSTVWY")})
                    eq_ratios.append(np.nan)

            df_counts = pd.DataFrame(aa_counts)
            merged["E/Q ratio"] = eq_ratios

            # --- Step 5. Combine and sort ---
            df_full = pd.concat([merged, df_counts], axis=1)
            df_full = df_full.sort_values(by=["Gene start (bp)"])

            st.subheader("Per-protein amino acid counts and E/Q ratios")
            st.dataframe(df_full)

            # --- Step 6. Average amino acid composition ---
            st.subheader(f"Average amino acid composition (Chromosome {chromosome_number})")
            avg_composition = df_counts.mean().sort_index()
            st.bar_chart(avg_composition)

            # --- Step 7. Sort by E/Q ratio and add progressive numbering ---
            eq_sorted = df_full.sort_values(by="E/Q ratio", ascending=True).dropna(subset=["E/Q ratio"])
            eq_sorted.insert(0, "No.", range(1, len(eq_sorted) + 1))

            st.subheader(f"Proteins sorted by E/Q ratio (Chromosome {chromosome_number})")
            st.dataframe(eq_sorted[["No.", "Gene name", "entry name", "E/Q ratio"]])

            # --- Step 8. Highlight key proteins (lowest, median, highest) ---
            min_row = eq_sorted.iloc[0]
            median_row = eq_sorted.iloc[len(eq_sorted)//2]
            max_row = eq_sorted.iloc[-1]

            # --- Step 9. Bar plot E/Q ratio ---
            fig, ax = plt.subplots(figsize=(12, 5))
            colors = sns.color_palette("coolwarm", len(eq_sorted))
            ax.bar(eq_sorted["No."], eq_sorted["E/Q ratio"], color=colors)

            for idx, row in [
                (0, min_row),
                (len(eq_sorted)//2, median_row),
                (len(eq_sorted)-1, max_row)
            ]:
                ax.text(
                    idx + 1,
                    row["E/Q ratio"] + 0.02,
                    f"{row['entry name']} ({row['E/Q ratio']:.2f})",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    rotation=45,
                    weight="bold",
                    color="black"
                )

            plt.ylabel("E/Q ratio (E / Q)")
            plt.xlabel("Protein index (No.)")
            plt.title(f"E/Q ratio distribution across Chromosome {chromosome_number}")
            plt.tight_layout()
            st.pyplot(fig)

            # --- Step 10. Summary info box ---
            st.info(f"""
**Summary for Chromosome {chromosome_number}:**
- Lowest E/Q: `{min_row['entry name']}` ({min_row['E/Q ratio']:.2f})
- Median E/Q: `{median_row['entry name']}` ({median_row['E/Q ratio']:.2f})
- Highest E/Q: `{max_row['entry name']}` ({max_row['E/Q ratio']:.2f})
""")

            # --- Step 11. Export Excel with plots ---
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
                df_full.to_excel(writer, index=False, sheet_name="ChromosomeData")
                eq_sorted.to_excel(writer, index=False, sheet_name="EQ_Sorted")

                # Save plot to Excel
                figfile = io.BytesIO()
                fig.savefig(figfile, format="png", bbox_inches="tight")
                figfile.seek(0)
                worksheet = writer.sheets["EQ_Sorted"]
                worksheet.insert_image("E2", "eq_plot.png", {"image_data": figfile})

            excel_data = output.getvalue()
            st.download_button(
                label="⬇️ Download Excel (with E/Q plot)",
                data=excel_data,
                file_name=f"chromosome_{chromosome_number}_AA_EQ_analysis.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            st.success("✅ Analysis completed successfully.")

        except Exception as e:
            st.error(f"❌ Error retrieving or processing data: {e}")

