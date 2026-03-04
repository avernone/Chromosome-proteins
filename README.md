# GPLAB - Gene and Protein Virtual LAb

In this Lab, the effect of amino acid availability as regulating factor of every protein synthesis is investigated

Contributors:

**Gianpiero Pescarmona**

Formerly Department of Oncology, University of Torino

**Francesca Silvagno**

Department of Oncology, University of Torino

**Annamaria Vernone**

Department of Neurosciences "Rita Levi Montalcini", University of Torino

Overview

This application analyzes amino acid composition and E/Q ratios of human proteins encoded on a selected chromosome (1–22, X, Y).

The app retrieves data from Ensembl BioMart (pybiomart) and performs protein-level and chromosome-level compositional analyses.

Data Source

Data are retrieved from:

Ensembl BioMart (hsapiens_gene_ensembl)

For the selected chromosome, the app collects:
 
Ensembl Gene ID

Gene name

Gene start and end positions (bp)

UniProt Swiss-Prot accession

Peptide sequence

Only entries with:

Swiss-Prot mapping

Available peptide sequence are included in the analysis.

Analysis Workflow
1. Full Protein Table (No Filtering)

The first table displays all proteins found on the selected chromosome.

Displayed columns include:

Gene stable ID

Gene name

Gene start (bp)

Gene end (bp)

UniProt accession

Protein length

Amino acid absolute counts (A–Y)

E/Q ratio

E/Q valid (True/False)

⚠️ No filtering is applied at this stage.

2. E/Q Usability Check

The application reports:

Total proteins found

Number of proteins with non-usable E/Q (Q = 0)

Number of proteins with valid E/Q

Rows with invalid E/Q can be inspected separately.

3. Average Amino Acid Composition

Average amino acid composition is computed using all proteins found, not only those with valid E/Q.

This provides a chromosome-wide compositional overview.

4. E/Q Ratio Summary Statistics

E/Q statistics are computed only on proteins with valid E/Q.

Reported metrics:

Mean

Median

Standard deviation

Minimum

Maximum

Count of valid E/Q proteins

Total proteins

Invalid E/Q count

5. Correlation Analysis

Pearson correlation between:

Protein length

E/Q ratio

Computed only on valid E/Q proteins.

6. E/Q Distribution

Histogram + KDE curve of E/Q ratios (valid proteins only).

Walking (Surfing) Analysis
Concept

The Walking (Surfing) analysis scans the chromosome using a sliding window of fixed size (number of proteins).

Proteins are ordered by Gene start (bp).

For each frame (window), the software counts how many proteins have:

Selected amino acid (AA) count > threshold T

Glutamate (E) count > T

Glutamine (Q) count > T

⚠️ E/Q ratio is intentionally excluded from walking because it requires a different threshold logic.

Walking Parameters

User-defined parameters:

Amino acid (AA) to analyze

Frame size (number of proteins per window)

Threshold T (absolute AA count > T)

Walking Output Table

For each sliding window, the table includes:

Frame number

Frame start bp (start position of first protein in frame)

Frame end bp (start position of last protein in frame)

Surfing_<AA> (count of proteins with AA > T)

Surfing_E

Surfing_Q

Protein_1 ... Protein_N (accessions within the frame)

Number of rows:

N − window + 1

Where N = total proteins.

Chromosome-wide Walking Visualization

A line plot provides an immediate visual interpretation of compositional clustering.

User options:

X-axis:

Frame start bp (genomic coordinate view)

Frame number (ordinal view)

Optional smoothing:

Rolling mean over selected number of frames

The plot shows:

Surfing_<AA>

Surfing_E

Surfing_Q

Peaks indicate local clusters of proteins enriched in the selected amino acid.

Biological Interpretation

The walking analysis is:

Protein-driven (fixed number of proteins per frame)

Not fixed genomic length

Therefore:

Dense genomic regions → narrower bp span

Sparse regions → wider bp span

Peaks in Surfing curves suggest:

Local compositional clustering

Possible functional grouping

Regional structural bias

Potential genomic organization effects

Export Functionality

The Excel export includes:

All proteins table

Valid E/Q subset

Average AA composition (all proteins)

E/Q summary statistics

Walking table

E/Q histogram

Correlation plot

Walking chromosome plot

Requirements

Install dependencies:

streamlit
pandas
numpy
matplotlib
seaborn
biopython
scipy
pybiomart
xlsxwriter

Run locally:

streamlit run app.py
Key Design Choices

Ensembl used instead of UniProt chromosome filtering (more stable and genome-driven)

Full dataset shown before filtering

Clear separation between:

Compositional analysis (all proteins)

Ratio-based analysis (valid subset)

Walking excludes E/Q due to threshold inconsistency

Visualization included for immediate chromosome-scale interpretation

The application accesses UniProt (https://www.uniprot.org/) and Ensembl (https://www.ensembl.org/index.html) databases programmatically, using python interfaces in order to obtain data always updated. By this automated procedure, the amino acid content inserted in the canonical data table is always computed starting from the latest data available online from the two databases.


The Uniprot Table and the Ensembl table were merged based on Uniprot/Swissprot id accession number (AC), length, and sequence in order to obtain a Canonical table (https://www.uniprot.org/help/canonical_and_isoforms). 

The procedure and applications are described in:

**Chromosome Walking: A Novel Approach to Analyse Amino Acid Content of Human Proteins Ordered by Gene Position** Appl. Sci. 2021, 11(8), 3511, Published online 2021 Apr 14, https://doi.org/10.3390/app11083511 (https://www.mdpi.com/2076-3417/11/8/3511).

Enter a human chromosome number (1–22, X, or Y) to retrieve all protein-coding genes.
The app will:

Retrieve protein sequences from Ensembl and UniProt
Compute amino acid composition for each gene
Calculate the E/Q ratio (glutamate/glutamine)
Display summary tables and visualizations
Export results and plots to Excel
