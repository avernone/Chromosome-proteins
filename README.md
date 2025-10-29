# GPLAB - Gene and Protein Virtual LAb

In this Lab, the effect of amino acid availability as regulating factor of every protein synthesis is investigated

Contributors:

**Gianpiero Pescarmona**

Formerly Department of Oncology, University of Torino

**Francesca Silvagno**

Department of Oncology, University of Torino

**Annamaria Vernone**

Department of Neurosciences "Rita Levi Montalcini", University of Torino

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
