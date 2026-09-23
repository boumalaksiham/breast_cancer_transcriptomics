# Breast Cancer Transcriptomic Signature Analysis — GSE42568

## Project question
Which gene-expression signals differ between breast-cancer tissue and normal breast tissue, and what biological processes are represented among the strongest signals?

## Public dataset
**Primary source:** NCBI Gene Expression Omnibus (GEO), accession **GSE42568**  
**Study title:** *Breast Cancer Gene Expression Analysis*  
**Organism:** *Homo sapiens*  
**Platform:** Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)  
**Original study design:** 104 breast-cancer biopsies and 17 normal breast-tissue samples.

GEO:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568

For a lightweight reproducible workflow, the notebook downloads a processed CSV copy archived on Zenodo:

**Zenodo DOI:** 10.5281/zenodo.4846212  
https://doi.org/10.5281/zenodo.4846212

The processed copy is used for analysis convenience; GSE42568 remains the biological source dataset.

## Why this project matters
Cancer changes which genes are more or less active in a cell. Gene-expression profiling measures those differences at large scale. This project asks whether breast-tumor tissue has a reproducible expression signature that distinguishes it from normal breast tissue.

## What the project does
1. Downloads a public processed copy of GSE42568.
2. Inspects the sample labels and expression matrix.
3. Performs quality checks.
4. Uses PCA to visualize major variation across samples.
5. Tests each expression feature for differences between tumor and normal tissue.
6. Corrects for thousands of simultaneous tests using Benjamini–Hochberg FDR.
7. Maps Affymetrix probe IDs to gene symbols using the official GPL570 annotation.
8. Generates a volcano plot and a heatmap of top signals.
9. Optionally performs pathway enrichment with Enrichr through `gseapy`.
10. Saves reproducible tables and figures.

## Methods

### Gene expression
A gene is a region of DNA. When a cell uses a gene, it produces RNA. Gene-expression assays measure how much RNA-associated signal is present for many genes at once.

### Differential expression
For every measured feature, we compare the tumor group with the normal group. A large difference alone is not enough; we also estimate how likely the difference is to appear by chance.

### Multiple-testing correction
Because thousands of features are tested simultaneously, some tiny p-values occur by chance. Benjamini–Hochberg correction controls the false discovery rate (FDR).

### PCA
Principal Component Analysis compresses thousands of expression measurements into a few axes that capture major sources of variation. It is useful for checking whether tumor and normal samples show broad separation and for spotting unusual samples.

### Pathway enrichment
A list of genes is often more interpretable as groups of related biological processes. Enrichment analysis asks whether particular pathways or biological functions occur more often than expected among the strongest genes.

## Important scientific limitation
GSE42568 is an **expression-microarray** dataset, not RNA-seq. It is still transcriptomics and is useful for learning differential-expression workflows, but the measurement technology differs from modern sequencing. This distinction should be stated accurately on a CV and in interviews.

## Suggested CV wording — only after you run and verify the results
**Breast Cancer Transcriptomic Signature Analysis — NCBI GEO GSE42568**  
*Python | GEO | Transcriptomics | Differential Expression | FDR | PCA | Pathway Enrichment*

- Analyzed a publicly available human breast-cancer transcriptomics cohort from NCBI GEO, comparing tumor and normal breast tissue through reproducible quality-control, PCA, and differential-expression workflows.
- Applied Welch's t-tests with Benjamini–Hochberg false-discovery-rate correction across high-dimensional expression features and mapped Affymetrix probes to gene-level annotations for biological interpretation.
- Visualized tumor-associated expression patterns with volcano plots and heatmaps and investigated enriched biological processes among significant genes.

Do **not** add exact counts of significant genes or specific pathways until the notebook has been run and those results have been verified.

## Repository structure
```text
breast_cancer_transcriptomics_GSE42568/
├── data/
│   ├── raw/
│   └── processed/
├── notebooks/
│   └── GSE42568_analysis.ipynb
├── scripts/
│   └── run_analysis.py
├── results/
│   ├── figures/
│   └── tables/
├── .gitignore
├── requirements.txt
└── README.md
```

## Run
```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python scripts/run_analysis.py
```

Or open `notebooks/GSE42568_analysis.ipynb` in Jupyter/VS Code/Colab.

## Data provenance
The project intentionally does not commit the downloaded expression matrix. The script retrieves it from the public archive so that the repository remains small and the source remains explicit.
