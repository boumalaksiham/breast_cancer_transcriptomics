# Learning Guide

## 1. What exactly is this dataset?
Each sample is breast tissue from a person. The original GEO study contains breast tumors and normal breast tissue. Each expression feature corresponds to an Affymetrix probe designed to measure transcript abundance.

Think of the matrix like this:

| sample | class | probe_1 | probe_2 | ... |
|---|---|---:|---:|---|
| person A | normal | 7.1 | 4.3 | ... |
| person B | tumor | 9.6 | 3.8 | ... |

Rows are samples. Columns are expression measurements.

## 2. What does "expression" mean?
Expression is a proxy for how active a gene is. Higher RNA-associated signal generally means more transcript from that gene was present in the sample.

## 3. Why PCA before hypothesis testing?
PCA is a sanity check. We want to see the dominant structure in the data without telling the algorithm which samples are tumors. Strong group separation is informative; unexpected isolated points can reveal outliers or technical variation.

## 4. Why not just compare means?
With tens of thousands of measurements, random differences are guaranteed. We therefore combine:
- an **effect size**: how large the tumor-normal difference is;
- a **p-value**: how compatible the observed difference is with random variation under the null;
- an **FDR-adjusted p-value**: correction for the huge number of tests.

## 5. Why Welch's t-test?
The groups have different sample sizes and may have different variances. Welch's test does not require equal variance.

## 6. What does FDR < 0.05 mean?
Roughly, among results called significant under that threshold, the procedure aims to limit the expected proportion of false discoveries to about 5%, under its assumptions. It does **not** mean there is a 95% probability the gene is truly involved in cancer.

## 7. Why map probes to genes?
Microarrays measure probe sets. Biological interpretation is usually done at the gene level, so we use the GPL570 platform annotation to connect probe identifiers to gene symbols.

## 8. Why pathway enrichment?
A single-gene list is hard to interpret. Pathway enrichment asks whether genes involved in the same biological process appear more often than expected.

## 9. What can you legitimately say on a CV?
You can say you analyzed a public breast-cancer transcriptomics dataset, performed PCA, differential-expression testing, FDR correction, probe-to-gene annotation, visualization, and pathway analysis **if you actually run those steps**.

You should not call the data RNA-seq. GSE42568 is microarray transcriptomics.

## 10. What makes this useful for biomedical-informatics outreach?
It shows that you can move from:
public biological data → preprocessing/QC → statistical testing → biological annotation → interpretation → reproducible documentation.

That is a real bioinformatics workflow, even though it is intentionally scoped as a first project.
