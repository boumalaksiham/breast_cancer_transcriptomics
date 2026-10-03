# Raw data

`GSE42568.csv` is included in this repository. The analysis script uses this cached file and downloads it only when it is missing.

Source used by the script: [GSE42568.csv](https://zenodo.org/records/4846212/files/GSE42568.csv?download=1).

`GPL570.annot.gz` is also tracked as an annotation artifact. The current Python analysis uses the MyGene service for mapping; it does not download or consume this GPL archive.

Run `python scripts/run_analysis.py` from the repository root, after installing its requirements. Preserve the input file and record its checksum when comparing runs; cached data and a later remote download may not be identical.
