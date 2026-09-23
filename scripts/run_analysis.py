from pathlib import Path
import gzip
import io
import re
import requests
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"

for d in (RAW, PROCESSED, FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

DATA_URL = "https://zenodo.org/records/4846212/files/GSE42568.csv?download=1"

def download(url: str, destination: Path):
    if destination.exists() and destination.stat().st_size > 1000:
        print(f"Using existing {destination.name}")
        return
    print(f"Downloading {destination.name}...")
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    destination.write_bytes(r.content)

def infer_label_column(df):
    candidates = []
    for c in df.columns:
        vals = df[c].dropna().astype(str).str.lower().unique()
        if 2 <= len(vals) <= 10:
            joined = " ".join(vals)
            if any(k in joined for k in ["normal", "tumor", "tumoral", "cancer", "control"]):
                candidates.append(c)
    if not candidates:
        raise ValueError(
            "Could not infer the tumor/normal label column. "
            "Print df.head() and df.columns, then set LABEL_COLUMN manually."
        )
    return candidates[0]

def normalize_label(x):
    s = str(x).strip().lower()
    if "normal" in s or "control" in s:
        return "Normal"
    if "tumor" in s or "tumoral" in s or "cancer" in s:
        return "Tumor"
    return str(x)

def benjamini_hochberg(pvals):
    return multipletests(pvals, method="fdr_bh")[1]

def map_entrez_to_symbols(entrez_ids):
    """
    Map Entrez Gene IDs to human gene symbols using MyGene.info.
    """
    entrez_ids = [str(x) for x in entrez_ids]

    results = []

    batch_size = 500

    for i in range(0, len(entrez_ids), batch_size):
        batch = entrez_ids[i:i + batch_size]

        response = requests.post(
            "https://mygene.info/v3/query",
            data={
                "q": ",".join(batch),
                "scopes": "entrezgene",
                "fields": "symbol",
                "species": "human"
            },
            timeout=60
        )

        response.raise_for_status()
        results.extend(response.json())

    mapping = {}

    for item in results:
        query = str(item.get("query"))

        if "symbol" in item:
            mapping[query] = item["symbol"]

    return mapping

def main():
    data_path = RAW / "GSE42568.csv"
    download(DATA_URL, data_path)

    df = pd.read_csv(data_path)
    print("Raw processed-mirror shape:", df.shape)
    print("First columns:", df.columns[:8].tolist())

    label_col = infer_label_column(df)
    print("Detected label column:", label_col)
    labels = df[label_col].map(normalize_label)

    # Keep only tumor/normal rows.
    keep = labels.isin(["Tumor", "Normal"])
    df = df.loc[keep].reset_index(drop=True)
    labels = labels.loc[keep].reset_index(drop=True)

    # Determine sample identifier column if present.
    non_numeric = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]
    id_candidates = [c for c in non_numeric if c != label_col]
    sample_id_col = id_candidates[0] if id_candidates else None

    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    X = df[numeric_cols].apply(pd.to_numeric, errors="coerce")

    # Remove unusable features.
    X = X.loc[:, X.notna().mean() >= 0.95]
    X = X.fillna(X.median())
    X = X.loc[:, X.var() > 0]

    print("Samples retained:", len(X))
    print(labels.value_counts())
    print("Expression features retained:", X.shape[1])

    # Save basic sample table.
    sample_table = pd.DataFrame({
        "sample_id": df[sample_id_col].astype(str) if sample_id_col else [f"sample_{i+1}" for i in range(len(df))],
        "group": labels
    })
    sample_table.to_csv(PROCESSED / "sample_metadata.csv", index=False)

    # QC summary.
    qc = pd.DataFrame({
        "sample_id": sample_table["sample_id"],
        "group": labels,
        "mean_expression": X.mean(axis=1),
        "median_expression": X.median(axis=1),
        "sd_expression": X.std(axis=1)
    })
    qc.to_csv(TAB / "qc_summary.csv", index=False)

    # PCA on the most variable features to reduce noise/computation.
    n_top = min(5000, X.shape[1])
    top_var = X.var().nlargest(n_top).index
    Z = StandardScaler().fit_transform(X[top_var])
    pca = PCA(n_components=2, random_state=42)
    pcs = pca.fit_transform(Z)
    pca_df = pd.DataFrame({
        "PC1": pcs[:, 0],
        "PC2": pcs[:, 1],
        "group": labels,
        "sample_id": sample_table["sample_id"]
    })
    pca_df.to_csv(TAB / "pca_coordinates.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 5))
    for group in ["Normal", "Tumor"]:
        m = pca_df["group"] == group
        ax.scatter(pca_df.loc[m, "PC1"], pca_df.loc[m, "PC2"], label=group, alpha=0.75)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
    ax.set_title("PCA of GSE42568 expression profiles")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "pca_tumor_vs_normal.png", dpi=200)
    plt.close(fig)

    # Differential expression at probe/feature level.
    tumor = X.loc[labels.eq("Tumor").values]
    normal = X.loc[labels.eq("Normal").values]

    rows = []
    for feature in X.columns:
        tvals = tumor[feature].to_numpy(dtype=float)
        nvals = normal[feature].to_numpy(dtype=float)
        stat, p = ttest_ind(tvals, nvals, equal_var=False, nan_policy="omit")
        rows.append({
            "probe_id": str(feature),
            "mean_tumor": np.nanmean(tvals),
            "mean_normal": np.nanmean(nvals),
            "mean_difference_tumor_minus_normal": np.nanmean(tvals) - np.nanmean(nvals),
            "t_statistic": stat,
            "p_value": p
        })

    de = pd.DataFrame(rows)
    de["fdr_bh"] = benjamini_hochberg(de["p_value"].fillna(1).values)
    # The processed Zenodo dataset already uses Entrez Gene IDs.
    de["entrez_id"] = (
    de["probe_id"]
    .astype(str)
    .str.replace("gene_", "", regex=False)
    )
    print("Mapping Entrez IDs to gene symbols...")
    gene_mapping = map_entrez_to_symbols(
    de["entrez_id"].unique()
    )
    de["gene_symbol"] = de["entrez_id"].map(gene_mapping)
    print(
    "Gene symbols mapped:",
    de["gene_symbol"].notna().sum(),
    "of",
    len(de)
    )    

    de["minus_log10_fdr"] = -np.log10(de["fdr_bh"].clip(lower=1e-300))
    de = de.sort_values(["fdr_bh", "p_value"])
    de.to_csv(TAB / "differential_expression_all_probes.csv", index=False)

    # A transparent significance rule for exploratory analysis.
    sig = de[(de["fdr_bh"] < 0.05) & (de["mean_difference_tumor_minus_normal"].abs() >= 1.0)].copy()
    sig.to_csv(TAB / "significant_probes_fdr05_absdiff1.csv", index=False)

    # Collapse to one strongest probe per gene for interpretation.
    gene_de = (
        de.dropna(subset=["gene_symbol"])
          .assign(abs_diff=lambda d: d["mean_difference_tumor_minus_normal"].abs())
          .sort_values(["gene_symbol", "fdr_bh", "abs_diff"], ascending=[True, True, False])
          .drop_duplicates("gene_symbol")
          .sort_values(["fdr_bh", "abs_diff"], ascending=[True, False])
    )
    gene_de.to_csv(TAB / "differential_expression_gene_level.csv", index=False)

    # Volcano plot.
    fig, ax = plt.subplots(figsize=(8, 6))
    is_sig = (de["fdr_bh"] < 0.05) & (de["mean_difference_tumor_minus_normal"].abs() >= 1.0)
    ax.scatter(
        de.loc[~is_sig, "mean_difference_tumor_minus_normal"],
        de.loc[~is_sig, "minus_log10_fdr"],
        s=8, alpha=0.35
    )
    ax.scatter(
        de.loc[is_sig, "mean_difference_tumor_minus_normal"],
        de.loc[is_sig, "minus_log10_fdr"],
        s=10, alpha=0.7
    )
    ax.axvline(-1, linestyle="--", linewidth=1)
    ax.axvline(1, linestyle="--", linewidth=1)
    ax.axhline(-np.log10(0.05), linestyle="--", linewidth=1)
    ax.set_xlabel("Mean expression difference (Tumor − Normal)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("Differential expression: breast tumor vs normal tissue")
    fig.tight_layout()
    fig.savefig(FIG / "volcano_plot.png", dpi=200)
    plt.close(fig)

    # Heatmap of top genes using matplotlib only.
    top_genes = gene_de.head(25)["probe_id"].tolist()
    heat = X[top_genes].copy()
    heat = (heat - heat.mean(axis=0)) / heat.std(axis=0).replace(0, 1)
    order = np.argsort(labels.map({"Normal":0, "Tumor":1}).values)
    heat = heat.iloc[order]
    ordered_labels = labels.iloc[order].tolist()
    gene_names = (
        gene_de.set_index("probe_id").loc[top_genes, "gene_symbol"]
        .fillna(pd.Series(top_genes, index=top_genes))
        .tolist()
    )

    fig, ax = plt.subplots(figsize=(11, 7))
    im = ax.imshow(heat.T.values, aspect="auto", interpolation="nearest")
    ax.set_yticks(range(len(gene_names)))
    ax.set_yticklabels(gene_names, fontsize=7)
    ax.set_xlabel("Samples (Normal first, then Tumor)")
    ax.set_title("Top differential-expression signals (z-scored by feature)")
    fig.colorbar(im, ax=ax, label="Standardized expression")
    fig.tight_layout()
    fig.savefig(FIG / "top_gene_heatmap.png", dpi=200)
    plt.close(fig)

    # Optional pathway enrichment using the public Enrichr API.
    # This avoids requiring gseapy/Rust and keeps the project lightweight.
    try:
        gene_list = gene_de.loc[
            (gene_de["fdr_bh"] < 0.05) &
            (gene_de["mean_difference_tumor_minus_normal"].abs() >= 1.0),
            "gene_symbol"
        ].dropna().astype(str).drop_duplicates().head(500).tolist()

        if len(gene_list) >= 10:
            add_url = "https://maayanlab.cloud/Enrichr/addList"
            payload = {
                "list": (None, "\n".join(gene_list)),
                "description": (None, "GSE42568 tumor vs normal significant genes")
            }
            add_resp = requests.post(add_url, files=payload, timeout=60)
            add_resp.raise_for_status()
            user_list_id = add_resp.json()["userListId"]

            libraries = ["GO_Biological_Process_2023", "KEGG_2021_Human"]
            enrichment_frames = []

            for library in libraries:
                enrich_url = "https://maayanlab.cloud/Enrichr/enrich"
                params = {
                    "userListId": user_list_id,
                    "backgroundType": library
                }
                enrich_resp = requests.get(enrich_url, params=params, timeout=60)
                enrich_resp.raise_for_status()
                rows = enrich_resp.json().get(library, [])

                if rows:
                    frame = pd.DataFrame(
                        rows,
                        columns=[
                            "rank",
                            "term",
                            "p_value",
                            "z_score",
                            "combined_score",
                            "overlapping_genes",
                            "adjusted_p_value",
                            "old_p_value",
                            "old_adjusted_p_value"
                        ]
                    )
                    frame.insert(0, "library", library)
                    enrichment_frames.append(frame)

            if enrichment_frames:
                enrichment = pd.concat(enrichment_frames, ignore_index=True)
                enrichment.to_csv(TAB / "pathway_enrichment.csv", index=False)
                print("Pathway enrichment saved.")
            else:
                print("Enrichr returned no enrichment rows.")
        else:
            print("Not enough significant mapped genes for enrichment.")
    except Exception as e:
        print("Pathway enrichment skipped (analysis still completed):", e)

    summary = {
        "samples_analyzed": int(len(X)),
        "tumor_samples": int((labels == "Tumor").sum()),
        "normal_samples": int((labels == "Normal").sum()),
        "features_tested": int(X.shape[1]),
        "significant_probes_fdr05_absdiff1": int(len(sig)),
        "mapped_genes": int(gene_de["gene_symbol"].notna().sum()),
        "pca_variance_pc1": float(pca.explained_variance_ratio_[0]),
        "pca_variance_pc2": float(pca.explained_variance_ratio_[1])
    }
    pd.Series(summary).to_csv(TAB / "analysis_summary.csv", header=["value"])
    print("\nAnalysis complete:")
    for k, v in summary.items():
        print(f"  {k}: {v}")

if __name__ == "__main__":
    main()
