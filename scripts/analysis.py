import scanpy as sc
import matplotlib.pyplot as plt

# Access Snakemake inputs/outputs/params
input_file = str(snakemake.input.h5ad)       # the filtered AnnData from preprocess step
umap_plot_file = str(snakemake.output.umap_plot)   # path to save the UMAP image
output_h5ad = str(snakemake.output.adata_final)    # path to save final AnnData
n_top_genes = int(snakemake.params.n_top_genes)    # number of HVGs from config

# 1. Load the filtered AnnData
adata = sc.read_h5ad(input_file)
print(f"Read AnnData for analysis: {adata.n_obs} cells, {adata.n_vars} genes")

adata = sc.read_h5ad(input_file)
print(f"Read AnnData for analysis: {adata.n_obs} cells, {adata.n_vars} genes")

adata.uns["log1p"] = {"base": None}

# 2. Identify highly variable genes (HVGs)
sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True, flavor='seurat')

# 2. Identify highly variable genes (HVGs)
# This focuses on genes with highest variance across cells, which will be used for PCA.
sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True, flavor='seurat')
# The above finds the top n_top_genes and subsets adata to only those genes (subset=True).
print(f"Selected {n_top_genes} highly variable genes for downstream analysis.")

# 3. Scale the data (optional but common before PCA).
# This step z-score normalizes each gene (so that each gene has mean~0 and variance~1).
sc.pp.scale(adata, max_value=10)
# (Clipping values to max 10 helps reduce the influence of extreme outliers.)

# 4. Principal Component Analysis (PCA)
sc.tl.pca(adata, svd_solver='arpack')
print("PCA complete. Explained variance of first 5 PCs:", adata.uns['pca']['variance_ratio'][:5])

# 5. Compute the k-nearest neighbors graph (for clustering/UMAP)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# 6. Clustering – use the Leiden algorithm (graph-based community detection)
sc.tl.leiden(adata, resolution=0.2)  # resolution controls how many clusters to find (higher -> more clusters)
# The resulting cluster labels are stored in adata.obs['leiden']
num_clusters = adata.obs['leiden'].nunique()
print(f"Leiden clustering found {num_clusters} clusters.")

# 7. UMAP dimensionality reduction
sc.tl.umap(adata)
# The UMAP coordinates are added to adata.obsm['X_umap']
print("Computed UMAP embedding for visualization.")

# 8. Plot UMAP with clusters
# Color points by the cluster assignment. We'll use Scanpy's plotting for convenience.
sc.pl.umap(adata, color='leiden', legend_loc='on data', title='PBMC 3k Leiden Clusters',
           show=False)  # don't show, as we're running in headless mode
plt.savefig(umap_plot_file)
plt.close()
print(f"UMAP plot saved to {umap_plot_file}")

# 9. Save the final AnnData with all results (PCA, UMAP, clusters, etc.)
adata.write(output_h5ad)
print(f"Final AnnData saved to {output_h5ad}")