import scanpy as sc

# Access Snakemake parameters (from Snakefile)
input_dir = str(snakemake.input[0])            # the "data/{sample}/" directory
output_file = str(snakemake.output[0])         # path for output AnnData .h5ad
min_genes = int(snakemake.params.min_genes)    # cell filter threshold
min_cells = int(snakemake.params.min_cells)    # gene filter threshold

# 1. Read 10x Genomics data
# Scanpy can read 10x formatted data (mtx + tsv files) using sc.read_10x_mtx.
# We specify the directory containing matrix.mtx, features.tsv, barcodes.tsv.
print(f"Reading 10x data from {input_dir} ...")
adata = sc.read_10x_mtx(input_dir, var_names='gene_symbols', cache=True)  
# var_names='gene_symbols' will use gene names (if features.tsv present) instead of gene IDs.
# The AnnData 'adata' now contains the raw count matrix.
print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes.")

# 2. Filter cells and genes based on thresholds
# Filter out cells with fewer than min_genes expressed genes
sc.pp.filter_cells(adata, min_genes=min_genes)
# Filter out genes expressed in fewer than min_cells cells
sc.pp.filter_genes(adata, min_cells=min_cells)
print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes remain.")

# (Optional) You could add additional QC filtering here, e.g., remove cells with high mitochondrial gene percent, etc.
# But for simplicity, we only apply the above basic filters.

# 3. Normalize the data
# Normalize each cell’s total counts to the same target (e.g., 1e4) and log-transform.
sc.pp.normalize_total(adata, target_sum=1e4)   # library-size normalize each cell to 10,000 counts [oai_citation:12‡training.galaxyproject.org](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html#:~:text=match%20at%20L998%20,for%20the%20computation%20of%20the)
sc.pp.log1p(adata)  # log-transform the data (logarithm of (counts+1))

# You may store the raw counts if needed for later use (e.g., for finding DE genes after clustering)
# adata.raw = adata  # save raw (log-transformed) data for later, if desired

# 4. Save the preprocessed AnnData to file
adata.write(output_file)
print(f"Saved filtered & normalized data to {output_file}")