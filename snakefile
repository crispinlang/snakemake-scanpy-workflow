# Snakefile

# Load the configuration file
configfile: "config.yaml"

# Get list of samples from config
SAMPLES = config["samples"]

# Define the final outputs for all samples. 
# Here we expect a UMAP plot and a processed AnnData file for each sample.
rule all:
    input:
        # Collect outputs for each sample using list comprehension
        expand("results/{sample}/umap_{sample}.png", sample=SAMPLES),
        expand("results/{sample}/adata_{sample}.h5ad", sample=SAMPLES)
    # No shell or script for the 'all' rule; it's just a summary of final targets.

# Rule: Preprocess data (load raw 10x, filter cells/genes, normalize).
# This rule takes the raw data files and produces a filtered+normalized AnnData.
rule preprocess:
    # Input: path to the 10x Genomics matrix directory for this sample.
    # We assume the 10x data (matrix.mtx, features.tsv, barcodes.tsv) are in data/{sample}/
    input:
        "data/{sample}/"   # directory with 10x data for the sample
    output:
        # Save intermediate AnnData after filtering & normalization
        "results/{sample}/adata_{sample}_filtered.h5ad"
    params:
        # Pass filtering parameters from config to the script
        min_genes=config["min_genes"],
        min_cells=config["min_cells"]
    threads: 1
    conda:
        "envs/scanpy.yaml"   # Use Scanpy conda environment for this rule
    script:
        "scripts/preprocess.py"   # This script will read input and apply preprocessing

# Rule: Analysis (HVGs, PCA, neighbors, clustering, UMAP, plotting).
# This uses the preprocessed AnnData to do PCA, find neighbors, cluster, and generate UMAP.
rule analyze:
    input:
        # Input is the filtered AnnData from the previous step
        h5ad="results/{sample}/adata_{sample}_filtered.h5ad"
    output:
        # Final outputs: (1) UMAP plot image, (2) final AnnData with all results
        umap_plot="results/{sample}/umap_{sample}.png",
        adata_final="results/{sample}/adata_{sample}.h5ad"
    params:
        # Pass HVG and other parameters from config
        n_top_genes=config["n_top_genes"]
    threads: 1
    conda:
        "envs/scanpy.yaml"   # Same environment (Scanpy) for this analysis step
    script:
        "scripts/analysis.py"  # This script performs PCA, clustering, UMAP, etc.