import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, BoundaryNorm
import scanpy as sc
from scipy.stats import spearmanr
from io import StringIO

def plot_expression_diff_vs_exon_count(nano, ill, exon_count_file, output_prefix):
    # Load exon count data
    exon_df = pd.read_csv(exon_count_file, sep='\t', header=None, names=["transcript_id", "exon_count"])
    exon_df["exon_count"] = pd.to_numeric(exon_df["exon_count"], errors="coerce")
    exon_df.dropna(subset=["exon_count"], inplace=True)
    exon_df["exon_count"] = exon_df["exon_count"].astype(int)

    # Calculate mean expression across cells for each transcript
    nano_mean = pd.Series(np.asarray(nano.X.mean(axis=0)).flatten(), index=nano.var.index)
    ill_mean = pd.Series(np.asarray(ill.X.mean(axis=0)).flatten(), index=ill.var.index)

    # Combine into a DataFrame
    df = pd.DataFrame({
        "nano": nano_mean,
        "ill": ill_mean
    })

    # Calculate log2 expression difference
    df["expression_diff"] = np.log2(df["ill"] + 1) - np.log2(df["nano"] + 1)

    # Merge with exon count
    df = df.merge(exon_df, left_index=True, right_on="transcript_id")

    # Group by exon count and compute mean and std dev
    grouped = df.groupby("exon_count")["expression_diff"].agg(["mean", "std"]).reset_index()

    # Sort by exon count
    grouped = grouped.sort_values("exon_count")

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(grouped["exon_count"], grouped["mean"], label="Mean log2 diff", color="darkblue")
    plt.fill_between(grouped["exon_count"],
                     grouped["mean"] - grouped["std"],
                     grouped["mean"] + grouped["std"],
                     color="lightblue", alpha=0.5, label="±1 SD")

    # Set x-ticks: only label 20 evenly spaced values
    xtick_locs = np.linspace(grouped["exon_count"].min(), grouped["exon_count"].max(), 20, dtype=int)
    xtick_labels = [str(x) for x in xtick_locs]
    plt.xticks(ticks=xtick_locs, labels=xtick_labels, rotation=45)

    plt.xlabel("Exon Count")
    plt.ylabel("Log2 Expression Difference (Illumina - ONT)")
    plt.title("Expression Difference vs Exon Count")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.png", dpi=300)
    plt.close()
    print(f"Saved plot: {output_prefix}.png")


def analyze_single_cell_expression_correlation(nano_dir, ill_dir, top_n=100, output_prefix='correlation_plot', plot_range=(0,13), priming='both'):
    """
    Compare single-cell expression profiles between Nanopore and Illumina data using top N transcripts.

    Parameters:
        nano_dir (str): Path to Nanopore data directory.
        ill_dir (str): Path to Illumina data directory.
        top_n (int): Number of top-expressed transcripts to use.
        output_prefix (str): Prefix for saved plots.
        plot_range (tuple): Axis limits for the scatter plots.
    """
    def load_data(directory):
        matrix = sc.read_mtx(f"{directory}/matrix.abundance.mtx")
        matrix.obs.index = pd.read_csv(f"{directory}/count.barcodes.txt", header=None)[0].values
        with open(f"{directory}/transcripts.txt", 'r') as f:
            data = f.read().replace('\n', ',').replace(',,', ',')
        matrix.var.index = [str(t) for t in list(pd.read_csv(StringIO(data), header=None, sep=",").values[0]) if 'ENSM' in str(t)]
        return matrix
    
    '''
    def plot_expression_diff_vs_exon_count(nano, ill, exon_count_file, output_prefix='diff_vs_exons'):
        """
        Plot difference in expression (Nanopore - Illumina) vs exon count.

        Parameters:
            nano (AnnData): Nanopore expression data.
            ill (AnnData): Illumina expression data.
            exon_count_file (str): TSV with transcript_id and exon_count.
            output_prefix (str): Prefix for the saved plot.
        """
        # Load exon count file
        exon_df = pd.read_csv(exon_count_file, sep='\t', header=None, names=['transcript_id', 'exon_count'])
        exon_df = exon_df.set_index('transcript_id')

        # Match transcripts in both datasets
        common_genes = nano.var_names.intersection(ill.var_names).intersection(exon_df.index)
        if len(common_genes) == 0:
            raise ValueError("No overlapping transcripts found among Nanopore, Illumina, and exon count file.")

        # Extract mean expression per transcript
        nano_expr = np.asarray(nano[:, common_genes].X.mean(axis=0)).flatten()
        ill_expr = np.asarray(ill[:, common_genes].X.mean(axis=0)).flatten()
        diff_expr = nano_expr - ill_expr  # Nanopore - Illumina

        # Build DataFrame
        df = pd.DataFrame({
            'transcript_id': common_genes,
            'nano_mean': nano_expr,
            'ill_mean': ill_expr,
            'expression_diff': diff_expr,
            'exon_count': exon_df.loc[common_genes, 'exon_count'].values
        })

        df['exon_count'] = pd.to_numeric(df['exon_count'], errors='coerce')
        df = df.dropna(subset=['exon_count'])  # Drop any transcripts with missing/invalid exon count
        df['exon_count'] = df['exon_count'].astype(int)
        print(df)

        # Scatter Plot: Expression Diff vs Exon Count
        plt.figure(figsize=(8, 6))
        plt.scatter(df['exon_count'], df['expression_diff'], alpha=0.6, s=10)
        plt.axhline(0, linestyle='--', color='gray')
        plt.xlabel("Exon count")
        plt.ylabel("Expression difference (Nanopore - Illumina)")
        plt.title("Expression difference vs Exon count per transcript")
    
        # Get unique exon counts
        unique_exons = np.sort(df['exon_count'].unique())

        # Show only 20 x-tick labels
        if len(unique_exons) > 20:
            xtick_vals = np.linspace(unique_exons.min(), unique_exons.max(), 20, dtype=int)
            xtick_labels = [str(x) for x in xtick_vals]
            plt.xticks(xtick_vals, xtick_labels, rotation=180)
        else:
            plt.xticks(unique_exons, [str(x) for x in unique_exons], rotation=180)
            plt.tight_layout()
            plt.savefig(f"{output_prefix}.png", dpi=600)
            plt.close()

        print(f"Saved plot: {output_prefix}.png")
        return df  # return dataframe in case further analysis is needed
    '''

    # Load and preprocess
    nano = load_data(nano_dir)
    ill = load_data(ill_dir)
    sc.pp.filter_cells(nano, min_genes=100)
    sc.pp.filter_cells(nano, min_counts=500)
    sc.pp.filter_genes(nano, min_cells=3)
    sc.pp.filter_cells(ill, min_genes=100)
    sc.pp.filter_cells(ill, min_counts=500)
    sc.pp.filter_genes(ill, min_cells=3)
    sc.pp.normalize_total(nano, target_sum=1000000)
    sc.pp.normalize_total(ill, target_sum=1000000)

    plot_expression_diff_vs_exon_count(nano, ill, exon_count_file='mm39.exons_count.tsv', output_prefix=f'{output_prefix}_exp_diff_vs_exons')


    # Find common cells
    common_cells = nano.obs.index.intersection(ill.obs.index)
    if len(common_cells) == 0:
        raise ValueError("No common cells found between datasets.")

    nano_common = nano[common_cells, :]
    ill_common = ill[common_cells, :]

    # Get top expressed genes (from Nanopore only)
    mean_expr_nano = np.asarray(nano_common.X.mean(axis=0)).flatten()
    top_gene_indices = np.argsort(mean_expr_nano)[-top_n:]
    top_genes = nano.var.index[top_gene_indices]

    # Compute Spearman correlations
    spearman_corrs = []
    for cell in common_cells:
        x = nano[cell, top_genes].X.toarray().flatten()
        y = ill[cell, top_genes].X.toarray().flatten()
        corr, _ = spearmanr(x, y)
        spearman_corrs.append(np.nan_to_num(corr))
    spearman_corrs = np.array(spearman_corrs)

    # Compute mean expression per cell across top genes
    nano_mean = np.asarray(nano_common[:, top_genes].X.mean(axis=1)).flatten()
    ill_mean = np.asarray(ill_common[:, top_genes].X.mean(axis=1)).flatten()

    # Compute median expression
    nano_median = np.median(nano_common[:, top_genes].X.toarray(), axis=1)
    ill_median = np.median(ill_common[:, top_genes].X.toarray(), axis=1)

    # log2(TPM + 1)
    nano_mean_log = np.log2(nano_mean + 1)
    ill_mean_log = np.log2(ill_mean + 1)
    nano_median_log = np.log2(nano_median + 1)
    ill_median_log = np.log2(ill_median + 1)

    # === Discretize the colormap ===
    n_bins = 4  # Number of discrete colors
    cmap = cm.get_cmap('jet', n_bins)  # Get discretized version of 'jet'
    
    # Normalize your correlation values
    corr_min = 0.0 #np.min(spearman_corrs)
    corr_max = 1.0 #np.max(spearman_corrs)
    bounds = np.linspace(corr_min, corr_max, n_bins + 1)
    norm = BoundaryNorm(boundaries=bounds, ncolors=n_bins)
    
    # Compute bin counts and percentages
    counts, _ = np.histogram(spearman_corrs, bins=bounds)
    percentages = (counts / len(spearman_corrs)) * 100
    
    # Print out the bin ranges and percentages
    for i in range(n_bins):
        print(f"Bin {i+1}: {bounds[i]:.2f} to {bounds[i+1]:.2f} → {percentages[i]:.1f}%")

    # === Mean Plot ===
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(nano_mean_log, ill_mean_log, c=spearman_corrs, cmap=cmap, norm=norm, s=5)
    cbar = plt.colorbar(scatter, boundaries=bounds, ticks=bounds)
    cbar.set_label('Spearman correlation')
    
    # Add percentages as labels next to the ticks
    tick_labels = [f"{bounds[i]:.2f}–{bounds[i+1]:.2f}\n{percentages[i]:.1f}%" for i in range(n_bins)]
    cbar.set_ticks((bounds[:-1] + bounds[1:]) / 2)
    cbar.set_ticklabels(tick_labels)

    plt.xlabel("Nanopore log2(mean TPM + 1)")
    plt.ylabel("Illumina log2(mean TPM + 1)")
    plt.xlim(*plot_range)
    plt.ylim(*plot_range)
    plt.title(f"Mean expression per {priming}cell (Top {top_n} Nanopore transcripts)")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_mean.png", dpi=600)
    plt.close()
    
    # === Median Plot ===
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(nano_median_log, ill_median_log, c=spearman_corrs, cmap=cmap, norm=norm, s=5)
    cbar = plt.colorbar(scatter, boundaries=bounds, ticks=bounds)
    cbar.set_label('Spearman correlation')

    cbar.set_ticks((bounds[:-1] + bounds[1:]) / 2)
    cbar.set_ticklabels(tick_labels)

    plt.xlabel("Nanopore log2(median TPM + 1)")
    plt.ylabel("Illumina log2(median TPM + 1)")
    plt.xlim(*plot_range)
    plt.ylim(*plot_range)
    plt.title(f"Median expression per {priming}cell (Top {top_n} Nanopore transcripts)")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_median.png", dpi=600)
    plt.close()

    print(f"Plots saved as '{output_prefix}_mean.png' and '{output_prefix}_median.png'.")


analyze_single_cell_expression_correlation(
    nano_dir="13G_single_cell/b01_nanopore_13G_randO_single_cell",
    ill_dir="b01_next1_13G_single_cell_randO",
    top_n=100,
    output_prefix="13G_correlation_top100_randO",
    plot_range=(0, 7),
    priming='randO'
)

analyze_single_cell_expression_correlation(
    nano_dir="13G_single_cell/b01_nanopore_13G_polyT_single_cell",
    ill_dir="b01_next1_13G_single_cell_polyT",
    top_n=100,
    output_prefix="13G_correlation_top100_polyT",
    plot_range=(0, 7), 
    priming='polyT'
)

analyze_single_cell_expression_correlation(
    nano_dir="13G_single_cell/b01_nanopore_13G_single_cell",
    ill_dir="b01_next1_13G_single_cell",
    top_n=100,
    output_prefix="13G_correlation_top100_both",
    plot_range=(0, 7),
    priming=''
)




