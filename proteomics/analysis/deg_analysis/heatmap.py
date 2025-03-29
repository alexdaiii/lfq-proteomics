from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import Colormap
from pydantic import BaseModel


def create_category_colors(
    *,
    metadata_df: pd.DataFrame,
    sample_col: str,
    category_col: str,
    category_name: str,
    category_colors: list[str | tuple],
    col_name_regex: str,
):
    lut = dict(zip(metadata_df[category_col].unique(), category_colors))
    colors = metadata_df[category_col].map(lut).to_frame(name=category_name)
    colors[sample_col] = metadata_df[sample_col]

    # regex to replace the column names (like if they have _imputed) after the sample
    colors[sample_col] = (
        colors[sample_col].astype(str).replace(col_name_regex, "", regex=True)
    )
    colors.set_index(sample_col, inplace=True)

    return colors


class MakeHeatmapOtherKwargs(BaseModel):
    # counts_df: pd.DataFrame
    # metadata_df: pd.DataFrame
    # sig_limma_results: pd.DataFrame
    col_name_replace_regex: str = ""
    sample_col: str = "Sample"
    """In the limma metadata input, what is the column name for the samples?"""
    category_col: str = "Group"
    """In the limma metadata input, what is the column name for the categories?"""
    category_name: str = "Condition"
    category_colors: list[str | tuple] | None = None
    cmap: str | list | Colormap | None = None
    sample_direction: Literal["row", "col"] = "col"
    show_gene_labels: bool = False
    show_sample_labels: bool = True
    use_z_score: bool = True
    cluster_samples: bool = True
    cluster_genes: bool = True
    title: str | None = None
    cbar_pos: tuple[float, float, float, float] | None = None
    cbar_units: str = "Z-score"
    cbar_label: str = "Abundance"
    figsize: tuple[float, float] | None = None
    dendrogram_ratio: tuple[float, float] | None = None

    class Config:
        arbitrary_types_allowed = True


def make_heatmap(
    *,
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    sig_limma_results: pd.DataFrame,
    col_name_replace_regex: str = "",
    sample_col: str,
    category_col: str,
    category_name: str = "Condition",
    category_colors: list[str | tuple],
    cmap: str | list = "coolwarm",
    sample_direction: Literal["row", "col"] = "col",
    show_gene_labels: bool = False,
    show_sample_labels: bool = True,
    use_z_score: bool = True,
    cluster_samples: bool = True,
    cluster_genes: bool = True,
    title: str = "",
    cbar_pos: tuple[float, float, float, float] = None,
    cbar_units: str = "Z-score",
    cbar_label: str = "Abundance",
    figsize: tuple[float, float] = None,
    dendrogram_ratio: tuple[float, float] = None,
) -> sns.matrix.ClusterGrid:
    """
    Creates a heatmap of the DE genes from a limma/deseq2 analysis.

    Args:
        counts_df: The counts dataframe. The row index MUST be the gene names.
        metadata_df: A metadata dataframe in deseq2 format.
        sig_limma_results: A dataframe of the significant results from limma.
            The index should be the genes that are significant.
        col_name_replace_regex: Should the column names be replaced with a regex?
        sample_col: In the metadata, what is the column name for the samples?
        category_col: In the metadata, what is the column name for the categories?
        category_name: What to label the categories on the heatmap. Default is 'Condition'.
        category_colors: A list of colors to use for the categories.
        cmap: The colormap to use for the heatmap. Default is 'coolwarm'.
        sample_direction: Should samples be plotted as rows or columns? Default is 'col'.
        show_gene_labels: Should the gene labels be shown? Default is False.
        show_sample_labels: Should the sample labels be shown? Default is True.
        use_z_score: Should the data be z-scored? Default is True.
        cluster_samples: Should the samples be clustered? Default is True.
        cluster_genes: Should the genes be clustered? Default is True.
        title: The title of the heatmap.
        cbar_pos: The position of the colorbar. Default is None.
        cbar_units: What to display as the units on the colorbar. Default is 'Z-score'.
        cbar_label: The label for the colorbar. Default is 'Abundance'.
        figsize: The size of the figure. Default is None.
        dendrogram_ratio: The ratio of the dendrogram to the heatmap. Default is None.

    Returns: The heatmap figure.
    """
    colors = create_category_colors(
        metadata_df=metadata_df,
        sample_col=sample_col,
        category_col=category_col,
        category_name=category_name,
        category_colors=category_colors,
        col_name_regex=col_name_replace_regex,
    )

    de_counts = counts_df.loc[sig_limma_results.index]
    de_counts.columns = de_counts.columns.str.replace(
        col_name_replace_regex, "", regex=True
    )

    cbar_kwargs = {"label": cbar_units}

    if sample_direction == "row":
        z_score = 1 if use_z_score else None
        figsize = (8, 4.5) if figsize is None else figsize
        dendrogram_ratio = (
            (0.05, 0.1) if dendrogram_ratio is None else dendrogram_ratio
        )
        cbar_pos = (0.3, 0.1, 0.3, 0.03) if cbar_pos is None else cbar_pos

        cbar_kwargs["orientation"] = "horizontal"

        fig = sns.clustermap(
            de_counts.T,
            cmap=cmap,
            z_score=z_score,
            col_cluster=cluster_genes,
            row_cluster=cluster_samples,
            figsize=figsize,
            yticklabels=show_sample_labels,
            xticklabels=show_gene_labels,
            row_colors=colors,
            dendrogram_ratio=dendrogram_ratio,
            cbar_pos=cbar_pos,
            cbar_kws=cbar_kwargs,
        )
    elif sample_direction == "col":
        z_score = 0 if use_z_score else None
        cbar_pos = (0.85, 0.5, 0.03, 0.3) if cbar_pos is None else cbar_pos
        dendrogram_ratio = (
            (0.1, 0.075) if dendrogram_ratio is None else dendrogram_ratio
        )
        figsize = (4.5, 6) if figsize is None else figsize

        fig = sns.clustermap(
            de_counts,
            cmap=cmap,
            z_score=z_score,
            col_cluster=cluster_samples,
            row_cluster=cluster_genes,
            figsize=figsize,
            yticklabels=show_gene_labels,
            xticklabels=show_sample_labels,
            dendrogram_ratio=dendrogram_ratio,
            cbar_pos=cbar_pos,
            col_colors=colors,
            cbar_kws=cbar_kwargs,
        )
    else:
        raise ValueError("Direction must be 'row' or 'col'")

    plt.title(cbar_label, fontweight="bold", loc="left", fontsize=10)
    fig.fig.suptitle(title, y=1.01, fontweight="bold")

    return fig


def make_heatmap_sample(
    *,
    counts_file: Path,
    metadata_file: Path,
    limma_results_file: Path,
    output_dir: Path,
    fc_cutoff: float = 1.5,
    pval_cutoff: float = 0.05,
    # kwargs for make_heatmap
    make_heatmap_kwargs: MakeHeatmapOtherKwargs,
) -> Path:
    """
    Creates the heatmap, returns the figure. Writes the figure, the DE matrix,
    and whatever seaborn plots to the output directory.

    Args:
        output_dir: The output directory.
        counts_file: The counts file.
        metadata_file: The metadata file.
        limma_results_file: The limma results file.
        fc_cutoff: The fold change cutoff.
        pval_cutoff: The p-value cutoff.
        make_heatmap_kwargs: The kwargs for :py:func:`make_heatmap`.

    Returns: The figure path.
    """

    if not output_dir.exists():
        raise FileNotFoundError(
            f"Output directory {output_dir} does not exist."
        )

    counts_df = pd.read_csv(counts_file, index_col=0)
    metadata_df = pd.read_csv(metadata_file)
    limma_results = pd.read_csv(limma_results_file, index_col=0)

    log_fc = np.log2(fc_cutoff)

    sig_limma_results = limma_results[
        (limma_results["logFC"].abs() >= log_fc)
        & (limma_results["adj.P.Val"] <= pval_cutoff)
    ]

    print(f"Found {len(sig_limma_results)} significant genes.")

    default_kwargs = {
        "category_colors": sns.color_palette(
            "Set1",
            n_colors=len(
                metadata_df[make_heatmap_kwargs.category_col].unique()
            ),
        ),
        "sample_col": "Sample",
        "category_col": "Group",
        "category_name": "Condition",
    }

    print("Creating heatmap...")
    cluster_fig = make_heatmap(
        counts_df=counts_df,
        metadata_df=metadata_df,
        sig_limma_results=sig_limma_results,
        **{
            **default_kwargs,
            **make_heatmap_kwargs.model_dump(exclude_unset=True),
        },
    )

    print("Saving heatmap...")
    for fmt in ["pdf", "png"]:
        cluster_fig.savefig(output_dir / f"heatmap.{fmt}", bbox_inches="tight")

    print("Saving DE matrix...")
    counts_df.loc[sig_limma_results.index].to_csv(output_dir / "de_matrix.csv")
    cluster_fig.data2d.to_csv(output_dir / "data2d.csv")

    plt.close()

    return output_dir / "heatmap.png"
