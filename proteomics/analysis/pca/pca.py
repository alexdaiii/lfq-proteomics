from pathlib import Path
from typing import NamedTuple

import pandas as pd
from matplotlib.figure import Figure
from sklearn.decomposition import PCA

from proteomics.analysis.io.load_metadata import MetadataMaps

from matplotlib import pyplot as plt
import seaborn as sns


class ScreeOutput(NamedTuple):
    figure: Figure
    explained_variance: float


def plot_scree(
    df: pd.DataFrame,
    *,
    plt_name: Path | None = None,
) -> ScreeOutput:
    pca = PCA()

    pca.fit(df.T)

    sns.set_style("ticks")

    # plot the scree plot
    fig = plt.figure()
    plt.plot(pca.explained_variance_ratio_)
    plt.xlabel("PC")
    plt.ylabel("Explained variance ratio")
    plt.title("Scree plot")

    explained_variance = pca.explained_variance_ratio_[:2].sum() * 100

    print(f"PC1 and PC2 explain {explained_variance:.2f}% of the variance")

    if plt_name:
        plt.savefig(plt_name)

    plt.close()

    return ScreeOutput(
        figure=fig,
        explained_variance=explained_variance,
    )


def plot_pca(
    df_pca: pd.DataFrame,
    fig_name: Path | None = None,
) -> Figure:
    # plot the PCA
    sns.set_style("ticks")

    fig, ax = plt.subplots()

    sns.scatterplot(
        data=df_pca,
        x="PC1",
        y="PC2",
        hue="Group",
        style="Group",
        s=100,
        ax=ax,
    )
    plt.title("PCA of the samples")

    # Add annotations with arrows
    for i in range(df_pca.shape[0]):
        pt_name = df_pca["Sample"].iloc[i].split("_imputed_intensity")[0]

        plt.annotate(
            pt_name,
            xy=(df_pca["PC1"].iloc[i], df_pca["PC2"].iloc[i]),
            xytext=(5, 5),
            textcoords="offset points",
            arrowprops=dict(arrowstyle="simple", lw=0.5, color="black"),
            fontsize=8,
        )

    sns.despine()

    if fig_name:
        plt.savefig(fig_name)

    plt.close()

    return fig


class PCAOutput(NamedTuple):
    pca_df: pd.DataFrame
    pca_fig: Figure
    scree_output: ScreeOutput


def run_pca(
    df_norm: pd.DataFrame,
    *,
    metadata_maps: MetadataMaps,
    scree_fig_name: Path,
    pca_fig_name: Path,
) -> PCAOutput:
    scree_output = plot_scree(df_norm, plt_name=scree_fig_name)

    pca = PCA(n_components=2)

    df_pca = pd.DataFrame(
        pca.fit_transform(df_norm.T),
        columns=["PC1", "PC2"],
    )

    df_col_to_condition = pd.DataFrame(
        [
            {"Sample": col, "Group": metadata_maps.sample_to_condition[col]}
            for col in df_norm.columns
        ]
    )

    df_pca = pd.concat(
        [
            df_pca,
            df_col_to_condition,
        ],
        axis=1,
    )

    pca_fig = plot_pca(df_pca, fig_name=pca_fig_name)

    return PCAOutput(
        pca_df=df_pca,
        pca_fig=pca_fig,
        scree_output=scree_output,
    )
