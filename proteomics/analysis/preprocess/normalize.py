from pathlib import Path
from typing import Literal, NamedTuple

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.figure import Figure


def center_median_normalize(df: pd.DataFrame) -> pd.DataFrame:
    center = df.median(axis=0)
    return df - center


def standardize(df: pd.DataFrame) -> pd.DataFrame:
    return (df - np.mean(df, axis=0)) / np.std(df, axis=0)


def log_normalize(
    df: pd.DataFrame, norm_method: Literal["center.median", "standardize"]
) -> pd.DataFrame:
    df_imputed_log2 = np.log2(df)

    if norm_method == "center.median":
        return center_median_normalize(df_imputed_log2)
    elif norm_method == "standardize":
        return pd.DataFrame(
            standardize(df_imputed_log2), columns=df_imputed_log2.columns
        )
    else:
        raise NotImplementedError(
            f"Normalization method {norm_method} not implemented"
        )


def plot_ms_abundances(
    df: pd.DataFrame,
    n_rows: int,
    n_cols: int,
    *,
    plt_name: Path | None = None,
    show: bool = False,
    title: str,
) -> Figure:
    if n_rows * n_cols < df.shape[1]:
        raise ValueError("n_rows * n_cols < number of columns in df")

    f1_fig, f1_ax = plt.subplots(
        n_rows, n_cols, figsize=(10, 10), constrained_layout=True
    )

    for i, ax in enumerate(f1_ax.flat):
        # if there are more subplots than columns in the dataframe, break
        if i >= df.shape[1]:
            break

        sns.histplot(df.iloc[:, i], ax=ax)
        ax.set_title(df.columns[i], fontsize=8)
        ax.set_xlabel("Intensity", fontsize=6)

    f1_fig.suptitle(title, fontsize=16)  # Set the title for the entire figure

    if plt_name:
        plt.savefig(plt_name)

    if show:
        plt.show()

    return f1_fig


class NormalizeRt(NamedTuple):
    df_norm: pd.DataFrame
    before_norm_fig: Figure
    after_norm_fig: Figure


def normalize(
    df: pd.DataFrame,
    *,
    plot_dir: Path,
    n_rows: int,
    n_cols: int,
) -> NormalizeRt:
    df_norm = log_normalize(df, "center.median")

    return NormalizeRt(
        df_norm=df_norm,
        before_norm_fig=plot_ms_abundances(
            df,
            n_rows,
            n_cols,
            plt_name=plot_dir / "before-normalization.png",
            title="Before normalization",
        ),
        after_norm_fig=plot_ms_abundances(
            df_norm,
            n_rows,
            n_cols,
            plt_name=plot_dir / "center-median.png",
            title="After normalization",
        ),
    )
