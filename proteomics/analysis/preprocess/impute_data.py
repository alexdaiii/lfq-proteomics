from typing import Literal

import numpy as np
import pandas as pd
from missforest import MissForest
import seaborn as sns

from proteomics.analysis.deg_analysis.heatmap import make_heatmap
from proteomics.analysis.io.load_metadata import MetadataMaps


__all__ = [
    "filter_data",
    "heatmap_missing_vals",
    "impute_missing_vals",
]

def filter_data(
        counts_org_df: pd.DataFrame,
        metadata_maps: MetadataMaps,
        *,
        min_count_per_group: int = 1,
        min_count_per_row: int | None = None,
        empty_type: Literal["NaN", "str"]
) -> pd.DataFrame:
    """
    Filters out genes that have too few observations
    before imputation. Each condition requires at
    least `min_count_per_group` observations. In
    addition, each row must have at least `min_count_per_row`
    observations. If `min_count_per_row` is None, then
    at least ceil(counts_df.shape[1] / 3) observations
    must not be empty.

    Args:
        counts_org_df: The counts dataframe. The row index MUST be the gene names.
        metadata_maps: A tuple with dictionaries that map the condition names
            to the group names.
        min_count_per_group: The minimum number of observations per group.
        min_count_per_row: The minimum number of observations per row.
        empty_type: The type of empty values in the counts_df.
            Can be "NaN" or "str".
    """
    counts = counts_org_df.copy()

    # replace str empty values with NaN
    if empty_type == "str":
        # display(counts.dtypes)
        print("Empty values are strings. Replacing with np.float64. "
              "Non numeric values are converted to NaN")
        counts = counts.apply(pd.to_numeric, errors='coerce')

        # display(counts.dtypes)
    initial_gene_count = counts.shape[0]

    # Filter genes based on minimum obs per group
    for condition, samples in metadata_maps.condition_to_sample.items():
        before_filter_count = counts.shape[0]
        counts = counts[counts[samples].notna().sum(axis=1) >= min_count_per_group]
        after_filter_count = counts.shape[0]
        print(f"Filtered out {before_filter_count - after_filter_count} "
              f"genes for condition '{condition}' that had less than "
              f"{min_count_per_group} observations")

    # Filter genes based on minimum obs per row
    total_samples = counts.shape[1]

    if min_count_per_row is None:
        min_count_per_row = int(np.ceil(total_samples / 3))
        print(f"Setting min_count_per_row to {min_count_per_row}")

    before_filter_count = counts.shape[0]
    counts = counts[counts.notna().sum(axis=1) >= min_count_per_row]
    after_filter_count = counts.shape[0]
    print(f"Filtered out {before_filter_count - after_filter_count} "
          f"genes based with less than {min_count_per_row} observations")

    final_gene_count = counts.shape[0]
    print(f"Total genes filtered out: {initial_gene_count - final_gene_count}")

    return counts

def heatmap_missing_vals(
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        *,
        sample_col: str,
        category_col: str,
        col_name_replace_regex: str = ""
):
    # Convert to boolean mask
    missing_mask = counts.isnull()
    # only get the rows with 1
    missing_mask = missing_mask[missing_mask.sum(axis=1) > 1]

    make_heatmap(
        counts_df=missing_mask,
        metadata_df=metadata,
        sig_limma_results=missing_mask,
        sample_col=sample_col,
        category_col=category_col,
        use_z_score=False,
        cmap=["black", "beige"],
        category_colors=sns.color_palette(
            "Set1",
        ),
        cbar_units="True/False",
        cbar_label="is_missing",
        col_name_replace_regex=col_name_replace_regex,
    )


def impute_missing_vals(
        counts: pd.DataFrame,

):
    # Initialize MissForest imputer
    mf = MissForest()

    # Fit and transform the entire dataset
    imputed_data = mf.fit_transform(counts)

    # Convert the imputed data back to a DataFrame
    imputed_df = pd.DataFrame(imputed_data,
                              columns=counts.columns,
                              index=counts.index)

    return imputed_df