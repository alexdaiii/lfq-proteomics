from pathlib import Path
from typing import TypedDict, NamedTuple

import pandas as pd


class ContrastInput(TypedDict):
    contrast: tuple[str, str]
    genes_sheet: str | None
    gene_list_col: str | None


class Contrast(ContrastInput):
    group_0_cols: list[str]
    group_1_cols: list[str]


class MetadataMaps(NamedTuple):
    condition_to_sample: dict[str, list[str]]
    sample_to_condition: dict[str, str]


def make_metadata(
    metadata_file: Path,
    metadata_sheet_name: str,
    metadata_index_col: int,
    condition_col: str,
) -> MetadataMaps:
    """
    Maps the conditions to sample names
    """
    df = pd.read_excel(
        metadata_file,
        sheet_name=metadata_sheet_name,
        index_col=metadata_index_col,
    )

    # make sure condition col exists and is not null
    if condition_col not in df.columns:
        raise ValueError(
            f"Condition column {condition_col} not found in metadata"
        )

    if df[condition_col].isnull().values.any():
        raise ValueError(
            f"Condition column {condition_col} has missing values"
        )

    condition_to_sample = (
        df.groupby(condition_col).apply(lambda x: x.index.tolist(), include_groups=False).to_dict()
    )
    sample_to_condition = df[condition_col].to_dict()

    return MetadataMaps(condition_to_sample, sample_to_condition)


def validate_metadata(
    counts_df: pd.DataFrame,
    metadata_maps: MetadataMaps,
) -> None:
    """
    Makes sure that all samples (columns) in the counts_df are in the metadata
    """

    counts_samples = set(counts_df.columns.tolist())
    metadata_samples = set(metadata_maps.sample_to_condition.keys())

    if counts_samples != metadata_samples:
        print(counts_samples)
        print(metadata_samples)
        raise ValueError(
            "Samples in the counts data do not match samples in the metadata"
        )

    print("Metadata validated")


def create_contrast_from_metadata(
    metadata_maps: MetadataMaps,
    contrasts: list[ContrastInput],
) -> list[Contrast]:
    contrast_list = []

    for contrast in contrasts:
        # check to make sure the conditions exist in the metadata
        if contrast["contrast"][0] not in metadata_maps.condition_to_sample:
            raise ValueError(
                f"Condition {contrast['contrast'][0]} not found in metadata"
            )
        if contrast["contrast"][1] not in metadata_maps.condition_to_sample:
            raise ValueError(
                f"Condition {contrast['contrast'][1]} not found in metadata"
            )

        group_0_cols = metadata_maps.condition_to_sample[
            contrast["contrast"][0]
        ]
        group_1_cols = metadata_maps.condition_to_sample[
            contrast["contrast"][1]
        ]

        contrast_list.append(
            Contrast(
                contrast=contrast["contrast"],
                genes_sheet=contrast["genes_sheet"],
                group_0_cols=group_0_cols,
                group_1_cols=group_1_cols,
                gene_list_col=contrast["gene_list_col"],
            )
        )

    return contrast_list
