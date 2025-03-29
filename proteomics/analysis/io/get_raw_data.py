from pathlib import Path
from typing import Annotated

import pandas as pd

__all__ = [
    "load_imputed_counts",
    "ImputedIntensity",
]

from pydantic import AfterValidator, FilePath

from proteomics.analysis.io.load_metadata import MetadataMaps
from proteomics.utils.base_params import BaseParams, is_xlsx_file


class ImputedIntensity(BaseParams):
    input_file: Annotated[FilePath, AfterValidator(is_xlsx_file)]
    index_col: int = 0
    sheet_name: str

    # metadata
    metadata_file: Annotated[FilePath, AfterValidator(is_xlsx_file)]
    metadata_sheet_name: str
    metadata_index_col: int = 0
    condition_col: str = "condition"


def load_imputed_counts(
    input_file: Path,
    sheet_name: str,
    index_col: int,
) -> pd.DataFrame:
    """
    Loads a full rank counts matrix from an excel file.
    """
    df = pd.read_excel(input_file, index_col=index_col, sheet_name=sheet_name)

    # check to make sure there are NO missing values
    if df.isnull().values.any():
        raise ValueError("There are missing values in the data.")

    return df


def load_orig_counts(
        input_file: Path,
        sheet_name: str,
        genes_col: str,
        metadata_maps: MetadataMaps,
) -> pd.DataFrame:
    """
    Loads the original counts matrix from an excel file.
    """
    df = pd.read_excel(input_file, sheet_name=sheet_name)

    # drop duplicates from genes_col
    print(f"Loaded {df.shape[0]} genes and {df.shape[1]} columns")
    df = df.drop_duplicates(subset=[genes_col], keep=False)
    print(f"Dropped duplicates. Remaining: {df.shape[0]}")

    df = df.set_index(genes_col)

    df = df[
        [
            x for x in df.columns if x in metadata_maps.sample_to_condition
        ]
    ]
    print(f"There are {df.shape[1]} observations with metadata")

    return df