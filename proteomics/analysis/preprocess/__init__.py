from typing import Annotated, Literal, TypedDict

from pydantic import FilePath, AfterValidator

from .normalize import normalize
from proteomics.utils.base_params import BaseParams, is_xlsx_file
from ..io.load_metadata import ContrastInput


class PreprocessParams(BaseParams):
    norm_method: Literal["center.median", "standardize"] = "center.median"

    # plots
    n_rows: int
    n_cols: int

    # contrasts
    gene_input_file: (
        Annotated[FilePath, AfterValidator(is_xlsx_file)] | None
    ) = None
    contrasts: list[ContrastInput]


__all__ = ["normalize", "PreprocessParams", "PreprocessParams"]
