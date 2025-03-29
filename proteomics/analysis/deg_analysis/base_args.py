import abc
from typing import TypeVar

from pydantic import BaseModel, Field, FilePath, DirectoryPath

from proteomics.utils.base_params import BaseParams
from proteomics.utils.run_r import RunRMixin

__all__ = [
    "DegAnalysisArgs",
    "DegPlotArgs",
    "RDegPlotArgs",
    "RConfig",
    "RunRDegAnalysis",
]


class DegAnalysisArgs(BaseModel):
    """
    Mixin for DEG analysis arguments.
    """
    # params for filtering the dataframe
    fc_threshold: float = Field(
        1.5, description="The fold change threshold for significance"
    )
    pval_threshold: float = Field(
        0.05, description="The p-value threshold for significance"
    )


class DegPlotArgs(DegAnalysisArgs):
    """
    The base arguments for all plots created after DEG analysis with limma.
    """

    # params for selecting columns in the datasframe
    gene_column: str = Field(
        ..., description="The column containing the gene symbols."
    )
    lfc_column: str = Field(
        "logFC", description="The column containing the log fold change values"
    )
    pval_column: str = Field(
        "adj.P.Val", description="The column containing the p-values"
    )

    # params for the plot aesthetics
    width: float = Field(..., description="The width of the plot")
    height: float = Field(..., description="The height of the plot")


class RDegPlotArgs(DegPlotArgs):
    """
    Mixin for all plots created after DEG analysis with limma, but for R scripts.
    """
    gene_column: str = Field(
        "X",
        description="The column containing the gene symbols. "
        "In the case of limma, this is usually 'X' - a blank column.",
    )


class RConfig(BaseParams):
    rscript_bin: FilePath = "/opt/conda/bin/Rscript"


RResultType = TypeVar("RResultType")

class RunRDegAnalysis(RunRMixin[RResultType], abc.ABC):
    """
    The base class for running R scripts for DEG analysis.
    """
    input_file: FilePath = Field(..., description="The input file")
    output_dir: DirectoryPath = Field(..., description="The output directory")
    experiment: str = Field(..., description="The experiment name")
