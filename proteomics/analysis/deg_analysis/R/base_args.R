library(optparse)

get_thresholds <- function() {
  args_list <- list(
    make_option(
      c("--fc_threshold"),
      type = "numeric",
      help = "The log fold change threshold",
      default = 1.5
    ),
    make_option(
      c("--pval_threshold"),
      type = "numeric",
      help = "The p-value threshold",
      default = 0.05
    )
  )

  return(args_list)
}

# Returns a list of base Args for this:
#
# class DegPlotArgs(BaseModel):
#     """
#     The base arguments for all plots created after DEG analysis with limma.
#     """
#
#     # params for selecting columns in the datasframe
#     gene_column: str = Field(
#         ..., description="The column containing the gene symbols."
#     )
#     lfc_column: str = Field(
#         "logFC", description="The column containing the log fold change values"
#     )
#     pval_column: str = Field(
#         "adj.P.Val", description="The column containing the p-values"
#     )
#
#     # params for filtering the dataframe
#     fc_threshold: float = Field(1.5, description="The fold change threshold for significance")
#     pval_threshold: float = Field(0.05, description="The p-value threshold for significance")
#
#     # params for the plot aesthetics
#     width: float = Field(..., description="The width of the plot")
#     height: float = Field(..., description="The height of the plot")
#
#
# class RDegPlotArgs(DegPlotArgs):
#     gene_column: str = Field(
#         "X",
#         description="The column containing the gene symbols. "
#         "In the case of limma, this is usually 'X' - a blank column.",
#     )
#
# class RunRMixin(RConfig, abc.ABC, Generic[RResultType]):
#     input_file: FilePath = Field(..., description="The input file")
#     output_dir: DirectoryPath = Field(..., description="The output directory")
#     experiment: str = Field(..., description="The experiment name")
get_base_args <- function() {

  args_list <- list(
    # RunRMixin
    make_option(
      c("--input_file"),
      type = "character",
      help = "The file containing the differentially expressed genes"
    ),
    make_option(
      c("--output_dir"),
      type = "character",
      help = "The directory to save the results"
    ),
    make_option(
      c("--experiment"),
      type = "character",
      help = "The name of the experiment"
    ),
    # DegPlotArgs
    make_option(
      c("--gene_column"),
      type = "character",
      help = "The column containing the gene symbols. Assumes output from limma which leaves the gene column blank",
      default = "X"
    ),
    make_option(
      c("--lfc_column"),
      type = "character",
      help = "The column containing the log fold change values",
      default = "logFC"
    ),
    make_option(
      c("--pval_column"),
      type = "character",
      help = "The column containing the p-values",
      default = "adj.P.Val"
    ),
    # make_option(
    #   c("--fc_threshold"),
    #   type = "numeric",
    #   help = "The log fold change threshold",
    #   default = 1.5
    # ),
    # make_option(
    #   c("--pval_threshold"),
    #   type = "numeric",
    #   help = "The p-value threshold",
    #   default = 0.05
    # ),
    # PlotArgs
    make_option(
      c("--width"),
      type = "numeric",
      help = "The width of the plot"
    ),
    make_option(
      c("--height"),
      type = "numeric",
      help = "The height of the plot"
    )
  )

  args_list <- c(args_list, get_thresholds())

  return(args_list)
}
