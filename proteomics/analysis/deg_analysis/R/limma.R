library(limma)
library(optparse)
source("base_args.R")

load_count_matrix <- function(count_filename) {
  # Load the count matrix
  count_matrix <- read.csv(count_filename, row.names = 1)

  # cat(
  #   "Loaded count matrix with",
  #   nrow(count_matrix),
  #   "genes and",
  #   ncol(count_matrix),
  #   "samples\n"
  # )

  return(count_matrix)
}

load_metadata <- function(metadata_filename) {
  # Load the metadata
  metadata <- read.csv(metadata_filename)
  return(metadata)
}

check_count_metadata <- function(count_df,
                                 meta_df,
                                 meta_sample_col = "Sample",
                                 meta_group_col = "Group") {
  # ensure all the columns in the metadata are in the counts
  col_cols <- colnames(count_df)
  meta_samples <- meta_df[[meta_sample_col]]

  if (!all(meta_samples %in% col_cols)) {
    stop("Not all metadata samples are in the count matrix")
  }

  # ensure that there are only 2 groups
  if (length(unique(meta_df[[meta_group_col]])) != 2) {
    stop("Number of groups in metadata is not equal to 2")
  }

}

create_design_matrix <- function(meta_df,
                                 contrast,
                                 group_col = "Group") {
  # From the metadata, create the design matrix and contrasts

  groups <- factor(meta_df[[group_col]])

  design <- model.matrix(~0 + groups)
  colnames(design) <- levels(groups)

  contrasts <- makeContrasts(
    contrasts = contrast,
    levels = design
  )

  return(list(design = design, contrasts = contrasts))
}


fit_limma <- function(
  count_filename,
  metadata_filename,
  contrast
) {
  count_df <- load_count_matrix(count_filename)
  meta_df <- load_metadata(metadata_filename)

  check_count_metadata(count_df, meta_df)

  # create the design matrix
  design_contrast <- create_design_matrix(meta_df, contrast)

  design <- design_contrast$design
  contrasts <- design_contrast$contrasts

  # Fit the expression matrix to a linear model
  fit <- lmFit(count_df, design)
  # Compute contrast
  fit_contrast <- contrasts.fit(fit, contrasts)
  # Bayes statistics of differential expression
  # *There are several options to tweak!*
  fit_contrast <- eBayes(fit_contrast)

  return(topTable(fit_contrast, number = Inf))

}

get_args <- function() {
  option_list <- list(
    make_option(c("--counts"), type = "character", help = "Path to the counts file"),
    make_option(c("--metadata"), type = "character", help = "Path to the metadata file"),
    make_option(c("--output_dir"), type = "character", help = "Output directory for results"),
    make_option(c("--sig_output_dir"), type = "character", help = "Output directory for significant results"),
    make_option(c("--contrast_name"), type = "character", help = "Name of the contrast"),
    make_option(c("--contrast_1"), type = "character", help = "First group for contrast"),
    make_option(c("--contrast_2"), type = "character", help = "Second group for contrast")
    # make_option(c("--p_sig"), type = "numeric", default = 0.05, help = "Significance p-value threshold"),
    # make_option(c("--fc"), type = "numeric", default = 1.5, help = "Fold change threshold")
  )

  option_list <- c(option_list, get_thresholds())

  # Create OptionParser object
  opt <- parse_args(OptionParser(option_list = option_list))

  return(opt)
}


main <- function() {

  args <- get_args()

  # check if output dir exists - if not throw an error
  if (!dir.exists(args$output_dir) || !dir.exists(args$sig_output_dir)) {
    stop("Output directories does not exist")
  }

  contrast_formula <- paste0(args$contrast_name, " = ", args$contrast_1, " - ", args$contrast_2)

  print("Running limma with the following contrast formula:")
  print(contrast_formula)

  result <- fit_limma(
    args$counts,
    args$metadata,
    contrast_formula
  )

  deg_file <- paste0(args$output_dir, "/", args$contrast_name, "_deg_limma.csv")
  write.csv(result, deg_file)

  # filter for padj < 0.05 and FC > 1.5
  fc_sig <- log2(args$fc_threshold)

  sig_genes <- result[result$adj.P.Val < args$pval_threshold & abs(result$logFC) >= fc_sig,]

  num_sig <- paste0(
    "There are ",
    nrow(sig_genes),
    " significant genes with p-value < ",
    args$pval_threshold,
    " and fold change > ",
    args$fc_threshold,
    " for contrast ",
    args$contrast_name,
    " with formula '",
    contrast_formula,
    "'"
  )

  print(num_sig)

  print("Writing significant genes to file")
  writeLines(
    num_sig,
    con = paste0(args$sig_output_dir, "/", args$contrast_name, "_deg_limma_sig.txt")
  )
  write.csv(sig_genes, paste0(args$sig_output_dir, "/", args$contrast_name, "_deg_limma_sig.csv"))

  cat(
    deg_file,
    "\n"
  )
}

main()
