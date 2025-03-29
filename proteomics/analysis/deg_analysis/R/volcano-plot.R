library(ggplot2)
library(EnhancedVolcano)
# Load the necessary library
library(gridExtra)
library(optparse)
source("base_args.R")

get_args <- function() {
  option_list <- list(
    make_option(
      c("--upregulated_color"),
      type = "character",
      default = "red",
      help = "The color for upregulated genes. Default is red"
    ),
    make_option(
      c("--downregulated_color"),
      type = "character",
      default = "royalblue",
      help = "The color for downregulated genes. Default is royalblue"
    ),
    make_option(
      c("--not_sig_color"),
      type = "character",
      default = "black",
      help = "The color for non-significant genes. Default is black"
    )
  )

  base_args_list <- get_base_args()
  # merge the two lists
  option_list <- c(base_args_list, option_list)

  # Create OptionParser object
  opt <- parse_args(OptionParser(option_list = option_list))

  return(opt)
}


run_volcano <- function(
  input_file,
  output_dir,
  genes_column = "X",
  # the column name for the gene names from Limma. By default its blank (X)
  log2_fc_col = "logFC",
  p_val_col = "adj.P.Val",
  p_cutoff = 0.05,
  fc_cutoff = 1.5,
  upregulated_color = "red",
  downregulated_color = "royalblue",
  not_sig_color = "black",
  plot_title = "Volcano Plot"
) {

  lfc_cutoff <- log2(fc_cutoff)

  # Load in the data (limma results)
  lfc_data <- read.csv(input_file)

  # Create a dataframe from the data and set the row names to the Gene column
  data_df <- as.data.frame(lfc_data, row.names = lfc_data[[genes_column]])

  num_de <- nrow(
    # must be both greater than the lfc_cutoff and less than the p_cutoff
    data_df[which(abs(data_df[[log2_fc_col]]) > lfc_cutoff & data_df[[p_val_col]] < p_cutoff),]
  )

  # create a colormap
  keyvals <- ifelse(
    data_df[[log2_fc_col]] < -lfc_cutoff & data_df[[p_val_col]] < p_cutoff, downregulated_color,
    ifelse(data_df[[log2_fc_col]] > lfc_cutoff & data_df[[p_val_col]] < p_cutoff, upregulated_color,
           not_sig_color
    )
  )

  caption <- paste0("Total: ", nrow(data_df), ". Significant: ", num_de)

  print(caption)

  keyvals[is.na(keyvals)] <- not_sig_color
  names(keyvals)[keyvals == upregulated_color] <- "Up"
  names(keyvals)[keyvals == not_sig_color] <- "Not significant"
  names(keyvals)[keyvals == downregulated_color] <- "Down"

  print("Creating plot")

  plot <- EnhancedVolcano(
    data_df,
    lab = rownames(data_df),
    selectLab = c(""),
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    x = log2_fc_col,
    y = p_val_col,
    title = plot_title,
    FCcutoff = lfc_cutoff,
    pCutoff = p_cutoff,
    caption = caption,
    subtitle = "",
    borderWidth = 1.0,
    borderColour = "black",
    border = "full",
    colCustom = keyvals,
    subtitleLabSize = -1,
    legendPosition = "right",
  )

  # Set the limits of the plot
  maxy <- ceiling(-log10(min(data_df[[p_val_col]], na.rm = TRUE)))
  minx <- floor(min(data_df[[log2_fc_col]], na.rm = TRUE))
  maxx <- ceiling(max(data_df[[log2_fc_col]], na.rm = TRUE))
  plot <- plot + coord_cartesian(ylim = c(NA, maxy), xlim = c(minx, maxx))

  print("Saving plot")
  for (format in c("svg", "pdf", "png")) {
    plot_name <- file.path(output_dir, paste0("volcano_plot.", format))
    ggsave(plot_name, plot, dpi = 400, width = 9, height = 8.5)
  }

  print(
    file.path(output_dir, paste0("volcano_plot.", format))
  )
}

main <- function() {
  args <- get_args()

  if (!dir.exists(args$output_dir)) {
    stop("Output directory does not exist")
  }

  run_volcano(
    input_file = args$input_file,
    output_dir = args$output_dir,
    genes_column = args$gene_column,
    log2_fc_col = args$lfc_column,
    p_val_col = args$pval_column,
    fc_cutoff = args$fc_threshold,
    p_cutoff = args$pval_threshold,
    upregulated_color = args$upregulated_color,
    downregulated_color = args$downregulated_color,
    not_sig_color = args$not_sig_color,
    plot_title = args$experiment
  )
}

main()
