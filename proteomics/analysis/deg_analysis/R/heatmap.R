create_heatmap_df <- function(
  counts_file,
  limma_results_file,
  remove_from_colname = "",
  fc_threshold = 1.5,
  pval_threshold = 0.05,
  lfc_column = "logFC",
  pval_column = "adj.P.Val"
) {
  df_counts <- read.csv(
    counts_file,
    row.names = 1
  )

  df_de <- read.csv(
    limma_results_file,
    row.names = 1
  )

  fc_threshold <- log2(fc_threshold)

  df_de <- df_de[
    abs(df_de[[lfc_column]]) >= fc_threshold &
      df_de[[pval_column]] <= pval_threshold,
  ]

  print(
    paste0(
      "Creating a heatmap for ",
      dim(df_de)[1],
      " significant genes"
    )
  )

  df_counts <- df_counts[rownames(df_de),]

  if (remove_from_colname != "") {
    print(
      paste0(
        "Replacing column names with '",
        remove_from_colname,
        "' using regex"
      )
    )

    colnames(df_counts) <- gsub(remove_from_colname, "", colnames(df_counts))
  }

  return(df_counts)

}

make_z_score <- function(df, direction = 'row') {
  # if direction is row, then the genes are the rows
  if (direction == 'row') {
    df <- t(scale(t(df)))
  } else {
    df <- scale(df)
  }

  return(df)
}

get_sample_annotation <- function(
  df,
  remove_from_name = "",
  row_or_col = "col"
  # are the samples in the columns or rows?
) {
  if (
    row_or_col == "col"
  ) {
    annotation <- data.frame(
      row.names = colnames(df),
      condition = gsub(remove_from_name, "", colnames(df))
    )
  } else {
    annotation <- data.frame(
      row.names = rownames(df),
      condition = gsub(remove_from_name, "", rownames(df))
    )
  }

  annotation$condition <- as.factor(annotation$condition)
  return(annotation)

}

get_sample_colors <- function(
  annotation_df
) {
  cats <- levels(annotation_df$condition)

  # Create a color palette with as many colors as there are unique conditions
  color_palette <- brewer.pal(length(cats), "Set1")

  # if there are more colors than conditions, shorten the color palette
  if (length(color_palette) > length(cats)) {
    color_palette <- color_palette[1:length(cats)]
  }

  # Create a named vector mapping conditions to colors
  condition_colors <- setNames(color_palette, cats)

  return(list(condition = condition_colors))
}


plot_and_save <- function(
  df_heat,  # a STANDARDIZED dataframe - Must be standatdized beforehand
  annotation_df,  # a dataframe with sample annotations
  my_color,  # a list of colors for the sample annotations
  output_dir,
  plot_title,
  # pheatmap options
  cluster_cols,
  cluster_rows,
  show_colnames,
  show_rownames,
  fontsize = 8
) {
  my_heatmap <- pheatmap(df_heat,
                         show_colnames = TRUE,
                         cluster_cols = cluster_cols,
                         cluster_rows = cluster_rows,
                         show_colnames = show_colnames,
                         show_rownames = show_rownames,
                         show_rownames = FALSE,
                         border_color = FALSE,
                         col = colorRampPalette(c("royalblue", "white", "red"))(100),
                         main = plot_title,
                         fontsize = fontsize,
                         annotation_col = annotation_df,
                         annotation_colors = my_color
  )

  formats <- c("pdf", "png", "svg")

  for (format in formats) {
    output_file <- file.path(heat_out_dir, paste0(experiment, ifelse(cluster_col, "_cluster", ""), "_heatmap.", format))

    ggsave(output_file, my_heatmap, width = 4.5, height = 6)

  }
}
