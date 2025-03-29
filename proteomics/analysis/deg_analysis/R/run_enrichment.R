library(clusterProfiler)
library(org.Mm.eg.db)
library(GOplot)
library(enrichplot)
library(stringr)
library(org.Hs.eg.db)
library(jsonlite)
library(optparse)
source("base_args.R")

get_gene_ids <- function(df, gene_column, org_db) {
  all_gene_id <- df[[gene_column]]
  # convert gene symbols to ENTREZID
  all_gene_id_to_ids <- bitr(
    all_gene_id,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org_db
  )

  return(
    all_gene_id_to_ids
  )
}

select_significant <- function(
  df,
  lfc_column,
  pval_column,
  fc_threshold = 1.5,
  pval_threshold = 0.05
) {
  filtered <- df[
    abs(df[[lfc_column]]) > log2(fc_threshold) &
      df[[pval_column]] < pval_threshold,
  ]

  return(filtered)
}

label_up_down <- function(
  df_filtered,
  lfc_column,
  gene_column,
  org_db
) {
  mydf <- as.data.frame(get_gene_ids(df_filtered, gene_column, org_db))

  # select the upregulated and downregulated genes
  df_down <- df_filtered[df_filtered[[lfc_column]] < 0,]

  mydf$group <- "upregulated"
  mydf$group[mydf$SYMBOL %in% df_down[[gene_column]]] <- "downregulated"

  return(mydf)

}

# Returns a json string of type:
#
# type EnrichResult = {
#     file_path: string;
#     ont: "BP" | "CC" | "MF" | "KEGG" | "None";
#     result_type: "enrichResult" | "compareResult" | "plot" | "geneIds";
# }
make_location <- function(
  file_path,
  ont,
  result_type
) {

  data <- list(
    file_path = file_path,
    ont = ont,
    result_type = result_type
  )

  json_data <- toJSON(data, auto_unbox = TRUE)

  return(json_data)
}


plot_ontology <- function(ont_result,
                          output_dir,
                          experiment,
                          ont,
                          width = 6.5,
                          height = 6.5) {

  p1 <- dotplot(ont_result)

  # add the title
  p1 <- p1 + ggtitle(paste0(ont, " for ", experiment))

  formats <- c("png", "pdf")

  for (format in formats) {
    output_file <- file.path(output_dir, paste0(ont, "_plot.", format))
    ggsave(output_file,
           p1,
           width = width,
           height = height
    )
  }

  print(
    paste0("Saved ", ont, " plot")
  )

  return(
    make_location(file.path(output_dir, paste0(ont, "_plot.", "png")), ont, "plot")
  )
}


ego <- function(ont,
                gene_ids,
                go_output_dir,
                experiment,
                org_db,
                organism,
                width = 6.5,
                height = 6.5
) {
  print(paste0("Getting ", ont))

  df_with_up_down <- NULL
  formula_res <- NULL

  if (ont == "KEGG") {
    formula_res <- enrichKEGG(
      gene = gene_ids$ENTREZID,
      organism = organism,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )
    df_with_up_down <- compareCluster(ENTREZID ~ group,
                                      data = gene_ids,
                                      fun = "enrichKEGG",
                                      organism = organism
    )

  } else if (ont == "BP" || ont == "CC" || ont == "MF") {

    formula_res <- enrichGO(
      gene = gene_ids$ENTREZID,
      OrgDb = org_db,
      ont = ont,
      readable = TRUE
    )
    df_with_up_down <- compareCluster(ENTREZID ~ group,
                                      data = gene_ids,
                                      fun = "enrichGO",
                                      pvalueCutoff = 0.05,
                                      pAdjustMethod = "BH",
                                      OrgDb = org_db,
                                      ont = ont
    )
  } else {
    print("Invalid ontology")
    return()
  }

  # try catch in case there are no terms
  if (is.null(formula_res)) {
    print(paste0("No ", ont))
    return()
  }

  print("Writing results")
  enrich_csv <- file.path(go_output_dir, paste0(ont, "_enrich.csv"))
  contrast_csv <- file.path(go_output_dir, paste0(ont, "_contrast.csv"))

  write.csv(formula_res, enrich_csv, row.names = FALSE)
  write.csv(df_with_up_down, contrast_csv, row.names = FALSE)

  print(paste0("Plotting ", ont))

  locations <- list()

  locations <- tryCatch(
  {
    c(locations, plot_ontology(
      ont_result = formula_res,
      ont = ont,
      output_dir = go_output_dir,
      experiment = experiment,
      width = width,
      height = height
    )
    )
  },
    error = function(e) {
      print(paste("Error in plot_ontology function:", e))
      return(NULL)
    }
  )

  locations <- c(locations, make_location(enrich_csv, ont, "enrichResult"))
  locations <- c(locations, make_location(contrast_csv, ont, "compareResult"))

  return(locations)

}

# Run enrichment analysis
#
# deg_file: path to the file containing the differentially expressed genes
# gene_column: the column containing the gene symbols
# lfc_column: the column containing the log fold change values
# pval_column: the column containing the p-values
# fc_threshold: the log fold change threshold
# pval_threshold: the p-value threshold
# output_dir: the directory to save the results
# organism: the organism to use for the enrichment analysis. Must be one of "mmu" or "hsa"
run_enrichment <- function(
  deg_file,
  gene_column,
  lfc_column,
  pval_column,
  fc_threshold,
  pval_threshold,
  output_dir,
  organism,
  experiment,
  width,
  height
) {

  org_db <- NULL
  if (organism == "mmu") {
    org_db <- org.Mm.eg.db
  } else if (organism == "hsa") {
    org_db <- org.Hs.eg.db
  } else {
    print("Invalid organism")
    return()
  }

  df <- read.csv(deg_file, header = TRUE)
  # select the significant genes

  # get the gene ids
  print("Getting gene ids")
  gene_ids <- label_up_down(
    select_significant(df,
                       lfc_column = lfc_column,
                       pval_column = pval_column,
                       fc_threshold = fc_threshold,
                       pval_threshold = pval_threshold
    ),
    lfc_column = lfc_column,
    gene_column = gene_column,
    org_db = org_db
  )

  # write the gene ids to a file
  gene_ids_file <- file.path(output_dir, "gene_ids.csv")
  write.csv(gene_ids, gene_ids_file, row.names = FALSE)

  ontologies <- c("BP", "CC", "MF", "KEGG")

  all_results <- list()

  print(
    list(
      make_location(gene_ids_file, "None", "geneIds")
    )
  )

  for (ont in ontologies) {
    # try catch this
    results <- tryCatch(
    {
      ego(
        ont = ont,
        gene_ids = gene_ids,
        go_output_dir = output_dir,
        experiment = experiment,
        org_db = org_db,
        organism = organism,
        width = width,
        height = height
      )
    },
      error = function(e) {
        print(paste("Error in ego function:", e))
        return(NULL)
      }

    )


    if (!is.null(results)) {
      all_results <- c(all_results, results)
    }
  }

  print("Finished enrichment analysis")

  # Flatten the list while keeping the list structure
  flattened_list <- do.call(c, all_results)

  for (val in flattened_list) {
    cat(
      val[1],
      "\n"
    )
  }
}


get_args <- function() {
  option_list <- list(
    make_option(
      c("--organism"),
      type = "character",
      help = "The organism to use for the enrichment analysis. Must be one of 'mmu' or 'hsa'"
    )
  )

  base_args_list <- get_base_args()
  # merge the two lists
  option_list <- c(base_args_list, option_list)

  # Create OptionParser object
  opt <- parse_args(OptionParser(option_list = option_list))

  return(opt)
}


main <- function() {

  opt <- get_args()

  # check if the output directory exists
  if (!dir.exists(opt$output_dir)) {
    stop("Output directory does not exist")
  }

  run_enrichment(
    deg_file = opt$input_file,
    gene_column = opt$gene_column,
    lfc_column = opt$lfc_column,
    pval_column = opt$pval_column,
    fc_threshold = opt$fc_threshold,
    pval_threshold = opt$pval_threshold,
    output_dir = opt$output_dir,
    organism = opt$organism,
    experiment = opt$experiment,
    width = opt$width,
    height = opt$height
  )

}


main()
