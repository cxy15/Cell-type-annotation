#!/usr/bin/env Rscript
# Step 2: GO enrichment (human: org.Hs.eg.db, mouse: org.Mm.eg.db)
options(warn = -1)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(tidyverse)
  library(jsonlite)
})

species <- tolower(Sys.getenv("SPECIES", "human"))
out_root <- Sys.getenv("CELLTYPE_PIPELINE_OUT", ".")
exchange <- file.path(out_root, "Exchange")
markers_path <- file.path(exchange, "all_markers.csv")
out_csv <- file.path(exchange, "enrichment_summary.csv")

if (!file.exists(markers_path)) {
  stop("Missing markers file: ", markers_path, " — run step 1 first.")
}

if (species %in% c("mouse", "mm", "mus")) {
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  org_db <- org.Mm.eg.db
  message("[02] Species: mouse (org.Mm.eg.db)")
} else {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  org_db <- org.Hs.eg.db
  message("[02] Species: human (org.Hs.eg.db)")
}

all_markers <- read.csv(markers_path, stringsAsFactors = FALSE)
clusters <- unique(all_markers$celltype)

result_table <- data.frame(Cluster = clusters, stringsAsFactors = FALSE)

for (i in seq_along(clusters)) {
  cluster_id <- clusters[i]

  cluster_genes <- all_markers %>%
    dplyr::filter(celltype == cluster_id) %>%
    dplyr::pull(gene)

  gene_ids <- tryCatch(
    bitr(
      cluster_genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org_db
    ),
    error = function(e) NULL
  )

  if (is.null(gene_ids) || nrow(gene_ids) == 0) {
    next
  }

  go_result <- tryCatch(
    enrichGO(
      gene = gene_ids$ENTREZID,
      OrgDb = org_db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.01,
      qvalueCutoff = 0.05,
      readable = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(go_result) || nrow(as.data.frame(go_result)) == 0) {
    next
  }

  go_df <- as.data.frame(go_result)
  go_df <- go_df %>%
    dplyr::mutate(GeneRatioNum = sapply(GeneRatio, function(x) eval(parse(text = x)))) %>%
    dplyr::arrange(dplyr::desc(GeneRatioNum)) %>%
    dplyr::slice_head(n = 5)

  for (j in 1:5) {
    col_name <- paste0("Top", j, "_Pathway")
    if (j <= nrow(go_df)) {
      result_table[i, col_name] <- go_df$Description[j]
    } else {
      result_table[i, col_name] <- NA_character_
    }
  }
}

result_table <- result_table %>% dplyr::arrange(as.numeric(as.character(Cluster)))
write.csv(result_table, out_csv, row.names = FALSE)
message("[02] Wrote ", out_csv)

prog <- file.path(exchange, "progress.json")
prev <- if (file.exists(prog)) jsonlite::fromJSON(prog, simplifyVector = TRUE) else list()
prev$step <- 2L
prev$status <- "done"
prev$enrichment_csv <- out_csv
jsonlite::write_json(prev, prog, auto_unbox = TRUE, pretty = TRUE)
message("[02] Step 2 complete.")
