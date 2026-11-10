# R 4.4+ with Matrix. Read frozen RDS slots without loading the Seurat package.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || !args[1] %in% c("skin", "follicle")) {
  stop("Usage: Rscript tools/export.R skin|follicle RDS ARCHIVE [demo|all]")
}
kind <- args[1]
archive <- args[3]
mode <- if (length(args) >= 4L) args[4] else "all"
stopifnot(mode %in% c("all", "demo"))
output <- file.path(if (mode == "demo") "preview" else "docs", "data", kind)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
read_table <- function(relative) read.csv(file.path(archive, relative), stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

obj <- readRDS(args[2])
meta <- attr(obj, "meta.data")
assay <- attr(obj, "assays")[["RNA"]]
counts <- attr(assay, "counts")
stopifnot(inherits(counts, "dgCMatrix"), !anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)))
embedding <- attr(attr(obj, "reductions")[[if (kind == "skin") "umap" else "hf_umap"]], "cell.embeddings")
rm(obj, assay)
invisible(gc())
stopifnot(setequal(colnames(counts), rownames(embedding)), setequal(colnames(counts), rownames(meta)))
embedding <- embedding[colnames(counts), 1:2, drop = FALSE]
meta <- meta[colnames(counts), , drop = FALSE]
if (kind == "skin") {
  check <- read_table("figures/Fig4/source_data/03_selected_celltype_marker_DotPlot_cluster_annotation.csv")
  labels <- check$cell_type[match(as.character(meta$seurat_clusters), check$cluster)]
  palette <- c("Suprabasal IFE"="#A6CEE3","Immune"="#1F78B4","Basal IFE"="#B2DF8A","ORS"="#33A02C",
    "Terminal IFE"="#FDBF6F","ORS/HFSC_like"="#E31A1C","uHF"="#FF7F00","CP"="#9C6BB5",
    "Cycling basal"="#FB9A99","FIB"="#B15928","Mixed"="#CAB2D6","IRS/CX"="#FFFF99")
} else {
  labels <- as.character(meta$hf_subtype_label)
  colors <- read_table("inputs/manual_annotations/T2_HF_celltype_colors.csv")
  palette <- setNames(colors$color, colors$cell_type)
  check <- read_table("figures/Fig5/source_data/04_HF_secondary_annotation_main_marker_DotPlot_cluster_annotation.csv")
}
types <- names(palette)
if (!setequal(types, labels)) stop("Atlas labels and palette differ: ", paste(setdiff(labels, types), collapse = ", "), " | ", paste(setdiff(types, labels), collapse = ", "))
stopifnot(all(meta$condition %in% c("X", "M")),
          nrow(embedding) < 65536L, all(is.finite(embedding)))
totals <- rowsum(check$n_cells, check$cell_type)
stopifnot(all(as.integer(table(factor(labels, levels = types))) == totals[types, 1]))
write.table(data.frame(cell_type = types, color = unname(palette)),
  file.path(output, "cell-types.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

coords <- file(file.path(output, "cells.bin"), "wb")
writeBin(as.double(embedding[, 1]), coords, size = 4, endian = "little")
writeBin(as.double(embedding[, 2]), coords, size = 4, endian = "little")
writeBin(as.raw(match(labels, types) - 1L), coords)
writeBin(as.raw(meta$condition == "M"), coords)
close(coords)

totals <- Matrix::colSums(counts)
stopifnot(all(totals > 0))
genes <- rownames(counts)
if (mode == "demo") genes <- intersect(c("Krt79", "Il11ra1", "Nr4a1"), genes)
stopifnot(length(genes) > 0)
by_gene <- methods::as(Matrix::t(counts[genes, , drop = FALSE]), "CsparseMatrix")
rm(counts)
invisible(gc())
scale <- 4096
chunk <- 24L
writeLines(genes, file.path(output, "genes.txt"), useBytes = TRUE)
ptr <- by_gene@p
for (start in seq(1L, length(genes), by = chunk)) {
  last <- min(start + chunk - 1L, length(genes))
  block <- gzfile(file.path(output, sprintf("%04d.bin", (start - 1L) %/% chunk)), "wb", compression = 5)
  writeBin(as.integer(diff(ptr[start:(last + 1L)])), block, size = 4, endian = "little")
  for (i in start:last) {
    if (ptr[i] == ptr[i + 1L]) next
    loc <- (ptr[i] + 1L):ptr[i + 1L]
    cell <- by_gene@i[loc]
    expression <- log1p(by_gene@x[loc] / totals[cell + 1L] * 10000)
    code <- as.integer(round(expression * scale))
    stopifnot(all(code >= 0L & code <= 65535L))
    writeBin(cell, block, size = 2, endian = "little")
    writeBin(code, block, size = 2, endian = "little")
  }
  close(block)
  if ((start - 1L) %% 2400L == 0L) message(kind, " genes ", last, "/", length(genes))
}
writeLines(paste(nrow(embedding), chunk, scale, sep = "\t"), file.path(output, "format.tsv"))
message(kind, " ", mode, " complete: ", length(genes), " genes; ", nrow(embedding), " cells")
