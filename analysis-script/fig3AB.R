if (dir.exists(".rlib")) .libPaths(c(".rlib", .libPaths()))   # optional project-local R library
suppressMessages({library(Seurat); library(Matrix); library(SingleR); library(scClassify); library(BiocParallel); library(parallel)
  library(ggplot2); library(ggrepel); library(dplyr)})
raw_dir <- "GSE178667_RAW"; out_dir <- "results"
cache <- file.path(out_dir, "cache"); dir.create(cache, showWarnings = FALSE, recursive = TRUE)
python <- Sys.getenv("CELLTYPIST_PYTHON", if (file.exists(".venv_celltypist/bin/python")) ".venv_celltypist/bin/python" else "python3"); n_cores <- as.integer(Sys.getenv("N_CORES", "4"))
set.seed(1)

## ---- 1. raw counts ----
read_sample <- function(gsm, tag) {
  f <- function(x) file.path(raw_dir, sprintf("%s_%s_%s", gsm, tag, x))
  m <- ReadMtx(mtx = f("matrix.mtx.gz"), cells = f("barcodes.tsv.gz"), features = f("features.tsv.gz"))
  colnames(m) <- paste(gsm, colnames(m), sep = "_"); m }
counts <- cbind(read_sample("GSM5395694", "lvpl1"), read_sample("GSM5395695", "lvpl3"))
cat("raw:", ncol(counts), "cells,", nrow(counts), "genes\n")

## ---- 2. reproduce the published clustering ----
# The paper does not state its QC cutoffs; these reproduce its 10,024 cells and its per-cluster statistics exactly.
so <- CreateSeuratObject(counts); so$sample <- sub("_.*", "", colnames(so))
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25 & nCount_RNA < 25000)
so <- NormalizeData(so, verbose = FALSE); so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE); so <- RunPCA(so, npcs = 50, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:25, verbose = FALSE); so <- FindClusters(so, resolution = 0.3, verbose = FALSE)
so <- RunUMAP(so, dims = 1:25, verbose = FALSE)
cat("after QC:", ncol(so), "cells,", nlevels(so$seurat_clusters), "clusters\n")
stopifnot(ncol(so) == 10024, nlevels(so$seurat_clusters) == 22)

## ---- 3. published cell type names ----
# Cell type of each published cluster, from Supplementary Table 1 of Youkilis & Bassnett (2021)
s1 <- data.frame(cluster = 1:22, celltype = c("Rods", "U - Fibroblast", "Amacrine", "Macrophage", "Cone Bipolar", "Melanocyte", "M\u00fcller Glia", "RPE",
  "PCE", "S - Fibroblast", "Rod Bipolar", "NPCE", "Cones", "Vascular Endothelium", "rM Macrophage", "N/K T-Cell", "Erythrocyte", "Smooth Muscle",
  "Neutrophil", "Schwann Cell", "B - Cell", "RGC"))
types <- s1$celltype
lab <- setNames(s1$celltype[match(as.integer(as.character(so$seurat_clusters)) + 1L, s1$cluster)], colnames(so))
stopifnot(!anyNA(lab))
# check against cluster sizes quoted in the paper (C2 1,622; C6 537; C9 355; C10 346; C12 304 cells)
stopifnot(all(table(factor(lab, levels = types))[types[c(2, 6, 9, 10, 12)]] == c(1622, 537, 355, 346, 304)))
smp <- setNames(so$sample, colnames(so))

## ---- 4. annotation tools on intact and perturbed data ----
cnt <- GetAssayData(so, layer = "counts"); cnt <- cnt[rowSums(cnt > 0) >= 10, ]                       # genes detected in at least 10 cells
lognorm <- function(m) { m <- as(m, "CsparseMatrix"); m@x <- log1p(m@x / rep(colSums(m), diff(m@p)) * 1e4); m }   # Seurat's default normalisation
tools <- list(                                                                                           # each takes counts C, training and test cell names
  "Seurat" = function(C, tr, te) {
    ref <- CreateSeuratObject(C[, tr]); ref <- NormalizeData(ref, verbose = FALSE); ref <- FindVariableFeatures(ref, nfeatures = 2000, verbose = FALSE)
    ref <- ScaleData(ref, verbose = FALSE); ref <- RunPCA(ref, npcs = 25, verbose = FALSE)
    qry <- NormalizeData(CreateSeuratObject(C[, te]), verbose = FALSE)
    an <- FindTransferAnchors(reference = ref, query = qry, dims = 1:25, reference.reduction = "pca", verbose = FALSE)
    data.frame(cell = te, label = TransferData(anchorset = an, refdata = unname(lab[tr]), dims = 1:25, verbose = FALSE)$predicted.id) },
  "SingleR" = function(C, tr, te) {
    L <- lognorm(C); p <- SingleR(test = L[, te], ref = L[, tr], labels = unname(lab[tr]), de.method = "wilcox", BPPARAM = SerialParam())
    data.frame(cell = te, label = ifelse(is.na(p$pruned.labels), "Unassigned", p$pruned.labels)) },
  "CellTypist" = function(C, tr, te) {                                                                  # Python; fixed seed so retraining is deterministic
    d <- tempfile("celltypist_"); dir.create(d); cells <- c(tr, te)
    writeMM(C[, cells], file.path(d, "counts.mtx")); writeLines(rownames(C), file.path(d, "genes.txt"))
    write.csv(data.frame(cell = cells, set = rep(c("train", "test"), c(length(tr), length(te))), label = unname(lab[cells])), file.path(d, "cells.csv"), row.names = FALSE)
    writeLines(c("import sys, scipy.io, pandas as pd, anndata as ad, scanpy as sc, celltypist",
                 "d = sys.argv[1]",
                 "a = ad.AnnData(scipy.io.mmread(d + '/counts.mtx').T.tocsr(), obs=pd.read_csv(d + '/cells.csv', index_col=0))",
                 "a.var_names = [g.strip() for g in open(d + '/genes.txt')]",
                 "sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)",
                 "tr = a[a.obs['set'] == 'train'].copy(); te = a[a.obs['set'] == 'test'].copy()",
                 "m = celltypist.train(tr, labels='label', n_jobs=2, feature_selection=True, check_expression=False, random_state=0)",
                 "celltypist.annotate(te, model=m, majority_voting=False).predicted_labels[['predicted_labels']].to_csv(d + '/pred.csv')"), file.path(d, "run.py"))
    if (system2(python, c(file.path(d, "run.py"), d), stdout = file.path(d, "log.txt"), stderr = file.path(d, "log.txt")) != 0) stop(paste(tail(readLines(file.path(d, "log.txt")), 8), collapse = "\n"))
    pr <- read.csv(file.path(d, "pred.csv"), row.names = 1, stringsAsFactors = FALSE); unlink(d, recursive = TRUE)
    data.frame(cell = te, label = pr[te, "predicted_labels"]) },
  "scClassify" = function(C, tr, te) {
    L <- lognorm(C); r <- scClassify(exprsMat_train = L[, tr], cellTypes_train = unname(lab[tr]), exprsMat_test = list(q = L[, te]), tree = "HOPACH", algorithm = "WKNN",
                                     selectFeatures = "limma", similarity = "pearson", returnList = FALSE, verbose = FALSE)
    data.frame(cell = te, label = as.character(r$testRes$q$pearson_WKNN_limma$predRes[te])) })

# reference/query split by experiment, both directions
e1 <- names(smp)[smp == "GSM5395694"]; e2 <- names(smp)[smp == "GSM5395695"]
split <- list(list(train = e1, test = e2), list(train = e2, test = e1))
perts <- c(thin90 = 0.9, thin80 = 0.8, genes90 = 0.9, genes80 = 0.8)          # share of UMIs (thin) or genes (genes) kept
pert_code <- c(genes80 = 3, thin90 = 4, thin80 = 5, genes90 = 7); n_rep <- 5  # fixed seed code per perturbation
runs <- list()
for (i in 1:2) runs[[sprintf("split_%d", i)]] <- c(split[[i]], list(perturb = "none", seed = i))
for (p in names(perts)) for (r in seq_len(n_rep)) for (i in 1:2)
  runs[[sprintf("perturb_%s_r%d_%d", p, r, i)]] <- c(split[[i]], list(perturb = p, seed = 1000 * pert_code[[p]] + 10 * r + i))
make_counts <- function(run) { set.seed(run$seed); tr <- run$train; te <- run$test
  if (run$perturb == "none") return(cnt[, c(tr, te)])
  keep <- perts[[run$perturb]]
  if (startsWith(run$perturb, "thin")) { Q <- cnt[, te]; Q@x <- as.numeric(rbinom(length(Q@x), Q@x, keep)); return(cbind(cnt[, tr], drop0(Q))) }   # each UMI of a query cell kept with probability `keep`
  cnt[sort(sample(nrow(cnt), round(keep * nrow(cnt)))), c(tr, te)] }                                                                             # same genes removed from reference and query

jobs <- expand.grid(run = names(runs), tool = names(tools), stringsAsFactors = FALSE); jobs$file <- file.path(cache, sprintf("%s__%s.rds", jobs$tool, jobs$run))
todo <- jobs[!file.exists(jobs$file), ]; todo <- todo[order(match(todo$tool, c("scClassify", "CellTypist", "Seurat", "SingleR"))), ]
cat(nrow(jobs), "tool runs,", nrow(todo), "to compute\n")
invisible(mclapply(seq_len(nrow(todo)), function(j) { run <- runs[[todo$run[j]]]
  r <- tryCatch(tools[[todo$tool[j]]](make_counts(run), run$train, run$test), error = function(e) { message(todo$tool[j], " ", todo$run[j], " FAILED: ", conditionMessage(e)); NULL })
  if (!is.null(r)) saveRDS(r, todo$file[j]) }, mc.cores = n_cores, mc.preschedule = FALSE))
if (any(!file.exists(jobs$file))) stop("some tool runs failed; rerun the script to retry them")

## ---- 5. uncertainty of each cell ----
clean <- function(l) ifelse(l %in% types, l, "No label")                      # unassigned and in-between answers count as one outcome
pred <- function(tool, id) { r <- do.call(rbind, lapply(1:2, function(i) readRDS(file.path(cache, sprintf("%s__%s_%d.rds", tool, id, i))))); clean(r$label[match(colnames(so), r$cell)]) }
changed <- do.call(cbind, lapply(names(tools), function(t) { base <- pred(t, "split")
  do.call(cbind, lapply(names(perts), function(p) sapply(seq_len(n_rep), function(r) pred(t, sprintf("perturb_%s_r%d", p, r)) != base))) }))
n_pred <- ncol(changed)                                                         # 4 tools x 4 perturbations x 5 replicates = 80
x <- data.frame(cell = colnames(so), sample = smp, type = factor(lab, levels = types), umap1 = Embeddings(so, "umap")[, 1], umap2 = Embeddings(so, "umap")[, 2],
                uncertainty = 100 * rowMeans(changed))
write.csv(transform(x[, c("cell", "sample", "type", "umap1", "umap2", "uncertainty")], uncertainty = round(uncertainty, 4)), file.path(out_dir, "uncertainty_per_cell.csv"), row.names = FALSE)   # also read by supp_fig_proportion_vs_uncertainty.R
cat(sprintf("%d predictions per cell | cells never changing label: %d (%.1f%%) | mean uncertainty: %.2f%%\n", n_pred, sum(x$uncertainty == 0), 100 * mean(x$uncertainty == 0), mean(x$uncertainty)))

## ---- 6. figures ----
ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"
type_cols <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948", "#86b6ef", "#f4a582", "#7fd1b3", "#f6cf6a",
                        "#b05c84", "#66b266", "#9085e9", "#a83232", "#184f95", "#a84a1f", "#0f7a55", "#8a6400", "#6d6c66", "#b9b8b0"), types)
bs <- 7; txt <- (bs - 1) / ggplot2::.pt
th <- theme_minimal(base_size = bs) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(colour = ink2),
  axis.text = element_text(colour = ink, size = bs - 0.5), axis.title = element_text(colour = ink2, size = bs),
  plot.title = element_text(colour = ink, face = "bold", size = bs + 1, margin = margin(b = 2)), plot.title.position = "plot",
  plot.subtitle = element_text(colour = ink2, size = bs - 0.3, lineheight = 1.1, margin = margin(b = 4)), legend.text = element_text(size = bs - 0.3),
  legend.title = element_text(size = bs - 0.3), plot.margin = margin(5, 7, 5, 5), plot.background = element_rect(fill = "white", colour = NA))
box <- theme(panel.border = element_rect(colour = ink2, fill = NA, linewidth = 0.4))

# Figure 4C: uncertainty of each cell
fig4c <- ggplot(x[order(x$uncertainty), ], aes(umap1, umap2)) +              # most uncertain cells drawn last, on top
  geom_point(aes(colour = uncertainty), size = 0.3, stroke = 0) +
  scale_colour_viridis_c(option = "D", direction = -1, end = 0.92, limits = c(0, 30), oob = scales::squish, breaks = c(0, 15, 30), labels = c("0%", "15%", "30%+"),
                         name = "Uncertainty", guide = guide_colourbar(barheight = 0.4, barwidth = 5, ticks = FALSE, title.position = "left", title.vjust = 1)) +
  labs(title = "Uncertainty of each cell", x = "UMAP 1", y = "UMAP 2") +
  th + theme(axis.text = element_blank(), legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, 0, 0, 0))

# Figure 4D: uncertainty by cell type
d <- x %>% group_by(type) %>% summarise(n = n(), mean = mean(uncertainty), .groups = "drop") %>% arrange(mean)
cat("\nmean uncertainty by cell type (%):\n"); print(as.data.frame(d[order(-d$mean), c("type", "n", "mean")]), row.names = FALSE, digits = 3)
d$type <- factor(d$type, levels = as.character(d$type))
fig4d <- ggplot(d, aes(mean, type)) + geom_col(aes(fill = type), width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", mean)), hjust = -0.2, size = txt * 1.1, colour = ink) +
  scale_fill_manual(values = type_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28)), labels = function(v) paste0(v, "%")) + scale_y_discrete(expand = expansion(add = 0.5)) +
  labs(title = "Uncertainty by cell type", x = "Mean uncertainty", y = NULL) +
  th + box + theme(axis.text.y = element_text(size = bs + 0.5))

ggsave(file.path(out_dir, "fig4c.pdf"), fig4c, width = 3.6, height = 3.8, useDingbats = FALSE)
ggsave(file.path(out_dir, "fig4d.pdf"), fig4d, width = 3.4, height = 3.6, useDingbats = FALSE)
cat("\nwrote fig4c.pdf, fig4d.pdf and uncertainty_per_cell.csv in", out_dir, "\n")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
