################################################################################
# 00.5_LD_decay.R  -  CENTROID-FILTERED PIPELINE (run after 00, before 06)
#
# Linkage disequilibrium decay with the four LDcorSV measures (Mangin et al. 2012):
#   r2    usual squared correlation between genotypes
#   r2s   corrected for population structure (S = ancestry proportions, K-1 columns)
#   r2v   corrected for relatedness (V = kinship / genomic relationship matrix)
#   r2vs  corrected for both
# for all SNP pairs on the same chromosome. The decay of each measure with physical
# distance is fitted with the Hill & Weir (1988) expectation, and the distance at which
# r2 falls to 0.2 is reported per chromosome, per sub-genome (A, B) and genome-wide.
#
# Implementation: the measures are computed in matrix form for all pairs of a chromosome
# at once. They are identical to LDcorSV's Measure.R2 / R2S / R2V / R2VS (checked to
# ~1e-15) but thousands of times faster than LDcorSV::LD.Measures, which loops over pairs.
# The kinship matrix is the VanRaden genomic relationship matrix plus a small ridge
# (CFG$kinship_ridge x mean diagonal), because a mean-centred relationship matrix is
# singular in the direction of the overall mean, which makes the r2v / r2vs mean undefined.
# Missing genotypes are imputed by the rounded SNP mean (LDcorSV drops them pairwise).
#
# Inputs : data_filtered/finger_millet_filtered.vcf (from 00); population structure from the
#          TESS3 Q matrix of step 03 if it exists (structure = "auto"), otherwise LEA::snmf.
# Progress is written to results/00.5_LD_decay/LD_decay_progress.log after every step, so if R
# crashes the last line shows where. Running from a terminal (Rscript 00.5_LD_decay.R) shows the
# full error message instead of an RStudio "session aborted" box.
# Check: run on the full 404-accession VCF, this script gives 2.29 M pairs and r2_VS = 0.2 at
#   82.4 kb (genome), 67.6 kb (A), 134.1 kb (B), against 2.32 M pairs and 78.7 / 64.7 / 127.8 kb
#   in the original LDcorSV analysis.
# Outputs: reference_data/LD_decay_summary.csv   <- read by 06_candidate_genes_GO.R
#          results/00.5_LD_decay/LD_decay_summary.csv, LD_decay_binned.csv,
#          structure_Q.csv, figures/*.pdf + 600-dpi *.png
################################################################################

CFG <- list(
  vcf            = "data_filtered/finger_millet_filtered.vcf",
  chr_map        = setNames(paste0(rep(1:9, each = 2), c("A", "B")), 1:18),  # TASSEL 1-18 -> 1A..9B; 19 = scaffolds (excluded)
  maf            = 0.01,        # as in the full-panel LD analysis (MAF 0.01, all SNPs): these settings
  max_snp_missing = 1,          # reproduce its pair count and r2_VS distances to within ~5% on the 404-accession VCF
  max_distance   = Inf,         # bp; Inf = all pairs on a chromosome (as in the full-panel analysis)
  max_snps_per_chr = Inf,       # subsample SNPs per chromosome for speed (Inf = all)
  structure      = "auto",      # "auto" = use q_file if it exists, else sNMF; or "snmf", "file", "none"
  K              = 4,           # ancestral populations for snmf
  snmf_reps      = 5,
  q_file         = "results/03_TESS3/tables/TESS3_Qmatrix_K4.csv",  # ancestry from step 03 (run 03 first to skip sNMF);
  # any CSV with the accession ID in column 1 and columns Q1, Q2, ...
  kinship_ridge  = 0.01,
  r2_threshold   = 0.2,
  measures       = c("r2", "r2s", "r2v", "r2vs"),
  summary_to     = "reference_data/LD_decay_summary.csv",   # copy used by script 06
  out_dir        = "results/00.5_LD_decay",
  dpi            = 600,
  seed           = 2026
)

need <- function(p) { m <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
if (length(m)) install.packages(m, repos = "https://cloud.r-project.org") }
need(c("data.table", "ggplot2", "MASS"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(CFG$seed)
dir.create(file.path(CFG$out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(CFG$summary_to), recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(CFG$out_dir, "LD_decay_progress.log"); cat("", file = LOG)
msg <- function(...) {                             # console + log file written immediately (survives a crash)
  line <- paste(format(Sys.time(), "%H:%M:%S"), paste(..., collapse = " "))
  cat(line, "\n"); cat(line, "\n", file = LOG, append = TRUE) }
save_fig <- function(name, p, w, h) {
  ggsave(file.path(CFG$out_dir, "figures", paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(CFG$out_dir, "figures", paste0(name, ".png")), p, width = w, height = h, dpi = CFG$dpi, bg = "white")
}

# ---- genotypes ---------------------------------------------------------------------------------
L <- readLines(CFG$vcf); hdr <- L[startsWith(L, "#CHROM")]; L <- L[!startsWith(L, "#")]
samples <- strsplit(hdr, "\t")[[1]][-(1:9)]
v <- do.call(rbind, strsplit(L, "\t", fixed = TRUE)); rm(L)
v <- v[!grepl(",", v[, 5]), , drop = FALSE]
gt <- substr(v[, -(1:9), drop = FALSE], 1, 3)
X <- (substr(gt, 1, 1) != "0") + (substr(gt, 3, 3) != "0"); X[substr(gt, 1, 1) == "."] <- NA
X <- t(X); storage.mode(X) <- "double"; rownames(X) <- samples
snp <- data.frame(chrom = ifelse(v[, 1] %in% names(CFG$chr_map), CFG$chr_map[v[, 1]], "Un"), pos = as.numeric(v[, 2]))
p <- colMeans(X, na.rm = TRUE) / 2
keep <- colMeans(is.na(X)) <= CFG$max_snp_missing & pmin(p, 1 - p) >= CFG$maf & snp$chrom != "Un"
X <- X[, keep]; snp <- snp[keep, ]
for (j in which(colSums(is.na(X)) > 0)) { x <- X[, j]; x[is.na(x)] <- round(mean(x, na.rm = TRUE)); X[, j] <- x }
msg("Accessions:", nrow(X), "| SNPs on chromosomes 1A-9B after filters:", ncol(X))

msg("Genotypes read and filtered")
# ---- kinship V (VanRaden GRM + ridge) and structure S --------------------------------------------
pp <- colMeans(X) / 2; Z <- sweep(X, 2, 2 * pp)
V <- tcrossprod(Z) / (2 * sum(pp * (1 - pp))); V <- V + diag(CFG$kinship_ridge * mean(diag(V)), nrow(V))
S <- NULL
if (CFG$structure == "auto") CFG$structure <- if (!is.na(CFG$q_file) && file.exists(CFG$q_file)) "file" else "snmf"
if (CFG$structure == "file" && (is.na(CFG$q_file) || !file.exists(CFG$q_file))) stop("q_file not found: ", CFG$q_file)
msg("Structure source:", CFG$structure)
if (CFG$structure == "snmf") {
  msg("Running sNMF (LEA) in a separate R process")
  if (!requireNamespace("LEA", quietly = TRUE)) { need("BiocManager"); BiocManager::install("LEA", update = FALSE, ask = FALSE) }
  # LEA's C code can crash the whole R session (e.g. if it cannot open a file), so sNMF runs in
  # its own Rscript process: a crash there stops this script with a message instead of aborting R.
  # The input file sits inside the working directory, because LEA builds its project path
  # relative to getwd() and fails for files elsewhere (e.g. tempdir()).
  sdir <- file.path(CFG$out_dir, "structure_sNMF"); dir.create(sdir, showWarnings = FALSE)
  gf <- file.path(sdir, "ld_genotypes.geno"); LEA::write.geno(X, gf)
  qout <- file.path(sdir, "Q_sNMF.csv"); unlink(qout)
  job <- file.path(sdir, "run_sNMF.R")
  writeLines(c(
    sprintf('invisible(capture.output(sn <- LEA::snmf("%s", K = %d, repetitions = %d, entropy = TRUE, project = "new", seed = %d, CPU = 1)))',
            gf, CFG$K, CFG$snmf_reps, CFG$seed),
    sprintf('best <- which.min(LEA::cross.entropy(sn, K = %d)); Q <- LEA::Q(sn, K = %d, run = best)', CFG$K, CFG$K),
    sprintf('write.csv(Q, "%s", row.names = FALSE)', qout)), job)
  status <- system2(file.path(R.home("bin"), "Rscript"), job, stdout = file.path(sdir, "sNMF_log.txt"), stderr = file.path(sdir, "sNMF_log.txt"))
  if (!file.exists(qout))
    stop("sNMF failed (exit status ", status, "; see ", file.path(sdir, "sNMF_log.txt"), "). ",
         "Run step 03 first and keep structure = \"auto\" (uses the TESS3 Q matrix), or set structure = \"none\".")
  Q <- as.matrix(read.csv(qout)); S <- Q[, -ncol(Q), drop = FALSE]
  msg("Structure: sNMF K =", CFG$K, "(best of", CFG$snmf_reps, "runs)")
} else if (CFG$structure == "file") {
  qd <- read.csv(CFG$q_file, check.names = FALSE); qd <- qd[match(rownames(X), trimws(qd[[1]])), ]
  if (anyNA(qd[[1]])) stop("Q file lacks some accessions")
  qc <- grep("^Q[0-9]+$|^V[0-9]+$|^cluster_?[0-9]+$", names(qd), value = TRUE, ignore.case = TRUE)
  if (length(qc) < 2) stop("q_file needs ancestry columns named Q1, Q2, ... (found: ", paste(names(qd), collapse = ", "), ")")
  Q <- as.matrix(qd[, qc]); S <- Q[, -ncol(Q), drop = FALSE]
  msg("Structure from", CFG$q_file)
}
if (!is.null(S)) { rownames(S) <- rownames(X)
write.csv(data.frame(SampleID = rownames(X), Q), file.path(CFG$out_dir, "structure_Q.csv"), row.names = FALSE) }

# ---- LD measures, vectorised (identical to LDcorSV Measure.R2 / R2S / R2V / R2VS) ------------------
pinv_sqrt <- function(M, tol = 1e-5) { e <- eigen(M, symmetric = TRUE); d <- ifelse(e$values < tol, 0, 1 / sqrt(e$values))
e$vectors %*% (d * t(e$vectors)) }
cos_sq <- function(M, eps = 1e-7) { G <- crossprod(M); d <- diag(G); R <- G^2 / outer(d, d)
bad <- d < eps; R[bad, ] <- 0; R[, bad] <- 0; R }
msg("Computing kinship inverse square root")
n <- nrow(X); one <- rep(1, n); W <- pinv_sqrt(V); Vi1 <- W %*% (W %*% one)
gls_center <- function(M) M - one %*% (crossprod(Vi1, M) / sum(one * Vi1))
if (!is.null(S)) { Sc <- scale(S, scale = FALSE); Zs <- W %*% gls_center(S) }
ld_chrom <- function(idx) {
  Xi <- X[, idx, drop = FALSE]; Xc <- scale(Xi, scale = FALSE); out <- list(r2 = cos_sq(Xc))
  Zi <- W %*% gls_center(Xi); out$r2v <- cos_sq(Zi)
  if (!is.null(S)) {
    out$r2s  <- cos_sq(Xc - Sc %*% MASS::ginv(crossprod(Sc)) %*% crossprod(Sc, Xc))
    out$r2vs <- cos_sq(Zi - Zs %*% MASS::ginv(crossprod(Zs)) %*% crossprod(Zs, Zi))
  }
  out
}

pairs <- list()
for (ch in unique(snp$chrom)) {
  idx <- which(snp$chrom == ch)
  if (length(idx) > CFG$max_snps_per_chr) idx <- sort(sample(idx, CFG$max_snps_per_chr))
  if (length(idx) < 3) next
  m <- ld_chrom(idx); ut <- which(upper.tri(m$r2), arr.ind = TRUE)
  d <- abs(snp$pos[idx][ut[, 2]] - snp$pos[idx][ut[, 1]])
  dt <- data.table(chrom = ch, dist = d)
  for (k in intersect(CFG$measures, names(m))) dt[[k]] <- m[[k]][ut]
  pairs[[ch]] <- dt[dist > 0 & dist <= CFG$max_distance]
  msg("  ", ch, ":", length(idx), "SNPs,", nrow(pairs[[ch]]), "pairs")
}
pairs <- rbindlist(pairs); pairs[, subgenome := substring(chrom, nchar(chrom))]
msg("Pairs:", nrow(pairs), "| approx. memory", format(object.size(pairs), units = "MB")); invisible(gc())
meas <- intersect(CFG$measures, names(pairs))

# ---- Hill & Weir (1988) decay fit and r2 = threshold distance -------------------------------------
hw <- function(d, C, n) { x <- C * d
((10 + x) / ((2 + x) * (11 + x))) * (1 + ((3 + x) * (12 + 12 * x + x^2)) / (n * (2 + x) * (11 + x))) }
fit_hw <- function(d, r) {
  # least squares in log10(C): coarse grid, then optimize() around the best grid point (pure R)
  sse <- function(lc) sum((r - hw(d, 10^lc, n))^2)
  grid <- seq(-8, 1, by = 0.25); g <- vapply(grid, sse, 0); b <- grid[which.min(g)]
  10^optimize(sse, c(b - 0.25, b + 0.25), tol = 1e-6)$minimum
}
thr_dist <- function(C, thr, dmax) {
  if (hw(0, C, n) <= thr) return(c(0, FALSE))
  up <- 1e9; if (hw(up, C, n) > thr) return(c(NA, TRUE))
  d <- uniroot(function(x) hw(x, C, n) - thr, c(0, up), tol = 1)$root
  c(d, d > dmax) }
chr_order <- unname(CFG$chr_map)
scopes <- c(list("Whole genome" = seq_len(nrow(pairs))),               # row indices only - no copies
            split(seq_len(nrow(pairs)), paste("Sub-genome", pairs$subgenome)),
            split(seq_len(nrow(pairs)), pairs$chrom))
scopes <- scopes[c("Whole genome", "Sub-genome A", "Sub-genome B", intersect(chr_order, names(scopes)))]
msg("Fitting Hill & Weir decay curves:", length(scopes), "scopes x", length(meas), "measures")
res <- rbindlist(lapply(names(scopes), function(sc) {
  ix <- scopes[[sc]]; dd <- list(dist = pairs$dist[ix])
  rbindlist(lapply(meas, function(mm) {
    dd[[mm]] <- pairs[[mm]][ix]
    C <- fit_hw(dd$dist, dd[[mm]]); td <- thr_dist(C, CFG$r2_threshold, max(dd$dist))
    hd <- tryCatch(uniroot(function(x) hw(x, C, n) - hw(0, C, n) / 2, c(0, 1e9), tol = 1)$root, error = function(e) NA)
    data.table(scope = sc, measure = mm, n_pairs = length(ix), C = C, d_threshold_kb = td[1] / 1000,
               extrapolated = as.logical(td[2]), r2_fit_at_0 = hw(0, C, n), d_half_decay_kb = hd / 1000,
               mean_r2 = mean(dd[[mm]]))
  }))
}))
msg("Fits done")
fwrite(res, file.path(CFG$out_dir, "LD_decay_summary.csv"))
fwrite(res, CFG$summary_to)                       # for 06_candidate_genes_GO.R
print(res[scope %in% c("Whole genome", "Sub-genome A", "Sub-genome B"), .(scope, measure, d_threshold_kb = round(d_threshold_kb, 1), extrapolated)])
msg("LD summary ->", CFG$summary_to, "and", file.path(CFG$out_dir, "LD_decay_summary.csv"))

# ---- binned means and figures ----------------------------------------------------------------------
brk <- unique(c(0, 10^seq(2, log10(max(pairs$dist)) + 0.1, by = 0.1)))
pairs[, bin := cut(dist, brk, include.lowest = TRUE)]
invisible(gc())
msg("Binning and plotting")
binned <- rbindlist(lapply(names(scopes)[1:3], function(sc) {
  dd <- pairs[scopes[[sc]]]
  b <- dd[, c(list(dist_kb = median(dist) / 1000, n_pairs = .N), lapply(.SD, mean)), by = bin, .SDcols = meas]
  b[, scope := sc] }))
fwrite(binned, file.path(CFG$out_dir, "LD_decay_binned.csv"))
bl <- melt(binned, id.vars = c("scope", "bin", "dist_kb", "n_pairs"), measure.vars = meas, variable.name = "measure", value.name = "r2")
curves <- rbindlist(lapply(seq_len(nrow(res)), function(i) if (res$scope[i] %in% names(scopes)[1:3]) {
  dk <- 10^seq(-1, log10(max(pairs$dist) / 1000), length.out = 200)
  data.table(scope = res$scope[i], measure = res$measure[i], dist_kb = dk, r2 = hw(dk * 1000, res$C[i], n)) }))
lab <- c(r2 = "r\u00b2", r2s = "r\u00b2 S (structure)", r2v = "r\u00b2 V (relatedness)", r2vs = "r\u00b2 VS (both)")
p1 <- ggplot(bl, aes(dist_kb, r2, colour = measure)) + geom_point(size = 1, alpha = .6) +
  geom_line(data = curves, linewidth = .8) + geom_hline(yintercept = CFG$r2_threshold, linetype = 2, colour = "grey40") +
  facet_wrap(~ scope, nrow = 1) + scale_x_log10() + scale_colour_brewer(palette = "Dark2", labels = lab, name = NULL) +
  labs(x = "Distance between SNPs (kb, log scale)", y = expression(LD ~ (r^2)),
       title = "Linkage disequilibrium decay (points = binned means; lines = Hill & Weir fit)",
       subtitle = paste0("Dashed line: r\u00b2 = ", CFG$r2_threshold, "; ", nrow(X), " accessions, ", ncol(X), " SNPs")) +
  theme_bw() + theme(legend.position = "bottom")
save_fig("LD_decay_curves", p1, 13, 4.8)
rc <- res[!scope %in% c("Whole genome", "Sub-genome A", "Sub-genome B")]
rc[, scope := factor(scope, levels = chr_order)]
p2 <- ggplot(rc, aes(scope, d_threshold_kb, fill = measure)) + geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Dark2", labels = lab, name = NULL) + scale_y_log10() +
  labs(x = "Chromosome", y = paste0("Distance at r\u00b2 = ", CFG$r2_threshold, " (kb, log scale)"),
       title = "LD decay distance by chromosome") + theme_bw() + theme(legend.position = "bottom")
save_fig("LD_decay_by_chromosome", p2, 11, 5)
msg("Done ->", normalizePath(CFG$out_dir))