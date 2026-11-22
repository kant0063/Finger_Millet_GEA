################################################################################
# finger_millet_pRDA_FST_LFMM.R
#
# Finger millet (Eleusine coracana) version of the yam scripts
#   05_partial_RDA_varpart.R   (Analysis Plan A5)  -> Section 2
#   06_Fst_outliers_revised.R  (Analysis Plan A6c) -> Section 3
#   09_LFMM.R                  (Analysis Plan A6e) -> Section 4
# combined into one standalone script (no 00_config_and_helpers.R needed). Methods and
# outputs are the same as in the yam scripts but with descriptive (unnumbered) file names;
# figures are written as PDF + PNG.
#
# Section 2  Partial RDA and variance partitioning (vegan)
#   Environment conditioned on ancestry (genetic PCs) and geography (lat, lon); VIF-pruned
#   variables; varpart with adjusted R2; 999 permutations by model, axis and term. RDA is
#   fitted to the full set of genetic PC scores (identical results to the SNP matrix,
#   far less memory). Run for the full panel and for each subset in CFG$subsets.
#   -> pRDA_summary_by_set.csv, pRDA_variance_partitioning_fractions.csv, pRDA_permutation_tests_by_axis.csv,
#      pRDA_permutation_tests_by_variable.csv, RDA_full_model_environment_geography_ancestry, RDA_partial_environment_given_ancestry_geography, rda_models.rds
# Section 3  Weir & Cockerham (1984) FST outliers
#   Per-locus WC FST computed directly from genotypes (not hierfstat::betas on dosages):
#   (1) global across the genetic groups, (2) pairwise among non-admixed groups. Private
#   alleles flagged. Outliers at the 95th / 99th / 99.5th percentiles.
#   -> FST_pairwise_between_groups_genomewide.csv, FST_outlier_loci_global.csv, FST_outlier_counts_by_chromosome.csv,
#      fst_per_locus_WC.csv, FST_per_locus_global_by_chromosome
# Section 4  LFMM2 genotype-environment association (LEA)
#   Per variable, genomic control, BH FDR; overlap with the eGWAS hits
#   (significant_markers.csv from environmental_gwas_complete.R).
#   -> LFMM2_significant_SNP_environment_associations.csv, lfmm_<set>.rds, geno_<set>.lfmm
#
# Finger millet specifics
#   - Groups: the Population column of the coordinate file (or Region if there is no
#     Population column), e.g. East Africa, Southern Africa, India, Nepal (Bancic et al. 2024).
#   - Subsets replace the yam species subsets: CFG$subsets (default Africa / Asia, the
#     primary genetic split in finger millet). Ancestry in every set = genetic PCs.
#   - Non-admixed samples (pairwise FST): a Q-matrix CSV if given (CFG$qmatrix_csv), else
#     LEA::snmf with K = number of groups; max ancestry >= CFG$nonadmixed_q and the sample's
#     main cluster = its group's main cluster.
#   - Chromosomes: TASSEL numbering 1-18 -> 1A,1B,...,9A,9B; 19 = unplaced scaffolds ("Un").
#   - Environment: CSV of SampleID + environmental variables (e.g. bioclim extracted at the
#     sites, as used in the eGWAS / genomic prediction scripts). VIF-pruned.
################################################################################

# ==============================================================================
# 0  CONFIGURATION
# ==============================================================================
CFG <- list(
  vcf          = "finger_millet_404.vcf",
  coords_csv   = "millet_passport_test.csv",          # SampleID, Latitude, Longitude, Population (or Region)
  env_csv      = "millet_extracted_bioclim_values.csv",  # SampleID + environmental variables
  egwas_csv    = "output/significant_markers.csv",    # eGWAS hits (optional; for "shared_with_EGWAS")
  qmatrix_csv  = NA,                                  # optional: SampleID + ancestry proportions (Q1..QK)
  env_vars     = NULL,                                # NULL = all numeric columns; or e.g. paste0("bio", c(1,4,12,15))
  subsets      = list(Africa = c("East Africa", "Southern Africa"), Asia = c("India", "Nepal")),
  min_subset_n = 30,                                  # subsets smaller than this are skipped
  chr_map      = setNames(paste0(rep(1:9, each = 2), c("A", "B")), 1:18),  # VCF CHROM -> chromosome
  # ---- genotype filters ----
  maf          = 0.05,
  max_snp_missing = 0.20,
  thin_bp      = 1000,                                # RDA: keep 1 SNP per thin_bp (0 = no thinning)
  # ---- RDA ----
  vif_max      = 10,
  n_pc         = 3,                                   # genetic PCs used as ancestry covariates
  perms        = 999,
  rda_max_snps = 0,                                   # 0 = all thinned SNPs
  rda_cores    = 1,
  # ---- FST ----
  nonadmixed_q = 0.8,
  # ---- LFMM2 ----
  K_lfmm_full  = 4,                                   # full panel (4 genetic clusters in Bancic et al. 2024)
  K_lfmm_within = list(Africa = 2, Asia = 2),
  fdr          = 0.05,
  # ---- output ----
  out_dir      = "pRDA_FST_LFMM_results",
  seed         = 2026
)

# ==============================================================================
# 1  PACKAGES, HELPERS, DATA
# ==============================================================================
need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
}
need(c("vegan", "dplyr", "ggplot2"))
if (!requireNamespace("LEA", quietly = TRUE)) { need("BiocManager"); BiocManager::install("LEA", update = FALSE, ask = FALSE) }
suppressPackageStartupMessages({ library(vegan); library(dplyr); library(ggplot2) })
set.seed(CFG$seed)
for (d in c("tables", "figures", "results")) dir.create(file.path(CFG$out_dir, d), recursive = TRUE, showWarnings = FALSE)
tab <- function(f) file.path(CFG$out_dir, "tables", f)
fig <- function(f) file.path(CFG$out_dir, "figures", f)
res <- function(f) file.path(CFG$out_dir, "results", f)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
save_fig <- function(name, p, width, height) {                  # PDF (vector) + PNG
  ggsave(fig(paste0(name, ".pdf")), p, width = width, height = height,
         device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf)
  ggsave(fig(paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
}

# ---- genotypes: individuals x SNPs, ALT dosage 0/1/2, NA = missing ----
read_vcf_dosage <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
  n_skip <- 0
  repeat { l <- readLines(con, n = 1); if (!length(l)) { close(con); stop("No #CHROM header") }
  if (startsWith(l, "#CHROM")) break; n_skip <- n_skip + 1 }
  close(con)
  smp <- strsplit(l, "\t")[[1]][-(1:9)]
  v <- read.delim(path, skip = n_skip + 1, header = FALSE, colClasses = "character", comment.char = "", quote = "")
  v <- v[!grepl(",", v[[5]]), , drop = FALSE]
  gt <- sub(":.*$", "", as.matrix(v[, -(1:9), drop = FALSE]))
  a1 <- sub("[/|].*$", "", gt); a2 <- sub("^[^/|]*[/|]", "", gt)
  d <- (a1 != "0") + (a2 != "0"); d[a1 == "." | a2 == "."] <- NA; storage.mode(d) <- "numeric"
  id <- make.unique(ifelse(v[[3]] %in% c(".", ""), paste(v[[1]], v[[2]], sep = "_"), v[[3]]))
  G <- t(d); rownames(G) <- trimws(smp); colnames(G) <- id
  list(G = G, snp = data.frame(SNP = id, CHROM_vcf = v[[1]], POS = as.integer(v[[2]]), stringsAsFactors = FALSE))
}
msg("Reading", CFG$vcf)
vc <- read_vcf_dosage(CFG$vcf)
snp <- vc$snp
snp$CHROM <- ifelse(snp$CHROM_vcf %in% names(CFG$chr_map), CFG$chr_map[snp$CHROM_vcf], "Un")
chr_levels <- c(unname(CFG$chr_map), "Un")

# ---- samples: coordinates, groups, environment ----
co <- read.csv(CFG$coords_csv, stringsAsFactors = FALSE, check.names = FALSE); nm <- tolower(names(co))
gcol <- c(grep("^pop", nm), grep("^region", nm), grep("^group", nm))[1]
if (is.na(gcol)) stop("Coordinate file needs a Population (or Region) column")
meta <- data.frame(SampleID = trimws(co[[1]]), Latitude = as.numeric(co[[grep("^lat", nm)[1]]]),
                   Longitude = as.numeric(co[[grep("^lon", nm)[1]]]), Population = as.character(co[[gcol]]),
                   stringsAsFactors = FALSE)
meta <- meta[!duplicated(meta$SampleID), ]
env <- read.csv(CFG$env_csv, stringsAsFactors = FALSE, check.names = FALSE); names(env)[1] <- "SampleID"
env$SampleID <- trimws(env$SampleID)
ev <- setdiff(names(env)[vapply(env, is.numeric, logical(1))], c("SampleID", grep("^lat|^lon", names(env), value = TRUE, ignore.case = TRUE)))
if (!is.null(CFG$env_vars)) ev <- intersect(CFG$env_vars, ev)
env <- env[complete.cases(env[, ev]), c("SampleID", ev)]
ids <- Reduce(intersect, list(rownames(vc$G), meta$SampleID[!is.na(meta$Latitude) & !is.na(meta$Population)], env$SampleID))
msg("Samples: VCF", nrow(vc$G), "| coordinates", nrow(meta), "| environment", nrow(env), "| all three", length(ids))
if (length(ids) < 30) stop("Fewer than 30 samples shared by VCF, coordinates and environment")
meta <- meta[match(ids, meta$SampleID), ]
G <- vc$G[ids, , drop = FALSE]
E_all <- as.matrix(env[match(ids, env$SampleID), ev, drop = FALSE]); rownames(E_all) <- ids
E_all <- E_all[, apply(E_all, 2, sd) > 0, drop = FALSE]

# SNP QC on the full panel (missingness, MAF); missing genotypes stay NA here
p <- colMeans(G, na.rm = TRUE) / 2
keep <- !is.na(p) & colMeans(is.na(G)) <= CFG$max_snp_missing & pmin(p, 1 - p) >= CFG$maf
G <- G[, keep, drop = FALSE]; snp <- snp[match(colnames(G), snp$SNP), ]
msg("SNPs after QC:", ncol(G), "| groups:", paste(names(table(meta$Population)), table(meta$Population), sep = "=", collapse = ", "))

# MAF filter within a subset and impute remaining missing calls (rounded mean): LFMM/RDA input
maf_filter <- function(X, maf = CFG$maf) {
  pp <- colMeans(X, na.rm = TRUE) / 2; X <- X[, !is.na(pp) & pmin(pp, 1 - pp) >= maf, drop = FALSE]
  for (j in which(colSums(is.na(X)) > 0)) { x <- X[, j]; x[is.na(x)] <- round(mean(x, na.rm = TRUE)); X[, j] <- x }
  X
}
# stepwise VIF pruning (drop the highest VIF until all < CFG$vif_max)
vif_prune <- function(E, thr = CFG$vif_max) {
  v <- colnames(E); vifs <- NULL
  repeat {
    if (length(v) < 2) break
    vifs <- sapply(v, function(x) 1 / (1 - summary(lm(E[[x]] ~ ., data = E[, setdiff(v, x), drop = FALSE]))$r.squared))
    if (max(vifs) < thr) break
    v <- setdiff(v, names(which.max(vifs)))
  }
  list(keep = v, vif = vifs)
}
vif_keep <- vif_prune(as.data.frame(scale(E_all)))$keep
msg("Environmental variables after VIF <", CFG$vif_max, ":", paste(vif_keep, collapse = ", "))
write.csv(data.frame(variable = colnames(E_all), kept_after_VIF = colnames(E_all) %in% vif_keep), tab("environmental_variables_VIF.csv"), row.names = FALSE)

# thinned SNP set for RDA (1 SNP per CFG$thin_bp per chromosome)
thin_snps <- function(s, bp) {
  if (bp <= 0) return(s$SNP)
  s <- s[order(s$CHROM_vcf, s$POS), ]
  unlist(lapply(split(s, s$CHROM_vcf), function(d) { kept <- d$SNP[1]; last <- d$POS[1]
  for (i in seq_len(nrow(d))[-1]) if (d$POS[i] - last >= bp) { kept <- c(kept, d$SNP[i]); last <- d$POS[i] }
  kept }))
}
thin <- thin_snps(snp, CFG$thin_bp)
msg("SNPs after thinning to 1 per", CFG$thin_bp, "bp (RDA):", length(thin))

# subsets: full panel + CFG$subsets (by group membership)
S <- list(full = ids)
for (s in names(CFG$subsets)) {
  i <- ids[meta$Population %in% CFG$subsets[[s]]]
  if (length(i) >= CFG$min_subset_n) S[[s]] <- i else msg("Subset", s, "skipped: only", length(i), "samples")
}
msg("Sets analysed:", paste(names(S), lengths(S), sep = "=", collapse = ", "))

# ==============================================================================
# 2  PARTIAL RDA AND VARIANCE PARTITIONING   (yam 05_partial_RDA_varpart.R)
# ==============================================================================
run_rda <- function(set_ids, label) {
  m <- meta[match(set_ids, meta$SampleID), ]
  Y <- maf_filter(G[set_ids, intersect(thin, colnames(G)), drop = FALSE])
  if (CFG$rda_max_snps > 0 && ncol(Y) > CFG$rda_max_snps) {
    set.seed(CFG$seed); Y <- Y[, sort(sample(ncol(Y), CFG$rda_max_snps)), drop = FALSE] }
  storage.mode(Y) <- "double"; n_snps <- ncol(Y)
  # RDA depends on Y only through its Euclidean geometry: the full set of PC scores gives
  # identical eigenvalues, R2, adj. R2, varpart, F, permutation P and site scores.
  pcY <- prcomp(Y, center = TRUE, scale. = FALSE); k <- pcY$sdev > 1e-8 * max(pcY$sdev)
  Y <- pcY$x[, k, drop = FALSE]; rm(pcY); invisible(gc())
  E <- as.data.frame(scale(E_all[set_ids, vif_keep, drop = FALSE])); vp <- vif_prune(E); E <- E[, vp$keep, drop = FALSE]
  Sg <- data.frame(lat = scale(m$Latitude)[, 1], lon = scale(m$Longitude)[, 1])
  A <- as.data.frame(Y[, 1:min(CFG$n_pc, ncol(Y) - 1), drop = FALSE])          # ancestry = genetic PCs
  names(A) <- paste0("genPC", seq_len(ncol(A)))
  A <- A[, vapply(A, function(v) sd(v) > 0, logical(1)), drop = FALSE]; if (!ncol(A)) A <- NULL
  D <- if (is.null(A)) cbind(E, Sg) else cbind(E, Sg, A)
  cond <- paste(c(names(Sg), if (!is.null(A)) names(A)), collapse = " + ")
  full <- rda(as.formula(paste("Y ~", paste(names(D), collapse = " + "))), data = D)
  penv <- rda(as.formula(paste("Y ~", paste(names(E), collapse = " + "), "+ Condition(", cond, ")")), data = D)
  vpart <- if (!is.null(A)) varpart(Y, E, A, Sg) else varpart(Y, E, Sg)
  msg(label, ": permutation tests (", CFG$perms, ")...")
  a_mod  <- anova.cca(penv, permutations = CFG$perms, parallel = CFG$rda_cores)
  a_axis <- anova.cca(penv, by = "axis", permutations = CFG$perms, parallel = CFG$rda_cores)
  a_term <- anova.cca(penv, by = "margin", permutations = CFG$perms, parallel = CFG$rda_cores)
  fr <- vpart$part$indfract
  list(full = full, penv = penv, vpart = vpart, E = E, ids = set_ids,
       summary = data.frame(set = label, n = length(set_ids), n_SNPs = n_snps, n_env = ncol(E),
                            env_vars = paste(names(E), collapse = "; "),
                            full_R2 = RsquareAdj(full)$r.squared, full_adjR2 = RsquareAdj(full)$adj.r.squared,
                            env_given_geo_anc_adjR2 = RsquareAdj(penv)$adj.r.squared, env_given_geo_anc_P = a_mod$`Pr(>F)`[1],
                            pure_env = fr$Adj.R.square[1], pure_anc = if (!is.null(A)) fr$Adj.R.square[2] else NA,
                            pure_geo = fr$Adj.R.square[if (!is.null(A)) 3 else 2],
                            n_sig_axes = sum(a_axis$`Pr(>F)` < 0.05, na.rm = TRUE),
                            sig_terms = paste(rownames(a_term)[which(a_term$`Pr(>F)` < 0.05)], collapse = "; ")),
       axis = data.frame(set = label, a_axis), term = data.frame(set = label, a_term),
       fractions = data.frame(set = label, fraction = rownames(fr), fr))
}
msg("== 2. Partial RDA ==")
out <- lapply(names(S), function(s) run_rda(S[[s]], s)); names(out) <- names(S)
write.csv(bind_rows(lapply(out, `[[`, "summary")), tab("pRDA_summary_by_set.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(out, `[[`, "fractions")), tab("pRDA_variance_partitioning_fractions.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(out, `[[`, "axis")), tab("pRDA_permutation_tests_by_axis.csv"), row.names = TRUE)
write.csv(bind_rows(lapply(out, `[[`, "term")), tab("pRDA_permutation_tests_by_variable.csv"), row.names = TRUE)
saveRDS(out, res("rda_models.rds")); print(bind_rows(lapply(out, `[[`, "summary")) %>% select(-env_vars))

# ---- Figures: full RDA vs partial RDA (environment | ancestry + geography) ----
plot_rda <- function(mod, set_ids, title) {
  ve <- round(100 * mod$CCA$eig[1:2] / mod$tot.chi, 1)
  si <- as.data.frame(scores(mod, 1:2, "sites")); si$group <- meta$Population[match(set_ids, meta$SampleID)]
  bp <- as.data.frame(scores(mod, 1:2, "bp")); sf <- 0.8 * max(abs(si[, 1:2])) / max(abs(bp)); bp <- bp * sf
  bp$var <- rownames(bp)
  ggplot(si, aes(RDA1, RDA2, colour = group)) + geom_point(alpha = .8) +
    geom_segment(data = bp, aes(0, 0, xend = RDA1, yend = RDA2), inherit.aes = FALSE, colour = "red3",
                 arrow = arrow(length = unit(.15, "cm"))) +
    geom_text(data = bp, aes(RDA1, RDA2, label = var), inherit.aes = FALSE, size = 2.6, vjust = -0.4) +
    labs(title = title, x = paste0("RDA1 (", ve[1], "% of total variance)"), y = paste0("RDA2 (", ve[2], "% of total variance)"),
         colour = "Group") + theme_bw()
}
p_full <- plot_rda(out$full$full, out$full$ids, "Full RDA: environment + geography + ancestry")
p_part <- plot_rda(out$full$penv, out$full$ids, "Partial RDA: environment conditioned on ancestry and geography")
save_fig("RDA_full_model_environment_geography_ancestry", p_full, 8, 6)
save_fig("RDA_partial_environment_given_ancestry_geography", p_part, 8, 6)
# Axis % are % of TOTAL genetic variance (state this in the caption).

# ==============================================================================
# 3  WEIR & COCKERHAM FST OUTLIERS   (yam 06_Fst_outliers_revised.R)
# ==============================================================================
# Per-locus WC (1984) FST from 0/1/2 genotypes (NA = missing), any number of groups.
wc_fst <- function(X, g) {
  g <- factor(g); r <- nlevels(g)
  ni <- sapply(levels(g), function(k) colSums(!is.na(X[g == k, , drop = FALSE])))
  pi <- sapply(levels(g), function(k) colSums(X[g == k, , drop = FALSE], na.rm = TRUE)) / (2 * ni)
  hi <- sapply(levels(g), function(k) colSums(X[g == k, , drop = FALSE] == 1, na.rm = TRUE)) / ni
  ni <- matrix(ni, ncol = r); pi <- matrix(pi, ncol = r); hi <- matrix(hi, ncol = r)
  nbar <- rowSums(ni) / r; nc <- (r * nbar - rowSums(ni^2) / (r * nbar)) / (r - 1)
  pbar <- rowSums(ni * pi, na.rm = TRUE) / (r * nbar)
  s2 <- rowSums(ni * (pi - pbar)^2, na.rm = TRUE) / ((r - 1) * nbar)
  hbar <- rowSums(ni * hi, na.rm = TRUE) / (r * nbar)
  a <- nbar / nc * (s2 - 1 / (nbar - 1) * (pbar * (1 - pbar) - (r - 1) / r * s2 - hbar / 4))
  b <- nbar / (nbar - 1) * (pbar * (1 - pbar) - (r - 1) / r * s2 - (2 * nbar - 1) / (4 * nbar) * hbar)
  cc <- hbar / 2
  ok <- is.finite(a + b + cc) & (a + b + cc) > 0
  fst <- ifelse(ok, a / (a + b + cc), NA)
  list(fst = fst, global = sum(a[ok]) / sum((a + b + cc)[ok]))
}
msg("== 3. WC FST ==")
grp <- meta$Population
glob <- wc_fst(G, grp)
fst <- data.frame(snp[, c("SNP", "CHROM", "POS")], FST_global = glob$fst)
msg("Genome-wide WC FST across", length(unique(grp)), "groups =", round(glob$global, 3))

# private alleles: alt allele present in exactly one group (strict set)
present <- sapply(split(seq_along(grp), grp), function(i) colSums(G[i, , drop = FALSE], na.rm = TRUE) > 0)
fst$private_in <- ifelse(rowSums(present) == 1, colnames(present)[max.col(present)], NA)

# non-admixed samples: Q-matrix supplied, else sNMF with K = number of groups
if (!is.na(CFG$qmatrix_csv)) {
  qd <- read.csv(CFG$qmatrix_csv, check.names = FALSE); Q <- as.matrix(qd[match(ids, trimws(qd[[1]])), -1])
} else {
  msg("Ancestry for non-admixed samples: LEA::snmf, K =", length(unique(grp)))
  gf <- res("geno_full_for_snmf.geno"); Gi <- maf_filter(G); LEA::write.geno(Gi, gf)
  invisible(capture.output(                        # sNMF prints a lot; keep the console readable
    sn <- LEA::snmf(gf, K = length(unique(grp)), repetitions = 3, entropy = TRUE, project = "new", seed = CFG$seed)))
  best <- which.min(LEA::cross.entropy(sn, K = length(unique(grp))))
  Q <- LEA::Q(sn, K = length(unique(grp)), run = best)
}
main_cl <- max.col(Q); maxq <- apply(Q, 1, max)
grp_cl <- tapply(main_cl, grp, function(x) as.integer(names(which.max(table(x)))))
pure <- !is.na(maxq) & maxq >= CFG$nonadmixed_q & main_cl == grp_cl[grp]
write.csv(data.frame(SampleID = ids, Population = grp, max_ancestry = round(maxq, 3), main_cluster = main_cl, non_admixed = pure),
          tab("non_admixed_samples.csv"), row.names = FALSE)
msg("Non-admixed (max ancestry >=", CFG$nonadmixed_q, "):", paste(names(table(grp[pure])), table(grp[pure]), sep = "=", collapse = ", "))

# pairwise contrasts among non-admixed groups
pg <- sort(unique(grp[pure])); pg <- pg[sapply(pg, function(k) sum(grp == k & pure)) >= 5]
pairs <- combn(pg, 2, simplify = FALSE)
pw <- sapply(pairs, function(pp) { i <- grp %in% pp & pure; wc_fst(G[i, ], grp[i])$fst })
colnames(pw) <- sapply(pairs, function(pp) paste(pp, collapse = " vs "))
fst <- cbind(fst, pw)
pw_global <- sapply(pairs, function(pp) { i <- grp %in% pp & pure; wc_fst(G[i, ], grp[i])$global })
write.csv(data.frame(contrast = colnames(pw), n_samples = sapply(pairs, function(pp) sum(grp %in% pp & pure)), WC_FST = pw_global),
          tab("FST_pairwise_between_groups_genomewide.csv"), row.names = FALSE)
print(data.frame(contrast = colnames(pw), WC_FST = round(pw_global, 3)))

q <- quantile(fst$FST_global, c(.95, .99, .995), na.rm = TRUE)
fst$outlier <- cut(fst$FST_global, c(-Inf, q, Inf), labels = c("none", "95th", "99th", "99.5th"))
write.csv(fst, res("fst_per_locus_WC.csv"), row.names = FALSE)
write.csv(fst %>% filter(outlier != "none"), tab("FST_outlier_loci_global.csv"), row.names = FALSE)
write.csv(fst %>% filter(outlier != "none") %>% count(CHROM, outlier), tab("FST_outlier_counts_by_chromosome.csv"), row.names = FALSE)
msg("Outliers with private alleles:", sum(fst$outlier != "none" & !is.na(fst$private_in), na.rm = TRUE),
    "of", sum(fst$outlier != "none", na.rm = TRUE))

# Figure: per-locus global WC FST along the chromosomes, contrast stated in the title
plot_fst_by_chromosome <- function(d) {
  d <- d %>% filter(!is.na(FST_global)) %>% mutate(CHROM = factor(CHROM, levels = chr_levels)) %>% arrange(CHROM, POS) %>%
    mutate(idx = row_number(), shade = factor(as.integer(CHROM) %% 2))
  mids <- d %>% group_by(CHROM) %>% summarise(mid = median(idx), .groups = "drop")
  ggplot(d, aes(idx, FST_global)) +
    geom_point(aes(colour = shade), size = 0.6, show.legend = FALSE) +
    geom_point(data = filter(d, outlier != "none"), aes(fill = outlier), shape = 21, size = 1.6, colour = "black", stroke = 0.2) +
    geom_hline(yintercept = q, linetype = c(3, 2, 1), colour = "grey30", linewidth = 0.3) +
    scale_colour_manual(values = c("0" = "grey35", "1" = "grey65")) +
    scale_fill_manual(values = c("95th" = "#FDB863", "99th" = "#E66101", "99.5th" = "#B2182B"), name = "Outlier\n(percentile)") +
    scale_x_continuous(breaks = mids$mid, labels = mids$CHROM, expand = c(0.005, 0)) +
    labs(x = "Chromosome (A and B sub-genomes; Un = unplaced scaffolds)", y = expression("Weir & Cockerham " * F[ST]),
         title = paste0("Per-locus Weir & Cockerham FST, global across ", length(unique(grp)), " groups (",
                        paste(sort(unique(grp)), collapse = ", "), "); genome-wide FST = ", round(glob$global, 3)),
         subtitle = "Lines: 95th (dotted), 99th (dashed) and 99.5th (solid) percentiles") +
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
}
save_fig("FST_per_locus_global_by_chromosome", plot_fst_by_chromosome(fst), 13, 4.5)

# ==============================================================================
# 4  LFMM2   (yam 09_LFMM.R)
# ==============================================================================
msg("== 4. LFMM2 ==")
hits <- if (file.exists(CFG$egwas_csv)) read.csv(CFG$egwas_csv, check.names = FALSE) else NULL
if (!is.null(hits)) {                        # environmental_gwas_complete.R names -> SNP / trait
  if ("marker" %in% names(hits)) hits$SNP <- hits$marker
  if ("variable" %in% names(hits)) hits$trait <- hits$variable
  if (!"set" %in% names(hits)) hits$set <- "full"
  msg("eGWAS hits for overlap:", nrow(hits))
} else msg("No eGWAS file (", CFG$egwas_csv, "): shared_with_EGWAS set to FALSE")
lf_out <- list()
for (s in names(S)) {
  set_ids <- S[[s]]; Gs <- maf_filter(G[set_ids, , drop = FALSE]); E <- scale(E_all[set_ids, vif_keep, drop = FALSE])
  f <- res(paste0("geno_", s, ".lfmm")); LEA::write.lfmm(Gs, f)
  K <- if (s == "full") CFG$K_lfmm_full else CFG$K_lfmm_within[[s]]
  if (is.null(K)) K <- 2
  # --- optional check of K (uncomment):
  # sn <- LEA::snmf(LEA::lfmm2geno(f), K = 1:8, entropy = TRUE, repetitions = 3, project = "new"); plot(sn)
  mod <- LEA::lfmm2(input = Gs, env = E, K = K, effect.sizes = TRUE)
  pv <- LEA::lfmm2.test(mod, input = Gs, env = E, full = FALSE, genomic.control = TRUE)
  pm <- pv$pvalues; if (is.null(dim(pm))) pm <- matrix(pm, ncol = 1)
  if (nrow(pm) != ncol(Gs)) pm <- base::t(pm)                     # orientation differs among LEA versions
  colnames(pm) <- colnames(E)
  qv <- apply(pm, 2, p.adjust, method = "BH")
  sig <- which(qv < CFG$fdr, arr.ind = TRUE)
  si <- unname(sig[, 1]); sj <- unname(sig[, 2]); sn <- colnames(Gs)[si]
  gif <- if (length(pv$gif) == ncol(E)) unname(pv$gif) else rep(NA_real_, ncol(E))
  tabs <- data.frame(set = rep(s, length(si)), K = rep(K, length(si)), SNP = sn,
                     CHROM = snp$CHROM[match(sn, snp$SNP)], POS = snp$POS[match(sn, snp$SNP)],
                     trait = colnames(E)[sj], p = pm[cbind(si, sj)], q = qv[cbind(si, sj)],
                     GIF = round(gif[sj], 3), stringsAsFactors = FALSE)
  tabs$shared_with_EGWAS <- if (is.null(hits)) FALSE else
    paste(tabs$SNP, tabs$trait) %in% paste(hits$SNP[hits$set == s], hits$trait[hits$set == s])
  lf_out[[s]] <- tabs
  saveRDS(list(mod = mod, snps = colnames(Gs), env = colnames(E), pvalues = pm, gif = pv$gif), res(paste0("lfmm_", s, ".rds")))
  msg("  ", s, ": K =", K, "|", ncol(Gs), "SNPs x", ncol(E), "variables |", nrow(tabs), "SNP-variable associations at FDR", CFG$fdr)
}
s10 <- bind_rows(lf_out); write.csv(s10, tab("LFMM2_significant_SNP_environment_associations.csv"), row.names = FALSE)
print(s10 %>% group_by(set) %>% summarise(n_SNP_trait = n(), n_SNPs = n_distinct(SNP), shared = sum(shared_with_EGWAS)))

writeLines(capture.output(sessionInfo()), res("sessionInfo.txt"))
msg("Done ->", normalizePath(CFG$out_dir))