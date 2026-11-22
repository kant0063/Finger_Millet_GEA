################################################################################
# finger_millet_genomic_prediction_V3.R
#
# V3 changes (from V2):
#  - Figures are written as vector PDF (CFG$fig_format; "png" still possible).
#  - Country and region are added to the GEAV heatmaps as annotation bars, and a second
#    heatmap clusters all samples together so GEAV clusters can be compared with
#    population, region and country. Cluster composition by country and region is
#    written to tables. Country/region come from the coordinate file; if there is no
#    country column, country is looked up from the coordinates (needs terra + geodata).
#  - GenoBridge section removed.
#
# Genomic prediction of environment of origin (environmental genomic selection)
# for finger millet (Eleusine coracana), at whole-genome and sub-genome (A, B) level.
#
# V2 changes, addressing the problems raised in review of the yam (Dioscorea) version:
#  1  Population structure. Random k-fold accuracy is inflated by structure, so every
#     accuracy is reported with (a) a structure-only baseline whose PCs are now computed
#     from the training samples only and projected onto the test samples, (b) the gain
#     over that baseline with a bootstrap CI, and (c) within-group accuracy (structure
#     removed). Leave-one-population-out and spatial blocks are the primary schemes.
#  2  Spatial autocorrelation. Spatial-block CV can add a buffer: training samples within
#     CFG$spatial_buffer_km of any test sample are dropped ("spatial blocks + buffer").
#  3  Uncertainty. 95% bootstrap CIs (individuals resampled within groups) for pooled r,
#     within-group r and gain over baseline, for every trait x model x SNP set x scheme.
#     Spearman rho and RMSE are reported alongside Pearson r.
#  4  Like-for-like model comparison. Slow models were run on a subset of traits, so
#     their means were not comparable with means over all traits. Summaries now give
#     means over all traits AND over the common trait set on which every model was run.
#  5  GEAV model choice and parent selection. V1 picked the model by spatial-block
#     accuracy, which was ~0 and no better than the baseline. V2 picks it by within-group
#     accuracy (the accuracy relevant to within-group GEAVs) and only ranks candidate
#     parents for variables whose within-group accuracy CI excludes zero; others are
#     reported as not predictable.
#  6  Core collections. The core used as a training set is chosen a priori by
#     CFG$core_rule, not by its own test accuracy (which is optimistic). A group-balanced
#     maximum-dissimilarity core is added, and cores missing a group are flagged (the yam
#     core lacking Cameroonian D. praehensilis performed much worse).
#  7  rrBLUP is fitted in its equivalent GBLUP form (additive VanRaden kernel), which is
#     much faster and gives identical predictions; test-set BLUPs come from the training
#     fit only.
#
# Whole-genome and sub-genome assessment:
#  - SNP sets: whole genome, each sub-genome, and the larger sub-genome subsampled to the
#    SNP number of the smaller (so A-vs-B differences are not just SNP number).
#  - Prediction accuracy per SNP set (all schemes), plus a multi-kernel A + B model.
#  - Genomic variance partitioned between sub-genomes (two-kernel Bayesian GBLUP).
#  - Agreement between A and B relationship matrices, and sub-genome-specific GEAVs.
#  - Optional removal of SNPs with excess heterozygosity (likely collapsed homeologs;
#    het > 2pq(1-F), F = 0.5, as in Bancic et al. 2024, Plant Genome 17:e20392).
#
# Inputs: a VCF and a coordinate CSV (SampleID, Latitude, Longitude, [Population],
# [Country], [Region]).
# Finger millet: the TASSEL VCF numbers chromosomes 1-19; 1-18 are mapped to 1A, 1B,
# 2A, 2B ... 9A, 9B (checked against Bancic et al. 2024 Table 1); 19 = unplaced
# scaffolds, used in the whole-genome set only. For a diploid set subgenome_regex = NULL.
# Run with:  Rscript finger_millet_genomic_prediction_V3.R
################################################################################

# ==============================================================================
# 0  CONFIGURATION  (edit this block)
# ==============================================================================
CFG <- list(
  species        = "Finger millet (Eleusine coracana)",
  
  # ---- inputs ----
  vcf            = "finger_millet_404.vcf",       # biallelic SNPs; sample IDs = SampleID
  coords_csv     = "sample_coordinates.csv",      # SampleID, Latitude, Longitude, [Population], [Country], [Region]
  env_csv        = "millet_extracted_bioclim_values.csv",   # optional: SampleID + environmental variables; if NA, WorldClim
  # bioclim + elevation are downloaded and extracted at the coordinates
  extra_env_csv  = NA,   # optional: SampleID + extra variables, merged in
  core_csv       = "millet_core_n100.csv",   # optional: an existing core collection (first column = SampleID)
  
  # ---- climate download (only when env_csv is NA) ----
  clim_cache     = "climate_cache",
  worldclim_res  = 2.5,
  add_elevation  = TRUE,
  add_soil_ph    = FALSE,
  
  # ---- traits ----
  env_vars       = NULL,                          # NULL = all; or e.g. c("bio1","bio5","bio12","bio14","elev")
  n_env_pcs      = 3,
  
  # ---- genotype filters ----
  maf            = 0.05,
  max_snp_missing = 0.20,
  max_ind_missing = 0.50,
  homeolog_het_F = 0.5,  # drop SNPs with observed het > 2pq(1-F); NA = keep all
  
  # ---- sub-genomes ----
  chr_map        = setNames(paste0(rep(1:9, each = 2), c("A", "B")), 1:18),  # VCF CHROM -> reference
  # chromosome; NULL if the VCF already uses 1A, 1B...
  subgenome_regex = "([AB])$",                    # NULL for diploids (sub-genome analyses skipped)
  matched_snp_control = TRUE,                     # add larger sub-genome subsampled to the smaller's SNP number
  
  # ---- populations / groups ----
  n_groups_if_no_population = 4,                  # k-means on genotype PCs if no Population column
  country_col    = NULL,   # NULL = auto-detect a column named Country in coords_csv
  region_col     = NULL,   # NULL = auto-detect Region / Area; if none, Population is used
  derive_country = TRUE,   # no country column: look up country from coordinates (terra + geodata)
  
  # ---- cross-validation ----
  n_pc           = 5,                             # PCs in the structure-only baseline
  cv_k           = 10,
  cv_reps        = 5,                             # random k-fold repeats (light models)
  heavy_reps     = 2,                             # repeats for slow models
  n_spatial_blocks = 5,
  spatial_buffer_km = 100,                        # 0 = no buffered scheme
  light_models   = c("GBLUP_add", "GBLUP_gauss", "GBLUP_exp"),   # GBLUP_add = rrBLUP
  heavy_models   = c("BayesB", "RKHS", "RF", "MK_AB"),           # MK_AB = multi-kernel A + B
  heavy_traits   = c("bio1", "bio5", "bio12", "bio14", "elev", "envPC1"),  # "all" = every trait
  bglr_iter      = 6000,
  bglr_burn      = 1000,
  n_boot         = 500,                           # bootstrap replicates for CIs
  conf           = 0.95,
  
  # ---- sub-genome variance partition ----
  vc_iter        = 12000,
  vc_burn        = 2000,
  
  # ---- core collection ----
  core_size      = NA,                            # NA = 30% of accessions
  core_reps      = 50,
  core_rule      = "balanced_max_dissimilarity",  # core used as training set, chosen a priori:
  # balanced_max_dissimilarity, max_dissimilarity,
  # corehunter_genotype or supplied
  
  # ---- candidate parents (within-population GEAVs) ----
  targets = list(
    drought    = list(var = "bio12", direction = "low",  label = "Low annual precipitation (BIO12)"),
    dry_season = list(var = "bio14", direction = "low",  label = "Low precipitation of driest month (BIO14)"),
    heat       = list(var = "bio5",  direction = "high", label = "High max temperature of warmest month (BIO5)"),
    acid_soil  = list(var = "phh2o", direction = "low",  label = "Low soil pH (skipped if not available)")
  ),
  n_candidates   = 5,
  n_subclusters  = 3,
  
  # ---- output / misc ----
  out_dir        = "genomic_prediction_results",
  fig_format     = "pdf",                          # "pdf" (vector) or "png"
  cores          = 4,                             # >1 = parallel over traits (Linux/macOS)
  seed           = 2026,
  run = c(cv = TRUE, subgenome_vc = TRUE, core = TRUE, geav = TRUE)
)

# ==============================================================================
# Packages and helpers
# ==============================================================================
need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
}
need(c("rrBLUP", "BGLR", "ranger", "dplyr", "tidyr", "ggplot2", "patchwork"))
if (is.na(CFG$env_csv)) need(c("terra", "geodata"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

set.seed(CFG$seed)
dir.create(file.path(CFG$out_dir, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CFG$out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
tab <- function(f) file.path(CFG$out_dir, "tables", f)
fig <- function(f) file.path(CFG$out_dir, "figures", f)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
save_fig <- function(name, plot, width, height) {
  f <- fig(paste0(name, ".", CFG$fig_format))
  if (CFG$fig_format == "pdf") {
    dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggsave(f, plot, width = width, height = height, device = dev, bg = "white", limitsize = FALSE)
  } else ggsave(f, plot, width = width, height = height, dpi = 600, bg = "white", limitsize = FALSE)
}
papply <- function(X, FUN) if (CFG$cores > 1 && .Platform$OS.type == "unix")
  parallel::mclapply(X, FUN, mc.cores = CFG$cores, mc.preschedule = FALSE) else lapply(X, FUN)

# VCF -> individuals x SNPs ALT-allele dosage (NA = missing), with chromosome and position
read_vcf_dosage <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
  n_skip <- 0
  repeat {
    l <- readLines(con, n = 1)
    if (!length(l)) { close(con); stop("No #CHROM header in ", path) }
    if (startsWith(l, "#CHROM")) break
    n_skip <- n_skip + 1
  }
  close(con)
  samples <- strsplit(l, "\t")[[1]][-(1:9)]
  v <- read.delim(path, skip = n_skip + 1, header = FALSE, colClasses = "character",
                  comment.char = "", quote = "")
  bi <- !grepl(",", v[[5]]) & nchar(v[[4]]) == 1 & nchar(v[[5]]) == 1
  v <- v[bi, , drop = FALSE]
  gt <- sub(":.*$", "", as.matrix(v[, -(1:9), drop = FALSE]))
  a1 <- sub("[/|].*$", "", gt); a2 <- sub("^[^/|]*[/|]", "", gt)
  d <- (a1 != "0") + (a2 != "0"); d[a1 == "." | a2 == "."] <- NA
  storage.mode(d) <- "numeric"
  G <- t(d); rownames(G) <- trimws(samples); colnames(G) <- paste(v[[1]], v[[2]], sep = "_")
  list(G = G, chrom = v[[1]], pos = as.integer(v[[2]]))
}
haversine_km <- function(lat, lon) {
  r <- pi / 180; la <- lat * r; lo <- lon * r
  dlat <- outer(la, la, "-"); dlon <- outer(lo, lo, "-")
  a <- sin(dlat / 2)^2 + outer(cos(la), cos(la)) * sin(dlon / 2)^2
  2 * 6371 * asin(pmin(sqrt(a), 1))
}
r_within <- function(p, y, g) {
  v <- vapply(split(seq_along(y), g), function(i)
    if (length(i) > 4 && sd(p[i]) > 0 && sd(y[i]) > 0) cor(p[i], y[i]) else NA_real_, 0)
  if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
}
safe_cor <- function(a, b, ...) if (sd(a) > 0 && sd(b) > 0) cor(a, b, ...) else NA_real_

# ==============================================================================
# 1  DATA
# ==============================================================================
msg("== 1. Data ==", CFG$species)
co <- read.csv(CFG$coords_csv, stringsAsFactors = FALSE, check.names = FALSE)
nm <- tolower(names(co))
lat_col <- names(co)[grep("^lat", nm)[1]]; lon_col <- names(co)[grep("^lon", nm)[1]]
pop_col <- names(co)[grep("^pop", nm)[1]]
if (is.na(lat_col) || is.na(lon_col)) stop("Coordinate file needs latitude and longitude columns")
ctry_col <- if (!is.null(CFG$country_col)) CFG$country_col else names(co)[grep("^country", nm)[1]]
reg_col  <- if (!is.null(CFG$region_col)) CFG$region_col else names(co)[grep("^region|^area|^geo.*area", nm)[1]]
meta <- data.frame(SampleID = trimws(co[[1]]), Latitude = co[[lat_col]], Longitude = co[[lon_col]],
                   Population = if (!is.na(pop_col)) as.character(co[[pop_col]]) else NA,
                   Country = if (!is.na(ctry_col)) as.character(co[[ctry_col]]) else NA,
                   Region = if (!is.na(reg_col)) as.character(co[[reg_col]]) else NA,
                   stringsAsFactors = FALSE)
meta <- meta[!duplicated(meta$SampleID), ]
if (all(is.na(meta$Country)) && CFG$derive_country) {
  meta$Country <- tryCatch({                       # country polygon containing each point
    need(c("terra", "geodata"))
    w <- geodata::world(resolution = 5, level = 0, path = CFG$clim_cache)
    pts0 <- terra::vect(meta, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
    as.character(terra::extract(w, pts0)$NAME_0)
  }, error = function(e) { msg("Country lookup from coordinates failed:", conditionMessage(e)); NA })
  if (!all(is.na(meta$Country))) msg("Country looked up from coordinates for", sum(!is.na(meta$Country)), "samples")
}
meta$Country[is.na(meta$Country) | meta$Country == ""] <- "Unknown"

vc <- read_vcf_dosage(CFG$vcf); G_raw <- vc$G
msg("VCF:", nrow(G_raw), "samples x", ncol(G_raw), "biallelic SNPs")

if (!is.na(CFG$env_csv)) {
  env <- read.csv(CFG$env_csv, check.names = FALSE, stringsAsFactors = FALSE)
  names(env)[1] <- "SampleID"; env$SampleID <- trimws(env$SampleID)
} else {
  pts <- terra::vect(meta, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
  wc <- geodata::worldclim_global("bio", res = CFG$worldclim_res, path = CFG$clim_cache)
  names(wc) <- paste0("bio", sub("^.*?([0-9]+)$", "\\1", names(wc), perl = TRUE))
  env <- data.frame(SampleID = meta$SampleID, terra::extract(wc, pts)[, names(wc)])
  if (CFG$add_elevation) {
    el <- geodata::elevation_global(res = CFG$worldclim_res, path = CFG$clim_cache)
    env$elev <- terra::extract(el, pts)[, 2]
  }
  if (CFG$add_soil_ph) {
    sp <- tryCatch(geodata::soil_world(var = "phh2o", depth = 5, path = CFG$clim_cache),
                   error = function(e) { msg("Soil pH download failed:", conditionMessage(e)); NULL })
    if (!is.null(sp)) env$phh2o <- terra::extract(sp, pts)[, 2]
  }
  write.csv(env, tab("environment_of_origin.csv"), row.names = FALSE)
}
if (!is.na(CFG$extra_env_csv)) {
  ex <- read.csv(CFG$extra_env_csv, check.names = FALSE, stringsAsFactors = FALSE)
  names(ex)[1] <- "SampleID"; ex$SampleID <- trimws(ex$SampleID)
  env <- left_join(env, ex, by = "SampleID")
}
if (!is.null(CFG$env_vars)) env <- env[, c("SampleID", intersect(CFG$env_vars, names(env)))]

env_ok <- env$SampleID[complete.cases(env)]
if (length(setdiff(env$SampleID, env_ok)))
  msg("Dropped", length(setdiff(env$SampleID, env_ok)), "sample(s) with missing environmental values")
ids <- Reduce(intersect, list(rownames(G_raw), meta$SampleID, env_ok))
msg("Samples in VCF / coordinates / environment / all three:", nrow(G_raw), "/", nrow(meta), "/", nrow(env), "/", length(ids))
if (length(ids) < 30) stop("Fewer than 30 samples shared by the VCF, coordinates and environment")
G <- G_raw[ids, , drop = FALSE]

# ---- genotype QC ----
ind_miss <- rowMeans(is.na(G))
if (any(ind_miss > CFG$max_ind_missing)) {
  msg("Dropped", sum(ind_miss > CFG$max_ind_missing), "sample(s) with >", CFG$max_ind_missing * 100, "% missing")
  ids <- ids[ind_miss <= CFG$max_ind_missing]; G <- G[ids, , drop = FALSE]
}
meta <- meta[match(ids, meta$SampleID), ]
snp_miss <- colMeans(is.na(G)); p <- colMeans(G, na.rm = TRUE) / 2
het <- colMeans(G == 1, na.rm = TRUE)
keep <- !is.na(p) & snp_miss <= CFG$max_snp_missing & pmin(p, 1 - p) >= CFG$maf
n_het_drop <- 0
if (!is.na(CFG$homeolog_het_F)) {
  het_bad <- het > 2 * p * (1 - p) * (1 - CFG$homeolog_het_F)
  n_het_drop <- sum(keep & het_bad, na.rm = TRUE); keep <- keep & !het_bad
}
G <- G[, keep, drop = FALSE]; chrom <- vc$chrom[keep]; p <- p[keep]
G <- apply(G, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x }); rownames(G) <- ids
msg("SNPs after missingness and MAF filters:", ncol(G), "| removed for excess heterozygosity (likely homeologs):", n_het_drop)

# ---- SNP sets: whole genome and sub-genomes ----
ref_chr <- if (!is.null(CFG$chr_map)) ifelse(chrom %in% names(CFG$chr_map), CFG$chr_map[chrom], chrom) else chrom
subg <- rep(NA_character_, length(ref_chr))
if (!is.null(CFG$subgenome_regex)) {
  m <- regmatches(ref_chr, regexec(CFG$subgenome_regex, ref_chr))
  subg <- vapply(m, function(z) if (length(z) >= 2) z[2] else NA_character_, "")
}
snp_sets <- list(whole = seq_len(ncol(G)))
sg_names <- sort(unique(na.omit(subg)))
DO_SUB <- length(sg_names) >= 2
if (DO_SUB) {
  for (s in sg_names) snp_sets[[s]] <- which(subg == s)
  if (CFG$matched_snp_control) {
    n_min <- min(lengths(snp_sets[sg_names])); small <- sg_names[which.min(lengths(snp_sets[sg_names]))]
    set.seed(CFG$seed)
    for (s in setdiff(sg_names, small)) if (length(snp_sets[[s]]) > n_min)
      snp_sets[[paste0(s, "_matched")]] <- sort(sample(snp_sets[[s]], n_min))
  }
} else if (!is.null(CFG$subgenome_regex)) msg("Fewer than two sub-genomes found: sub-genome analyses skipped")
snp_tab <- data.frame(snp_set = names(snp_sets), n_snps = lengths(snp_sets))
snp_tab$chromosomes <- vapply(snp_sets, function(i) paste(unique(ref_chr[i]), collapse = ","), "")
write.csv(snp_tab, tab("snp_sets.csv"), row.names = FALSE)
msg("SNP sets:", paste(snp_tab$snp_set, snp_tab$n_snps, sep = "=", collapse = ", "),
    "| unassigned (scaffolds):", sum(is.na(subg)))

# ---- kernels per SNP set ----
Gc <- scale(G, scale = FALSE)
make_kernels <- function(idx) {
  X <- Gc[, idx, drop = FALSE]
  K_add <- tcrossprod(X) / (2 * sum(p[idx] * (1 - p[idx])))           # VanRaden GRM
  D <- as.matrix(dist(X)); D <- D / max(D); th <- median(D[upper.tri(D)])
  list(add = K_add, gauss = exp(-(D / th)^2), exp = exp(-D / th))
}
KS <- lapply(snp_sets, make_kernels)

# ---- populations ----
eg <- eigen(KS$whole$add, symmetric = TRUE)
PCg <- eg$vectors[, 1:min(CFG$n_pc, nrow(G) - 1), drop = FALSE] %*% diag(sqrt(pmax(eg$values[1:min(CFG$n_pc, nrow(G) - 1)], 0)))
if (all(is.na(meta$Population))) {
  meta$Population <- paste0("G", kmeans(PCg, CFG$n_groups_if_no_population, nstart = 50)$cluster)
  msg("No Population column: groups from k-means on genotype PCs")
}
grp <- meta$Population
if (all(is.na(meta$Region))) meta$Region <- meta$Population
meta$Region[is.na(meta$Region) | meta$Region == ""] <- "Unknown"
msg("Populations:", paste(names(table(grp)), table(grp), sep = "=", collapse = ", "))
msg("Regions:", paste(names(table(meta$Region)), table(meta$Region), sep = "=", collapse = ", "),
    "| countries:", length(unique(meta$Country)))

# ---- traits ----
E <- as.matrix(env[match(ids, env$SampleID), setdiff(names(env), "SampleID"), drop = FALSE])
E <- E[, apply(E, 2, sd) > 0, drop = FALSE]
Y <- scale(E)
if (CFG$n_env_pcs > 0) {
  epc <- prcomp(Y)$x[, 1:min(CFG$n_env_pcs, ncol(Y)), drop = FALSE]
  colnames(epc) <- paste0("envPC", seq_len(ncol(epc))); Y <- cbind(Y, scale(epc))
}
rownames(Y) <- ids
heavy <- if (identical(CFG$heavy_traits, "all")) colnames(Y) else intersect(CFG$heavy_traits, colnames(Y))
msg("Traits:", ncol(Y), "| slow models run on:", paste(heavy, collapse = ", "))

# ---- folds (shared by all models so comparisons are paired) ----
set.seed(CFG$seed)
folds_random <- replicate(CFG$cv_reps, split(seq_along(ids), sample(rep(seq_len(CFG$cv_k), length.out = length(ids)))),
                          simplify = FALSE)
folds_lopo <- split(seq_along(ids), grp)
n_blocks <- min(CFG$n_spatial_blocks, nrow(unique(meta[, c("Longitude", "Latitude")])))
km_sp <- kmeans(meta[, c("Longitude", "Latitude")], n_blocks, nstart = 50)
folds_space <- split(seq_along(ids), km_sp$cluster)
meta$SpatialBlock <- km_sp$cluster
DIST_KM <- haversine_km(meta$Latitude, meta$Longitude)
schemes <- list("leave-one-population-out" = list(folds = list(folds_lopo), buffer = 0),
                "spatial blocks" = list(folds = list(folds_space), buffer = 0))
if (CFG$spatial_buffer_km > 0)
  schemes[[sprintf("spatial blocks + %g km buffer", CFG$spatial_buffer_km)]] <- list(folds = list(folds_space), buffer = CFG$spatial_buffer_km)
schemes <- c(setNames(list(list(folds = folds_random, buffer = 0)), paste0("random ", CFG$cv_k, "-fold")), schemes)
RANDOM <- names(schemes)[1]
train_idx <- function(te, buffer) {
  tr <- setdiff(seq_along(ids), te)
  if (buffer > 0) tr <- tr[apply(DIST_KM[tr, te, drop = FALSE], 1, min) > buffer]
  tr
}
if (CFG$spatial_buffer_km > 0) {
  nb <- sapply(folds_space, function(te) length(train_idx(te, CFG$spatial_buffer_km)))
  msg("Training samples per spatial block with buffer:", paste(nb, collapse = ", "),
      "(without:", paste(sapply(folds_space, function(te) length(ids) - length(te)), collapse = ", "), ")")
}
write.csv(meta, tab("samples_used.csv"), row.names = FALSE)

# ---- models ----
gblup_pred <- function(K, y, tr, te) {
  Ktr <- K[tr, tr]; f <- rrBLUP::mixed.solve(y[tr], K = Ktr)
  ridge <- 1e-6 * mean(diag(Ktr))
  as.numeric(f$beta)[1] + as.numeric(K[te, tr, drop = FALSE] %*% solve(Ktr + diag(ridge, length(tr)), f$u))
}
baseline_pred <- function(y, tr, te) {          # PCs from training samples only, test projected
  K <- KS$whole$add; e <- eigen(K[tr, tr], symmetric = TRUE)
  k <- min(CFG$n_pc, length(tr) - 2); keep <- e$values[1:k] > 1e-8; k <- sum(keep)
  U <- e$vectors[, 1:k, drop = FALSE]; lam <- e$values[1:k]
  Xtr <- U %*% diag(sqrt(lam), k); Xte <- K[te, tr, drop = FALSE] %*% U %*% diag(1 / sqrt(lam), k)
  f <- lm.fit(cbind(1, Xtr), y[tr]); as.numeric(cbind(1, Xte) %*% f$coefficients)
}
bglr_pred <- function(ETA, y, tr, te) {
  yNA <- rep(NA_real_, length(y)); yNA[tr] <- y[tr]
  f <- BGLR::BGLR(yNA, ETA = ETA, nIter = CFG$bglr_iter, burnIn = CFG$bglr_burn, verbose = FALSE,
                  saveAt = file.path(tempdir(), paste0("bglr_", Sys.getpid(), "_")))
  f$yHat[te]
}
fit_pred <- function(model, set, y, tr, te) {
  K <- KS[[set]]
  pred <- switch(model,
                 baseline    = baseline_pred(y, tr, te),
                 GBLUP_add   = gblup_pred(K$add, y, tr, te),
                 GBLUP_gauss = gblup_pred(K$gauss, y, tr, te),
                 GBLUP_exp   = gblup_pred(K$exp, y, tr, te),
                 BayesB      = bglr_pred(list(list(X = Gc[, snp_sets[[set]], drop = FALSE], model = "BayesB")), y, tr, te),
                 RKHS        = bglr_pred(list(list(K = K$gauss, model = "RKHS")), y, tr, te),
                 RF          = { f <- ranger::ranger(x = G[tr, snp_sets[[set]], drop = FALSE], y = y[tr], num.trees = 500,
                                                     seed = CFG$seed, num.threads = 1)
                 predict(f, G[te, snp_sets[[set]], drop = FALSE])$predictions },
                 MK_AB       = bglr_pred(lapply(sg_names, function(s) list(K = KS[[s]]$add, model = "RKHS")), y, tr, te),
                 stop("unknown model: ", model))
  as.numeric(pred)
}
run_cv <- function(model, set, y, folds, buffer) {
  pred <- rep(NA_real_, length(y))
  for (te in folds) pred[te] <- fit_pred(model, set, y, train_idx(te, buffer), te)
  pred
}
# Bootstrap: individuals resampled within groups; paired with the baseline
boot_stats <- function(pred, base, y, g) {
  ok <- !is.na(pred); pred <- pred[ok]; base <- base[ok]; y <- y[ok]; g <- g[ok]
  gi <- split(seq_along(y), g)
  B <- replicate(CFG$n_boot, {
    i <- unlist(lapply(gi, function(v) v[sample.int(length(v), length(v), replace = TRUE)]))
    c(safe_cor(pred[i], y[i]), r_within(pred[i], y[i], g[i]),
      safe_cor(pred[i], y[i]) - safe_cor(base[i], y[i]), r_within(pred[i], y[i], g[i]) - r_within(base[i], y[i], g[i]))
  })
  a <- (1 - CFG$conf) / 2
  q <- unname(apply(B, 1, quantile, probs = c(a, 1 - a), na.rm = TRUE))
  c(r_lo = q[1, 1], r_hi = q[2, 1], r_within_lo = q[1, 2], r_within_hi = q[2, 2],
    gain_lo = q[1, 3], gain_hi = q[2, 3], gain_within_lo = q[1, 4], gain_within_hi = q[2, 4])
}
point_stats <- function(pred, y, g) {
  ok <- !is.na(pred)
  c(r = safe_cor(pred[ok], y[ok]), r_within = r_within(pred[ok], y[ok], g[ok]),
    rho = safe_cor(pred[ok], y[ok], method = "spearman"), rmse = sqrt(mean((pred[ok] - y[ok])^2)))
}
combos_for <- function(trait) {
  cb <- data.frame(model = "baseline", set = "whole")
  cb <- rbind(cb, expand.grid(model = CFG$light_models, set = names(snp_sets), stringsAsFactors = FALSE))
  if (trait %in% heavy) {
    hm <- intersect(CFG$heavy_models, c("BayesB", "RKHS", "RF"))
    if (length(hm)) cb <- rbind(cb, data.frame(model = hm, set = "whole"))
    if (DO_SUB && "MK_AB" %in% CFG$heavy_models) cb <- rbind(cb, data.frame(model = "MK_AB", set = "A+B"))
  }
  cb
}
reps_for <- function(model) if (model %in% c("BayesB", "RKHS", "RF", "MK_AB")) CFG$heavy_reps else CFG$cv_reps
set_key <- function(set) if (set == "A+B") "whole" else set

# ==============================================================================
# 2  CROSS-VALIDATED PREDICTION ACCURACY (whole genome and sub-genomes)
# ==============================================================================
acc <- NULL; OOF <- list()
if (CFG$run["cv"]) {
  msg("== 2. Cross-validation ==")
  one_trait <- function(trait) {
    set.seed(CFG$seed)
    y <- Y[, trait]; cb <- combos_for(trait); rows <- list(); oof <- list()
    for (sch in names(schemes)) {
      S <- schemes[[sch]]; nrep <- if (sch == RANDOM) length(S$folds) else 1
      preds <- list()
      for (j in seq_len(nrow(cb))) {
        mo <- cb$model[j]; st <- cb$set[j]; nr <- min(nrep, reps_for(mo))
        pr <- sapply(seq_len(nr), function(r) run_cv(mo, set_key(st), y, S$folds[[r]], S$buffer))
        preds[[paste(mo, st)]] <- as.matrix(pr)
      }
      base1 <- preds[["baseline whole"]][, 1]
      for (j in seq_len(nrow(cb))) {
        key <- paste(cb$model[j], cb$set[j]); pr <- preds[[key]]
        ps <- apply(pr, 2, point_stats, y = y, g = grp)
        bs <- boot_stats(pr[, 1], base1, y, grp)
        rows[[paste(sch, key)]] <- data.frame(trait, model = cb$model[j], snp_set = cb$set[j], scheme = sch,
                                              n_reps = ncol(pr), r = mean(ps["r", ]), r_rep_sd = if (ncol(pr) > 1) sd(ps["r", ]) else NA,
                                              r_within = mean(ps["r_within", ], na.rm = TRUE), rho = mean(ps["rho", ]), rmse = mean(ps["rmse", ]),
                                              t(bs), check.names = FALSE)
        if (sch == RANDOM) oof[[key]] <- rowMeans(pr)
      }
    }
    msg("  ", trait, "done")
    list(acc = bind_rows(rows), oof = oof)
  }
  res <- papply(colnames(Y), one_trait); names(res) <- colnames(Y)
  acc <- bind_rows(lapply(res, `[[`, "acc"))
  OOF <- lapply(res, `[[`, "oof")
  base <- acc %>% filter(model == "baseline") %>% select(trait, scheme, r_base = r, r_within_base = r_within)
  acc <- acc %>% left_join(base, by = c("trait", "scheme")) %>%
    mutate(gain = r - r_base, gain_within = r_within - r_within_base,
           gain_significant = gain_lo > 0, r_within_significant = r_within_lo > 0)
  acc$scheme <- factor(acc$scheme, levels = names(schemes))
  write.csv(acc, tab("prediction_accuracy.csv"), row.names = FALSE)
  
  # Summaries: all traits and the common trait set (every model run) - like-for-like
  common <- acc %>% group_by(trait) %>% summarise(n_models = n_distinct(paste(model, snp_set)), .groups = "drop") %>%
    filter(n_models == max(n_models)) %>% pull(trait)
  summ <- function(d, lab) d %>% group_by(scheme, model, snp_set) %>%
    summarise(trait_set = lab, n_traits = n_distinct(trait), mean_r = mean(r, na.rm = TRUE),
              mean_r_within = mean(r_within, na.rm = TRUE), mean_gain = mean(gain, na.rm = TRUE),
              prop_traits_gain_significant = mean(gain_significant, na.rm = TRUE),
              prop_traits_r_within_significant = mean(r_within_significant, na.rm = TRUE), .groups = "drop")
  acc_sum <- bind_rows(summ(acc, "all traits"), summ(filter(acc, trait %in% common), "common trait set")) %>%
    arrange(trait_set, scheme, desc(mean_r))
  write.csv(acc_sum, tab("prediction_accuracy_summary.csv"), row.names = FALSE)
  msg("Common trait set (all models run):", paste(common, collapse = ", "))
  print(acc_sum %>% filter(trait_set == "all traits") %>% select(scheme, model, snp_set, n_traits, mean_r, mean_r_within, mean_gain) %>%
          mutate(across(where(is.numeric) & !n_traits, ~ round(.x, 3))), n = Inf)
  
  # Figure 1: whole genome, models by scheme
  d1 <- acc %>% filter(snp_set == "whole")
  p1 <- ggplot(d1, aes(model, r, fill = model)) + geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.size = .5) + facet_wrap(~scheme, nrow = 1) +
    labs(x = NULL, y = "Prediction accuracy (r)", title = paste0(CFG$species, ": environment of origin, whole genome"),
         subtitle = "One point per variable; baseline = genotype PCs only (population structure)") +
    theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig("accuracy_whole_genome_by_scheme", p1, width = 12, height = 4.8)
  
  # Figure 2: gain over baseline with CI, best GBLUP model per trait (whole genome)
  d2 <- acc %>% filter(snp_set == "whole", model %in% CFG$light_models) %>% group_by(trait, scheme) %>%
    slice_max(r, n = 1, with_ties = FALSE) %>% ungroup()
  p2 <- ggplot(d2, aes(gain, trait, colour = gain_significant)) + geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = gain_lo, xmax = gain_hi), height = 0) + geom_point() + facet_wrap(~scheme, nrow = 1) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey55"), name = "CI excludes 0") +
    labs(x = paste0("Gain in r over the structure-only baseline (", CFG$conf * 100, "% CI)"), y = NULL,
         title = "Accuracy beyond population structure (best GBLUP model per variable)") + theme_bw()
  save_fig("gain_over_baseline", p2, width = 13, height = 3 + 0.18 * ncol(Y))
  
  # Figure 3: SNP sets (whole vs sub-genomes)
  if (DO_SUB) {
    d3 <- acc %>% filter(model %in% CFG$light_models | model == "MK_AB")
    d3$snp_set <- factor(d3$snp_set, levels = c(names(snp_sets), "A+B"))
    p3 <- ggplot(d3, aes(snp_set, r, fill = snp_set)) + geom_hline(yintercept = 0, linetype = 2) +
      geom_boxplot(outlier.size = .5) + facet_grid(model ~ scheme) +
      labs(x = "SNP set", y = "Prediction accuracy (r)", title = "Whole genome vs sub-genomes",
           subtitle = "_matched = larger sub-genome subsampled to the SNP number of the smaller; A+B = multi-kernel model") +
      theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    save_fig("accuracy_by_subgenome", p3, width = 13, height = 2.2 + 2.2 * n_distinct(d3$model))
    # paired A vs B comparison (GBLUP_add), per scheme
    ab <- acc %>% filter(model == "GBLUP_add", snp_set %in% names(snp_sets)) %>%
      select(trait, scheme, snp_set, r, r_within) %>%
      tidyr::pivot_wider(names_from = snp_set, values_from = c(r, r_within))
    write.csv(ab, tab("subgenome_accuracy_paired.csv"), row.names = FALSE)
  }
}

# ==============================================================================
# 3  SUB-GENOME STRUCTURE AND GENOMIC VARIANCE PARTITION
# ==============================================================================
if (DO_SUB && CFG$run["subgenome_vc"]) {
  msg("== 3. Sub-genome variance partition ==")
  ut <- upper.tri(KS$whole$add)
  kc <- data.frame(comparison = c(combn(c("whole", sg_names), 2, paste, collapse = " vs "))) 
  kc$grm_correlation <- apply(combn(c("whole", sg_names), 2), 2, function(x) cor(KS[[x[1]]]$add[ut], KS[[x[2]]]$add[ut]))
  for (g in unique(grp)) {                       # within-group agreement between sub-genomes
    i <- which(grp == g); if (length(i) < 5) next; u <- upper.tri(diag(length(i)))
    kc[[paste0("within_", g)]] <- apply(combn(c("whole", sg_names), 2), 2, function(x)
      cor(KS[[x[1]]]$add[i, i][u], KS[[x[2]]]$add[i, i][u]))
  }
  write.csv(kc, tab("subgenome_GRM_agreement.csv"), row.names = FALSE)
  print(kc)
  
  vc_one <- function(trait) {
    set.seed(CFG$seed); y <- Y[, trait]
    f <- BGLR::BGLR(y, ETA = lapply(sg_names, function(s) list(K = KS[[s]]$add, model = "RKHS")),
                    nIter = CFG$vc_iter, burnIn = CFG$vc_burn, verbose = FALSE,
                    saveAt = file.path(tempdir(), paste0("vc_", Sys.getpid(), "_")))
    v <- sapply(f$ETA, `[[`, "varU"); names(v) <- sg_names
    w <- rrBLUP::mixed.solve(y, K = KS$whole$add)
    data.frame(trait, t(setNames(v, paste0("var_", sg_names))), var_residual = f$varE,
               t(setNames(v / sum(v), paste0("share_of_genomic_", sg_names))),
               h2_genomic_AplusB = sum(v) / (sum(v) + f$varE),
               h2_whole_REML = w$Vu / (w$Vu + w$Ve), check.names = FALSE)
  }
  vcd <- bind_rows(papply(colnames(Y), vc_one))
  for (s in sg_names) vcd[[paste0("share_of_SNPs_", s)]] <- length(snp_sets[[s]]) / sum(lengths(snp_sets[sg_names]))
  write.csv(vcd, tab("subgenome_variance_partition.csv"), row.names = FALSE)
  msg("Mean share of genomic variance:", paste(sg_names, round(colMeans(vcd[, paste0("share_of_genomic_", sg_names), drop = FALSE]), 3),
                                               sep = "=", collapse = ", "),
      "| share of SNPs:", paste(sg_names, round(unlist(vcd[1, paste0("share_of_SNPs_", sg_names)]), 3), sep = "=", collapse = ", "))
  vl <- tidyr::pivot_longer(vcd[, c("trait", paste0("var_", sg_names), "var_residual")], -trait,
                            names_to = "component", values_to = "variance")
  vl$component <- sub("^var_", "", vl$component)
  vl$trait <- factor(vl$trait, levels = vcd$trait[order(vcd[[paste0("share_of_genomic_", sg_names[1])]])])
  p4 <- ggplot(vl, aes(variance, trait, fill = component)) + geom_col(position = "fill") +
    geom_vline(xintercept = 0.5, linetype = 3) +
    labs(x = "Proportion of variance", y = NULL, fill = NULL,
         title = "Environment-of-origin variance explained by each sub-genome",
         subtitle = "Two-kernel Bayesian GBLUP (additive GRMs per sub-genome)") + theme_bw()
  save_fig("subgenome_variance_partition", p4, width = 8, height = 2.5 + 0.18 * ncol(Y))
}

# ==============================================================================
# 4  CORE COLLECTION AS TRAINING SET (core chosen a priori)
# ==============================================================================
if (CFG$run["core"]) {
  msg("== 4. Core collection as training set ==")
  core_size <- if (is.na(CFG$core_size)) round(0.3 * length(ids)) else CFG$core_size
  D_whole <- as.matrix(dist(Gc))
  read_core <- function(f) { L <- readLines(f, warn = FALSE); L <- L[nzchar(trimws(gsub(",", "", L)))]
  intersect(trimws(gsub('"', "", sapply(strsplit(L, ","), `[`, 1))), ids) }
  maxdiss <- function(idx, n) {
    Dm <- D_whole[idx, idx, drop = FALSE]; sel <- which.max(rowSums(Dm))
    while (length(sel) < min(n, length(idx))) sel <- c(sel, which.max(apply(Dm[, sel, drop = FALSE], 1, min)))
    idx[sel]
  }
  cores <- list()
  if (!is.na(CFG$core_csv) && file.exists(CFG$core_csv)) cores$supplied <- read_core(CFG$core_csv)
  if (requireNamespace("corehunter", quietly = TRUE)) {
    ch <- tryCatch({ set.seed(CFG$seed)
      gd <- corehunter::genotypes(round(G), format = "biparental")
      corehunter::sampleCore(corehunter::coreHunterData(genotypes = gd), size = core_size,
                             obj = list(corehunter::objective("EN", "MR"), corehunter::objective("AN", "MR"))) },
      error = function(e) { msg("Core Hunter failed:", conditionMessage(e)); NULL })
    if (!is.null(ch)) cores$corehunter_genotype <- intersect(ch$sel, ids)
  }
  cores$max_dissimilarity <- ids[maxdiss(seq_along(ids), core_size)]
  alloc <- table(grp); alloc[] <- pmax(2, round(core_size * alloc / length(ids)))   # every group represented
  cores$balanced_max_dissimilarity <- ids[unlist(lapply(names(alloc), function(g) maxdiss(which(grp == g), alloc[[g]])))]
  chosen <- if (CFG$core_rule %in% names(cores)) CFG$core_rule else "balanced_max_dissimilarity"
  if (chosen != CFG$core_rule) msg("Core rule", CFG$core_rule, "not available; using", chosen)
  
  comp_txt <- function(k) { tb <- table(factor(grp[match(cores[[k]], ids)], levels = sort(unique(grp))))
  paste(names(tb), tb, sep = "=", collapse = "; ") }
  missing_groups <- function(k) paste(setdiff(unique(grp), grp[match(cores[[k]], ids)]), collapse = ";")
  eval_core <- function(core_ids) {
    core_i <- match(core_ids, ids); nc <- length(core_i); comp <- table(grp[core_i])
    strat <- function() unlist(lapply(names(comp), function(g) { pool <- which(grp == g)
    pool[sample.int(length(pool), min(comp[[g]], length(pool)))] }))
    set.seed(CFG$seed)
    sets <- list(core = list(core_i),
                 random_matched = replicate(CFG$core_reps, sample.int(length(ids), nc), simplify = FALSE),
                 stratified_matched = replicate(CFG$core_reps, strat(), simplify = FALSE))
    out <- list()
    for (trait in colnames(Y)) for (mo in c("baseline", CFG$light_models)) for (sch in names(sets)) {
      st <- sapply(sets[[sch]], function(tr) { te <- setdiff(seq_along(ids), tr)
      point_stats(fit_pred(mo, "whole", Y[, trait], tr, te), Y[te, trait], grp[te]) })
      out[[paste(trait, mo, sch)]] <- data.frame(trait, model = mo, training = sch, n_train = nc,
                                                 r = mean(st["r", ]), r_sd = if (ncol(st) > 1) sd(st["r", ]) else NA, r_within = mean(st["r_within", ], na.rm = TRUE))
    }
    ca <- bind_rows(out)
    b <- ca %>% filter(model == "baseline") %>% select(trait, training, r_base = r)
    ca %>% left_join(b, by = c("trait", "training")) %>% mutate(gain = r - r_base)
  }
  res_core <- lapply(names(cores), function(k) { msg("Evaluating core:", k, "-", comp_txt(k)); eval_core(cores[[k]]) %>% mutate(core = k) })
  names(res_core) <- names(cores)
  cmp <- bind_rows(res_core) %>% filter(training == "core", model != "baseline") %>% group_by(core) %>%
    summarise(n = first(n_train), mean_r = mean(r, na.rm = TRUE), mean_r_within = mean(r_within, na.rm = TRUE),
              mean_gain = mean(gain, na.rm = TRUE), .groups = "drop") %>%
    mutate(composition = sapply(core, comp_txt), groups_missing = sapply(core, missing_groups),
           used_as_training_set = core == chosen)
  vs <- bind_rows(res_core) %>% filter(model != "baseline") %>% group_by(core, training) %>%
    summarise(mean_r = mean(r, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = training, values_from = mean_r, names_prefix = "mean_r_")
  cmp <- left_join(cmp, vs, by = "core")
  write.csv(cmp, tab("core_comparison_summary.csv"), row.names = FALSE)
  write.csv(bind_rows(res_core), tab("core_comparison_detail.csv"), row.names = FALSE)
  for (k in names(cores)) write.csv(meta[match(cores[[k]], meta$SampleID), c("SampleID", "Population", "Longitude", "Latitude")],
                                    tab(paste0("core_", k, ".csv")), row.names = FALSE)
  if (any(nzchar(cmp$groups_missing))) msg("Cores missing a group:", paste(cmp$core[nzchar(cmp$groups_missing)], collapse = ", "))
  print(cmp)
  ca <- res_core[[chosen]]
  lab_tr <- c(core = "Core", random_matched = "Random\n(size-matched)", stratified_matched = "Population-stratified\n(size-matched)")
  ca$training <- factor(ca$training, levels = names(lab_tr), labels = lab_tr)
  p5 <- ggplot(filter(ca, model != "baseline"), aes(training, r, fill = training)) + geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.size = .6) + facet_wrap(~model) +
    labs(x = NULL, y = "Accuracy on non-training samples (r)",
         title = paste0("Core '", chosen, "' (n = ", length(cores[[chosen]]), ", chosen a priori) vs size-matched training sets"),
         subtitle = "Whole genome; one point per environmental variable") +
    theme_bw() + theme(legend.position = "none", axis.text.x = element_text(size = 8))
  save_fig("core_vs_random_training", p5, width = 9, height = 4.5)
}

# ==============================================================================
# 5  WITHIN-POPULATION GEAVs AND CANDIDATE PARENTS (predictable variables only)
# ==============================================================================
if (CFG$run["geav"] && !is.null(acc)) {
  msg("== 5. GEAVs ==")
  # Model per variable: whole-genome GBLUP model with the highest within-group accuracy
  # under random CV. A variable is "predictable" if that accuracy's CI excludes zero.
  best <- acc %>% filter(scheme == RANDOM, snp_set == "whole", model %in% CFG$light_models) %>%
    group_by(trait) %>% slice_max(r_within, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(predictable = r_within_lo > 0)
  write.csv(best %>% select(trait, model, r, r_within, r_within_lo, r_within_hi, gain_within, predictable),
            tab("GEAV_model_per_variable.csv"), row.names = FALSE)
  msg("Variables predictable within groups:", sum(best$predictable), "of", nrow(best))
  geav_from <- function(set) {
    m <- sapply(colnames(Y), function(t) OOF[[t]][[paste(best$model[best$trait == t], set)]])
    rownames(m) <- ids; m <- apply(m, 2, function(v) v - ave(v, grp)); scale(m)
  }
  geav_within <- geav_from("whole")
  out <- data.frame(SampleID = ids, Population = grp, Region = meta$Region, Country = meta$Country, geav_within, check.names = FALSE)
  write.csv(out, tab("GEAV_within_population.csv"), row.names = FALSE)
  if (DO_SUB) for (s in sg_names) {
    gs <- geav_from(s)
    write.csv(data.frame(SampleID = ids, Population = grp, gs, check.names = FALSE),
              tab(paste0("GEAV_within_population_subgenome_", s, ".csv")), row.names = FALSE)
  }
  if (DO_SUB) {                                  # do A- and B-based GEAVs agree?
    agree <- data.frame(trait = colnames(Y), predictable = best$predictable[match(colnames(Y), best$trait)])
    sa <- geav_from(sg_names[1]); sb <- geav_from(sg_names[2])
    agree[[paste0("cor_GEAV_", sg_names[1], "_", sg_names[2])]] <- sapply(colnames(Y), function(t) safe_cor(sa[, t], sb[, t]))
    write.csv(agree, tab("GEAV_subgenome_agreement.csv"), row.names = FALSE)
  }
  
  cand <- list()
  for (tn in names(CFG$targets)) {
    tg <- CFG$targets[[tn]]
    if (!tg$var %in% colnames(Y)) { msg("Target", tn, "skipped:", tg$var, "not among the variables"); next }
    b <- best[best$trait == tg$var, ]
    if (!isTRUE(b$predictable)) {
      msg("Target", tn, "skipped:", tg$var, "is not predictable within groups (within-group r CI includes 0)")
      cand[[tn]] <- data.frame(target = tn, description = tg$label, variable = tg$var, SampleID = NA, Population = NA,
                               within_population_GEAV = NA, r_within = b$r_within, note = "not predictable: no parents ranked")
      next
    }
    o <- if (tg$direction == "low") order(out[[tg$var]]) else order(-out[[tg$var]])
    top <- head(o, CFG$n_candidates)
    cand[[tn]] <- data.frame(target = tn, description = tg$label, variable = tg$var, out[top, c("SampleID", "Population", "Region", "Country")],
                             within_population_GEAV = out[[tg$var]][top], r_within = b$r_within,
                             note = sprintf("model %s; within-group r %.2f (%.2f-%.2f)", b$model, b$r_within, b$r_within_lo, b$r_within_hi))
  }
  cand_df <- bind_rows(cand); rownames(cand_df) <- NULL
  write.csv(cand_df, tab("candidate_parents_within_population.csv"), row.names = FALSE)
  print(cand_df)
  
  pv <- best$trait[best$predictable]
  if (length(pv) >= 2) {
    gw <- geav_within[, pv, drop = FALSE]
    ord <- unlist(lapply(split(seq_along(ids), grp), function(i) if (length(i) < 3) i else i[hclust(dist(gw[i, , drop = FALSE]))$order]))
    subcl <- setNames(rep(NA_character_, length(ids)), ids)
    for (i in split(seq_along(ids), grp)) {
      s <- if (length(i) < 3) rep(1, length(i)) else cutree(hclust(dist(gw[i, , drop = FALSE])), min(CFG$n_subclusters, length(i)))
      subcl[i] <- paste0(grp[i[1]], "_", s)
    }
    cand_ids <- unique(na.omit(cand_df$SampleID))
    trait_ord <- colnames(gw)[hclust(dist(t(gw)))$order]
    
    # whole-panel clustering of within-population GEAVs (ignores population)
    hc_all <- hclust(dist(gw)); k_all <- length(unique(grp))
    gcl <- paste0("C", cutree(hc_all, k_all)); meta$GEAV_cluster <- gcl
    sc <- data.frame(SampleID = ids, Population = grp, Region = meta$Region, Country = meta$Country,
                     Subcluster_within_population = subcl[ids], Cluster_all_samples = gcl)
    write.csv(sc, tab("GEAV_subclusters.csv"), row.names = FALSE)
    comp <- bind_rows(
      sc %>% count(clustering = "within_population", cluster = Subcluster_within_population, Region, Country),
      sc %>% count(clustering = "all_samples", cluster = Cluster_all_samples, Region, Country))
    write.csv(comp, tab("GEAV_cluster_composition_by_country_region.csv"), row.names = FALSE)
    print(table(Cluster = gcl, Region = meta$Region))
    
    # colours: regions, then countries grouped by region
    regs <- sort(unique(meta$Region)); pops <- sort(unique(grp))
    reg_cols <- setNames(grDevices::hcl.colors(max(3, length(regs)), "Dark 3")[seq_along(regs)], regs)
    pop_cols <- setNames(grDevices::hcl.colors(max(3, length(pops)), "Set 2")[seq_along(pops)], pops)
    ctry <- unique(meta[order(meta$Region, meta$Country), c("Region", "Country")]$Country)
    ctry_cols <- setNames(grDevices::hcl.colors(length(ctry) + 1, "Spectral")[seq_along(ctry)], ctry)
    if ("Unknown" %in% ctry) ctry_cols["Unknown"] <- "grey80"
    
    geav_heatmap <- function(row_ids, facet, annots, title, subtitle, name) {
      fac <- if (facet) facet_grid(Population ~ ., scales = "free_y", space = "free_y", switch = "y") else NULL
      ylev <- row_ids
      base_df <- meta[match(row_ids, ids), c("SampleID", "Population", "Region", "Country", "GEAV_cluster")]
      base_df$SampleID <- factor(base_df$SampleID, levels = ylev)
      th_ann <- theme_minimal(base_size = 9) +
        theme(axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1), panel.spacing.y = unit(2, "pt"),
              legend.key.size = unit(8, "pt"))
      ann_plot <- function(v, cols, first) {
        d <- base_df; d$x <- v; d$val <- d[[v]]
        ggplot(d, aes(x, SampleID, fill = val)) + geom_tile() + fac +
          scale_fill_manual(values = cols, name = v, guide = guide_legend(ncol = if (length(cols) > 15) 2 else 1)) +
          labs(x = NULL, y = NULL) + th_ann +
          theme(strip.text.y.left = if (first) element_text(angle = 0, hjust = 1) else element_blank(),
                strip.placement = "outside")
      }
      cl_lev <- sort(unique(gcl))
      cols_for <- list(Population = pop_cols, Region = reg_cols, Country = ctry_cols,
                       GEAV_cluster = setNames(grDevices::hcl.colors(max(3, length(cl_lev)), "Dark 2")[seq_along(cl_lev)], cl_lev))
      ps <- lapply(seq_along(annots), function(i) ann_plot(annots[i], cols_for[[annots[i]]], first = facet && i == 1))
      hm <- data.frame(SampleID = factor(rep(row_ids, ncol(gw)), levels = ylev),
                       Population = rep(grp[match(row_ids, ids)], ncol(gw)),
                       trait = factor(rep(colnames(gw), each = length(row_ids)), levels = trait_ord),
                       GEAV = pmax(-3, pmin(3, as.vector(gw[match(row_ids, ids), , drop = FALSE]))))
      ph <- ggplot(hm, aes(trait, SampleID, fill = GEAV)) + geom_tile() + fac +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-3, 3),
                             name = "Within-population\nGEAV (z, capped \u00b13)") +
        scale_y_discrete(position = "right", labels = function(x) ifelse(x %in% cand_ids, paste0(x, " *"), x)) +
        labs(x = NULL, y = NULL) + theme_minimal(base_size = 9) +
        theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
              axis.text.y = element_text(size = max(2, min(7, 400 / length(ids)))),
              strip.text.y = element_blank(), panel.grid = element_blank(), panel.spacing.y = unit(2, "pt"))
      p <- patchwork::wrap_plots(c(ps, list(ph)), nrow = 1, widths = c(rep(0.6, length(ps)), max(4, ncol(gw)))) +
        patchwork::plot_layout(guides = "collect") +
        patchwork::plot_annotation(title = title, subtitle = subtitle)
      save_fig(name, p, width = 7 + 0.35 * ncol(gw) + 0.6 * length(ps) + if (length(ctry) > 15) 2 else 0,
               height = min(45, 3 + length(ids) * 0.09))
    }
    # (1) grouped by population, clustered within population
    geav_heatmap(ids[ord], facet = TRUE, annots = c("Region", "Country"),
                 title = paste0(CFG$species, ": within-population GEAVs"),
                 subtitle = "Rows clustered within population; variables predictable within groups only; * = candidate parent",
                 name = "GEAV_within_population_heatmap")
    # (2) all samples clustered together, to see whether GEAV clusters follow population, region or country
    geav_heatmap(ids[hc_all$order], facet = FALSE, annots = c("GEAV_cluster", "Population", "Region", "Country"),
                 title = paste0(CFG$species, ": within-population GEAVs, all samples clustered together"),
                 subtitle = "Rows ordered by hierarchical clustering of GEAVs across all samples (GEAV_cluster = tree cut into as many clusters as populations); * = candidate parent",
                 name = "GEAV_heatmap_clustered_all_samples")
  } else msg("Fewer than two predictable variables: GEAV heatmap skipped")
}

writeLines(capture.output(sessionInfo()), file.path(CFG$out_dir, "sessionInfo.txt"))
msg("Done. Results in", normalizePath(CFG$out_dir))