# ==== 04_environmental_GWAS.R - CENTROID-FILTERED PIPELINE (one accession per putative country centroid) ====
# Environmental GWAS (6 models; per variable the model with lambda closest to 1). Output significant_markers.csv is used by 05 and 06.
# Inputs are read from data_filtered/ (created by 00_filter_centroid_accessions.R).

################################################################################
# Environmental GWAS (eGWAS): analysis and visualization in one script
#
# Combines environmental_gwas_from_vcf.R and visualize_egwas_results.R.
#
# Inputs
#   1. A VCF of biallelic SNPs.
#   2. A CSV of environmental data: first column = sample ID (matching the VCF),
#      optional latitude / longitude columns, then one column per environmental variable.
#
# What it does
#   - Genotype QC (missingness, MAF), kinship (rrBLUP::A.mat), PCs.
#   - Six mixed models with different structure corrections (K, K + 3 PCs, and
#     latitude / longitude covariates), fitted with rrBLUP::GWAS (P3D).
#   - Genomic inflation (lambda) per variable x model; QQ plots with ALL models on
#     the same graph (one panel per variable).
#   - Significant markers by FDR (Storey q-value, BH fallback), with, for every
#     significant marker:
#       beta          additive effect per ALT allele (GLS under the mixed model)
#       PVE_percent   % of the environmental variable's variance explained by the marker:
#                     100 * beta^2 * var(genotype) / var(variable)
#       partial_R2_percent  % of the structure-adjusted residual variance explained
#                     (t^2 / (t^2 + df) from the mixed-model GLS fit)
#       Single-marker PVEs are not additive: markers in LD (or sharing structure) overlap.
#   - Independent loci: significant markers clumped by LD window / r2 (columns locus, lead_marker);
#     the lead markers are fitted jointly and variance_explained.csv reports their adjusted and
#     10-fold cross-validated R2 (% of the residual variance after structure and kinship). Fitting
#     ALL significant markers jointly overfits (in-sample values can exceed 100%).
#       R2_raw_percent      % variance explained by a simple regression, no structure
#                     correction (shown for comparison; inflated by structure)
#   - Joint variance partition per variable (structure, polygenic, markers).
#   - Manhattan plots, summary figures, allele-environment boxplots, and a
#     geographic allele map drawn on a world basemap (land, country borders,
#     ocean, continent labels).
#   - Figures as PDF (vector) and PNG; tables as CSV; an HTML summary report.
#
# Packages: rrBLUP, dplyr, tidyr, purrr, tibble, ggplot2, patchwork, maps
#           (qvalue optional; installed from Bioconductor if possible)
################################################################################

# ==============================================================================
# 0  CONFIGURATION
# ==============================================================================
vcf_file      <- "data_filtered/finger_millet_filtered.vcf"     # VCF (.vcf or .vcf.gz)
env_file      <- "data_filtered/millet_extracted_bioclim_values.csv"   # sample ID + lat/lon + variables
coords_file   <- "data_filtered/millet_passport_test.csv"                              # optional: separate CSV with sample ID + latitude +
                                                  # longitude (used if env_file has no coordinates)
lat_col       <- NULL                             # NULL = auto-detect; or the exact column name
lon_col       <- NULL
output_dir    <- "results/04_environmental_GWAS/output"                         # tables and RData
figures_dir   <- "results/04_environmental_GWAS/figures"                        # all figures + HTML report

min_maf        <- 0.05     # minimum minor allele frequency
missing_thresh <- 0.10     # maximum proportion missing per marker
fdr_threshold  <- 0.20     # FDR threshold for significance
# Joint variance of significant markers (avoids overfitting): markers are clumped into independent
# loci (within clump_window of a stronger marker on the same chromosome, or r2 >= clump_r2 with it);
# one lead marker per locus is fitted jointly (at most max_lead_frac x n), and adjusted and
# cross-validated R2 are reported instead of the in-sample value.
clump_window   <- c(A = 100200, B = 190200, genome = 121200)   # LD decay (r2_VS = 0.2), filtered panel
clump_r2       <- 0.5
max_lead_frac  <- 0.10
cv_folds       <- 10
sig_model      <- "best"   # model used to call significant markers and estimate PVE:
                           # "model1".."model6", or "best" = per-variable model with
                           # lambda closest to 1
n_top_markers  <- 15       # pleiotropic markers shown in allele-environment plots
n_map_markers  <- 5        # markers shown on the geographic map
fig_png_dpi    <- 300      # PNG copies (used in the HTML report); PDFs are vector

# Models: structure correction for each
egwas_models <- tibble::tribble(
  ~model,   ~n_PCs, ~covariates,                  ~description,
  "model1", 0,      character(),                  "K",
  "model2", 3,      character(),                  "K + 3 PCs",
  "model3", 0,      c("latitude"),                "K + lat",
  "model4", 3,      c("latitude"),                "K + 3 PCs + lat",
  "model5", 0,      c("latitude", "longitude"),   "K + lat + lon",
  "model6", 3,      c("latitude", "longitude"),   "K + 3 PCs + lat + lon"
)

# ==============================================================================
# 1  PACKAGES AND HELPERS
# ==============================================================================
need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
}
need(c("rrBLUP", "dplyr", "tidyr", "purrr", "tibble", "ggplot2", "patchwork", "maps"))
if (!requireNamespace("qvalue", quietly = TRUE)) {
  try({ need("BiocManager"); BiocManager::install("qvalue", update = FALSE, ask = FALSE) }, silent = TRUE)
}
HAS_QVALUE <- requireNamespace("qvalue", quietly = TRUE)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble); library(ggplot2); library(patchwork)
})
# Note: dplyr verbs are called with dplyr:: where other packages commonly mask them.

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
theme_set(theme_bw())
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

save_fig <- function(name, plot, width, height) {          # PDF (vector) + PNG
  dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(figures_dir, paste0(name, ".pdf")), plot, width = width, height = height,
         device = dev, bg = "white", limitsize = FALSE)
  ggsave(file.path(figures_dir, paste0(name, ".png")), plot, width = width, height = height,
         dpi = fig_png_dpi, bg = "white", limitsize = FALSE)
}
genomic_lambda <- function(p) {
  p <- p[!is.na(p) & p > 0 & p <= 1]
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
}
genotype_class <- function(g) factor(ifelse(g <= -0.5, "Ref/Ref", ifelse(g >= 0.5, "Alt/Alt", "Het")),
                                     levels = c("Ref/Ref", "Het", "Alt/Alt"))
geno_cols  <- c("Ref/Ref" = "#3498db", "Het" = "#95a5a6", "Alt/Alt" = "#e74c3c")
geno_shape <- c("Ref/Ref" = 16, "Het" = 17, "Alt/Alt" = 15)

# ==============================================================================
# 2  GENOTYPES
# ==============================================================================
msg("Reading VCF:", vcf_file)
read_vcf_numeric <- function(path) {               # samples x markers, coded -1/0/1 (rrBLUP)
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
  bi <- !grepl(",", v[[5]])                        # biallelic only
  v <- v[bi, , drop = FALSE]
  gt <- sub(":.*$", "", as.matrix(v[, -(1:9), drop = FALSE]))
  a1 <- sub("[/|].*$", "", gt); a2 <- sub("^[^/|]*[/|]", "", gt)
  d <- (a1 != "0") + (a2 != "0") - 1; d[a1 == "." | a2 == "."] <- NA
  storage.mode(d) <- "numeric"
  id <- ifelse(v[[3]] %in% c(".", ""), paste(v[[1]], v[[2]], sep = "_"), v[[3]])
  id <- make.unique(id)
  G <- t(d); rownames(G) <- trimws(samples); colnames(G) <- id
  list(G = G, info = data.frame(marker = id, chrom = v[[1]], pos = as.integer(v[[2]]), stringsAsFactors = FALSE))
}
vc <- read_vcf_numeric(vcf_file)
geno_mat_t <- vc$G; snp_info <- vc$info
poly <- apply(geno_mat_t, 2, function(x) length(unique(x[!is.na(x)])) > 1)
geno_mat_t <- geno_mat_t[, poly, drop = FALSE]; snp_info <- snp_info[poly, ]
msg("VCF:", ncol(geno_mat_t), "polymorphic biallelic markers x", nrow(geno_mat_t), "samples")

# ==============================================================================
# 3  ENVIRONMENTAL DATA
# ==============================================================================
env_data <- read.csv(env_file, stringsAsFactors = FALSE, check.names = FALSE)
names(env_data)[1] <- "sample_id"; env_data$sample_id <- trimws(env_data$sample_id)
# Latitude / longitude: accepts lat, Lat, LATITUDE, Latitude (deg), lat_dd, decimalLatitude,
# DECLATITUDE (genebank passport data), y/x are NOT assumed.
find_coord <- function(nms, type, user) {
  if (!is.null(user)) { if (!user %in% nms) stop("Column '", user, "' not found"); return(user) }
  n <- tolower(gsub("[^A-Za-z]", "", nms))
  pat <- if (type == "lat") "^(decimal|dec)?lat(itude)?(dd|deg|degrees|decimal)?$"
         else "^(decimal|dec)?(lon|long|longitude)(dd|deg|degrees|decimal)?$"
  nms[grepl(pat, n)][1]
}
standardise_coords <- function(d, lat_user, lon_user) {
  la <- find_coord(names(d), "lat", lat_user); lo <- find_coord(names(d), "lon", lon_user)
  if (!is.na(la)) names(d)[names(d) == la] <- "latitude"
  if (!is.na(lo)) names(d)[names(d) == lo] <- "longitude"
  for (cc in intersect(c("latitude", "longitude"), names(d)))      # text -> numeric ("12,5" -> 12.5)
    if (!is.numeric(d[[cc]])) d[[cc]] <- suppressWarnings(as.numeric(gsub(",", ".", trimws(d[[cc]]))))
  list(d = d, found = c(lat = la, lon = lo))
}
sc <- standardise_coords(env_data, lat_col, lon_col); env_data <- sc$d
if (!all(c("latitude", "longitude") %in% names(env_data)) && !is.na(coords_file)) {
  cf <- read.csv(coords_file, stringsAsFactors = FALSE, check.names = FALSE)
  names(cf)[1] <- "sample_id"; cf$sample_id <- trimws(cf$sample_id)
  cs <- standardise_coords(cf, lat_col, lon_col)
  if (!all(c("latitude", "longitude") %in% names(cs$d))) stop("No latitude/longitude columns found in ", coords_file)
  env_data <- env_data %>% dplyr::select(-any_of(c("latitude", "longitude"))) %>%
    dplyr::left_join(cs$d[, c("sample_id", "latitude", "longitude")], by = "sample_id")
  msg("Coordinates taken from", coords_file)
}
HAS_COORDS <- all(c("latitude", "longitude") %in% names(env_data))
if (HAS_COORDS) {
  msg("Coordinates: latitude and longitude found;", sum(is.na(env_data$latitude) | is.na(env_data$longitude)),
      "sample(s) without coordinates")
} else {
  msg("WARNING: no latitude/longitude columns recognised. Column names in", env_file, ":",
      paste(names(env_data), collapse = ", "))
  msg("WARNING: models with latitude/longitude covariates (model3-model6) cannot run.",
      "Set lat_col / lon_col, or give coords_file.")
}
env_vars <- setdiff(names(env_data), c("sample_id", "latitude", "longitude"))
env_vars <- env_vars[vapply(env_data[env_vars], is.numeric, logical(1))]
msg("Environmental variables:", length(env_vars), "|", paste(head(env_vars, 6), collapse = ", "),
    if (length(env_vars) > 6) "..." else "", "| coordinates:", HAS_COORDS)

common <- intersect(rownames(geno_mat_t), env_data$sample_id)
msg("Samples with genotype and environment data:", length(common))
if (!length(common)) stop("No common samples between VCF and environmental file; check IDs.")
geno_mat_t <- geno_mat_t[common, , drop = FALSE]
env_data <- env_data[match(common, env_data$sample_id), ]

# ==============================================================================
# 4  QC, IMPUTATION, KINSHIP, PCs
# ==============================================================================
miss <- colMeans(is.na(geno_mat_t))
keep <- miss <= missing_thresh
p_alt <- (colMeans(geno_mat_t, na.rm = TRUE) + 1) / 2
maf <- pmin(p_alt, 1 - p_alt)
keep <- keep & maf >= min_maf
msg("QC: removed", sum(miss > missing_thresh), "markers for missingness and",
    sum(miss <= missing_thresh & maf < min_maf), "for MAF <", min_maf)
geno_mat_final <- geno_mat_t[, keep, drop = FALSE]
snp_info_final <- snp_info[keep, ]
maf_final <- maf[keep]; p_alt_final <- p_alt[keep]
for (j in which(colSums(is.na(geno_mat_final)) > 0)) {
  x <- geno_mat_final[, j]; x[is.na(x)] <- mean(x, na.rm = TRUE); geno_mat_final[, j] <- x
}
msg("Final:", nrow(geno_mat_final), "samples x", ncol(geno_mat_final), "markers")

K <- rrBLUP::A.mat(geno_mat_final, min.MAF = 0, max.missing = 1)
eK <- eigen(K, symmetric = TRUE)
K_eigens <- data.frame(sample_id = rownames(K), PC1 = eK$vectors[, 1], PC2 = eK$vectors[, 2], PC3 = eK$vectors[, 3])
msg("Kinship", nrow(K), "x", ncol(K), "| PC1 explains", round(100 * eK$values[1] / sum(eK$values), 2), "%")

geno_df <- data.frame(snp_info_final, t(geno_mat_final), check.names = FALSE)
pheno <- env_data %>% dplyr::select(sample_id, all_of(env_vars), any_of(c("latitude", "longitude"))) %>%
  dplyr::rename(individual = sample_id) %>% as.data.frame()
for (v in env_vars) if (anyNA(pheno[[v]])) {
  msg("Missing values in", v, "- imputed with the median")
  pheno[[v]][is.na(pheno[[v]])] <- median(pheno[[v]], na.rm = TRUE)
}

# Phenotype as analysed by a model: residuals on that model's covariates
model_pheno <- function(var, covars) {
  y <- pheno[[var]]
  cv <- intersect(setdiff(covars, var), names(pheno))
  if (length(cv)) y <- as.numeric(residuals(lm(reformulate(cv, response = var), data = pheno,
                                                  na.action = na.exclude)))   # NA where a covariate is missing
  y
}

# ==============================================================================
# 5  ENVIRONMENTAL GWAS
# ==============================================================================
msg("Running eGWAS:", nrow(egwas_models), "models x", length(env_vars), "variables")
gwas_list <- list()
for (i in seq_len(nrow(egwas_models))) {
  mod <- egwas_models$model[i]; npc <- egwas_models$n_PCs[i]; cov <- egwas_models$covariates[[i]]
  if (length(cov) && !all(cov %in% names(pheno))) {
    msg("  ", mod, "skipped: covariates", paste(setdiff(cov, names(pheno)), collapse = ", "), "not available"); next
  }
  msg("  ", mod, "(", egwas_models$description[i], ")")
  for (v in env_vars) {
    ph <- data.frame(individual = pheno$individual, y = model_pheno(v, cov)); names(ph)[2] <- v
    r <- tryCatch(rrBLUP::GWAS(pheno = ph, geno = geno_df, n.PC = npc, min.MAF = 0, P3D = TRUE, plot = FALSE, K = K),
                  error = function(e) { msg("    error in", v, ":", conditionMessage(e)); NULL })
    if (!is.null(r)) gwas_list[[paste(mod, v)]] <- data.frame(model = mod, variable = v, marker = r$marker,
                                                              chrom = r$chrom, pos = r$pos, score = r[[v]])
  }
}
ran <- unique(sub(" .*$", "", names(gwas_list)))
not_ran <- setdiff(egwas_models$model, ran)
msg("Models run:", paste(ran, collapse = ", "))
if (length(not_ran)) msg("WARNING: models NOT run:", paste(not_ran, collapse = ", "), "- see messages above")
egwas_models <- egwas_models[egwas_models$model %in% ran, ]
if (sig_model != "best" && !sig_model %in% ran) stop("sig_model '", sig_model, "' did not run")
gwas_analyzed <- dplyr::bind_rows(gwas_list) %>%
  dplyr::mutate(score = ifelse(is.na(score), 0, score), p_value = 10^-score)
gwas_combined <- gwas_analyzed
model_labels <- setNames(paste0(egwas_models$model, ": ", egwas_models$description), egwas_models$model)

# ---- genomic inflation ----
p_lambda <- gwas_analyzed %>% dplyr::group_by(variable, model) %>%
  dplyr::summarise(lambda = genomic_lambda(p_value), .groups = "drop") %>%
  dplyr::group_by(variable) %>% dplyr::mutate(best_model = model[which.min(abs(1 - lambda))], is_best = model == best_model) %>%
  dplyr::ungroup()
msg("Best model per variable (lambda closest to 1):")
print(p_lambda %>% dplyr::filter(is_best) %>% dplyr::select(variable, model, lambda), n = Inf)

# ==============================================================================
# 6  QQ PLOTS: all models on the same graph, one panel per variable
# ==============================================================================
msg("QQ plots")
qq_df <- gwas_analyzed %>% dplyr::filter(!is.na(p_value), p_value > 0) %>%
  dplyr::group_by(variable, model) %>% dplyr::arrange(p_value, .by_group = TRUE) %>%
  dplyr::mutate(observed = -log10(p_value), expected = -log10(ppoints(dplyr::n()))) %>%
  # thin the dense bulk (p > 0.01) to keep the PDF small; tail kept in full
  dplyr::filter(p_value < 0.01 | dplyr::row_number() %% max(1, floor(dplyr::n() / 2000)) == 0) %>% dplyr::ungroup()
n_mark <- dplyr::n_distinct(gwas_analyzed$marker)
ci <- dplyr::tibble(i = seq_len(n_mark)) %>%
  dplyr::mutate(expected = -log10(ppoints(n_mark)),
         lo = -log10(qbeta(0.975, i, n_mark - i + 1)), hi = -log10(qbeta(0.025, i, n_mark - i + 1))) %>%
  dplyr::filter(dplyr::row_number() <= 2000 | dplyr::row_number() %% max(1, floor(n_mark / 2000)) == 0)
lam_lab <- p_lambda %>% dplyr::group_by(variable) %>% dplyr::arrange(model, .by_group = TRUE) %>%
  dplyr::summarise(label = paste0(sub("model", "m", model), " \u03bb=", sprintf("%.2f", lambda), ifelse(is_best, "*", ""),
                           collapse = "\n"), .groups = "drop")
n_vars <- length(env_vars); n_cols <- min(4, n_vars); n_rows <- ceiling(n_vars / n_cols)
p_qq <- ggplot(qq_df, aes(expected, observed, colour = model)) +
  geom_ribbon(data = ci, aes(x = expected, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey88") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_point(size = 0.6, alpha = 0.7) +
  geom_text(data = lam_lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.1, size = 2.7, lineheight = 0.9) +
  facet_wrap(~ variable, ncol = n_cols, scales = "free_y") +
  scale_colour_brewer(palette = "Dark2", labels = model_labels, name = "Model") +
  labs(x = expression(Expected ~ -log[10](italic(p))), y = expression(Observed ~ -log[10](italic(p))),
       title = "QQ plots: all models on the same graph",
       subtitle = "Grey band = 95% CI under the null; * = model with lambda closest to 1") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
save_fig("egwas_qqplots_all_models", p_qq, width = 3.6 * n_cols + 1, height = 3.3 * n_rows + 1.5)

# ---- lambda figures ----
p_lam_box <- ggplot(p_lambda, aes(variable, lambda)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "red") +
  geom_hline(yintercept = c(0.95, 1.05), linetype = 3, colour = "orange") +
  geom_point(aes(colour = model), size = 2.3, position = position_dodge(width = 0.6)) +
  scale_colour_brewer(palette = "Dark2", labels = model_labels, name = "Model") +
  labs(title = "Genomic inflation factor (\u03bb) by variable and model", x = NULL, y = "\u03bb") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
save_fig("lambda_by_model", p_lam_box, width = max(8, 0.45 * n_vars + 4), height = 6)
p_lam_heat <- ggplot(p_lambda, aes(model, variable, fill = lambda)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", lambda)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, name = "\u03bb") +
  labs(title = "Genomic inflation factor", x = NULL, y = NULL) + theme_minimal() + theme(panel.grid = element_blank())
save_fig("lambda_heatmap", p_lam_heat, width = 8, height = max(5, 0.35 * n_vars + 2))
p_best <- p_lambda %>% dplyr::filter(is_best) %>% dplyr::count(model) %>%
  ggplot(aes(model, n, fill = model)) + geom_col() + geom_text(aes(label = n), vjust = -0.4) +
  scale_fill_brewer(palette = "Dark2") + labs(title = "Model with \u03bb closest to 1", x = NULL, y = "Variables") +
  theme(legend.position = "none")
save_fig("best_model_summary", p_best, width = 7, height = 5)

# ==============================================================================
# 7  SIGNIFICANT MARKERS
# ==============================================================================
msg("Significant markers: FDR <=", fdr_threshold, "| model:", sig_model)
model_for <- function(v) if (sig_model == "best") p_lambda$best_model[p_lambda$variable == v][1] else sig_model
sig_all <- list()
for (v in env_vars) {
  d <- gwas_analyzed %>% dplyr::filter(variable == v, model == model_for(v), !is.na(p_value))
  if (!nrow(d)) next
  q <- if (HAS_QVALUE) tryCatch(qvalue::qvalue(d$p_value)$qvalues, error = function(e) NULL) else NULL
  d$q_value <- if (is.null(q)) p.adjust(d$p_value, "BH") else q
  d$fdr_method <- if (is.null(q)) "Benjamini-Hochberg" else "Storey q-value"
  s <- d %>% dplyr::filter(q_value <= fdr_threshold) %>% dplyr::arrange(p_value)
  if (nrow(s)) msg("  ", v, ":", nrow(s), "significant markers [", d$fdr_method[1], "]")
  sig_all[[v]] <- s
}
sig_markers <- dplyr::bind_rows(sig_all)
msg("Total significant marker-variable associations:", nrow(sig_markers))

# ==============================================================================
# 8  PERCENT VARIANCE EXPLAINED (per significant marker) AND JOINT PARTITION
# ==============================================================================
# Null mixed model y = Xb + u + e, u ~ N(0, Vu K), with the same fixed effects as the
# GWAS model (intercept + n PCs of K; covariates already regressed out of y).
# Rotating by H^(-1/2) (H = Vu K + Ve I) gives GLS estimates of each marker's effect.
null_fit <- function(v, mod) {
  i <- match(mod, egwas_models$model)
  y_an <- model_pheno(v, egwas_models$covariates[[i]]); ok <- !is.na(y_an)   # drop samples lacking covariates
  X <- cbind(1, eK$vectors[ok, seq_len(egwas_models$n_PCs[i]), drop = FALSE])
  Ko <- K[ok, ok]; e <- if (all(ok)) eK else eigen(Ko, symmetric = TRUE)
  f <- rrBLUP::mixed.solve(y_an[ok], K = Ko, X = X, method = "REML")
  w <- 1 / sqrt(f$Vu * pmax(e$values, 0) + f$Ve)
  rot <- function(M) e$vectors %*% (w * crossprod(e$vectors, M))       # H^(-1/2) M
  list(ok = ok, y_raw = pheno[[v]][ok], y_an = y_an[ok], X = X, fit = f, ys = rot(y_an[ok]), Xs = rot(X), rot = rot, Kd = mean(diag(Ko)))
}
marker_pve <- function(nf, markers) {
  out <- lapply(markers, function(m) {
    g <- geno_mat_final[nf$ok, m] + 1                  # ALT dosage 0/1/2
    xs <- nf$rot(g); A <- cbind(nf$Xs, xs)
    fit <- lm.fit(A, nf$ys); df <- length(nf$ys) - ncol(A)
    s2 <- sum(fit$residuals^2) / df
    se <- sqrt(s2 * diag(chol2inv(qr.R(fit$qr)))[ncol(A)])
    b <- fit$coefficients[ncol(A)]; tval <- b / se
    data.frame(marker = m, alt_allele_freq = (mean(g)) / 2, beta = b, beta_se = se,
               PVE_percent = 100 * b^2 * var(g) / var(nf$y_raw),
               partial_R2_percent = 100 * tval^2 / (tval^2 + df),
               R2_raw_percent = 100 * summary(lm(nf$y_raw ~ g))$r.squared)
  })
  dplyr::bind_rows(out)
}
# ---- independent loci (clumping) and honest joint variance for one variable ----
# d: significant markers of one variable (marker, chrom, pos, p_value); G: samples x markers (ALT dosage)
clump_loci <- function(d, G, windows, r2_max) {
  d <- d[order(d$p_value), ]; locus <- rep(NA_integer_, nrow(d)); k <- 0
  sg <- function(ch) { n <- suppressWarnings(as.integer(ch)); ifelse(is.na(n) | n > 18, "genome", ifelse(n %% 2 == 1, "A", "B")) }
  for (i in seq_len(nrow(d))) if (is.na(locus[i])) {
    k <- k + 1; locus[i] <- k
    same <- which(is.na(locus) & d$chrom == d$chrom[i])
    if (!length(same)) next
    near <- abs(d$pos[same] - d$pos[i]) <= windows[[sg(d$chrom[i])]]
    r2 <- suppressWarnings(cor(G[, d$marker[i]], G[, d$marker[same], drop = FALSE])^2); r2[is.na(r2)] <- 0
    locus[same[near | as.vector(r2) >= r2_max]] <- k
  }
  d$locus <- locus; d$lead_marker <- !duplicated(locus); d
}
joint_r2 <- function(ys, Xs, Gs, folds = 10, seed = 1) {
  n <- length(ys); p0 <- ncol(Xs); k <- ncol(Gs)
  r0 <- lm.fit(Xs, ys)$residuals; r1 <- lm.fit(cbind(Xs, Gs), ys)$residuals
  R2 <- 1 - sum(r1^2) / sum(r0^2); R2adj <- 1 - (1 - R2) * (n - p0) / (n - p0 - k)
  set.seed(seed); f <- sample(rep(seq_len(folds), length.out = n)); e0 <- e1 <- numeric(n)
  for (j in seq_len(folds)) { tr <- f != j; te <- !tr
    b0 <- lm.fit(Xs[tr, , drop = FALSE], ys[tr])$coefficients; b0[is.na(b0)] <- 0
    A <- cbind(Xs, Gs); b1 <- lm.fit(A[tr, , drop = FALSE], ys[tr])$coefficients; b1[is.na(b1)] <- 0
    e0[te] <- ys[te] - Xs[te, , drop = FALSE] %*% b0; e1[te] <- ys[te] - A[te, , drop = FALSE] %*% b1 }
  c(in_sample = 100 * R2, adjusted = 100 * R2adj, cross_validated = 100 * (1 - sum(e1^2) / sum(e0^2)))
}

variance_df <- dplyr::tibble()
if (nrow(sig_markers)) {
  pve_list <- list(); var_list <- list(); clump_list <- list()
  for (v in unique(sig_markers$variable)) {
    mod <- model_for(v); nf <- null_fit(v, mod)
    mk <- sig_markers$marker[sig_markers$variable == v]
    pve_list[[v]] <- data.frame(variable = v, marker_pve(nf, mk))
    # independent loci, then a joint GLS fit of one lead marker per locus (strongest loci first)
    dv <- as.data.frame(sig_markers[sig_markers$variable == v, c("marker", "chrom", "pos", "p_value")])
    dv$chrom <- as.character(dv$chrom)
    cl <- clump_loci(dv, geno_mat_final[nf$ok, , drop = FALSE] + 1, clump_window, clump_r2)
    lead <- cl$marker[cl$lead_marker]
    kmax <- max(1, floor(max_lead_frac * sum(nf$ok))); if (length(lead) > kmax) lead <- head(lead, kmax)
    Gs <- nf$rot(geno_mat_final[nf$ok, lead, drop = FALSE] + 1)
    jr <- joint_r2(as.numeric(nf$ys), nf$Xs, Gs, folds = cv_folds)
    clump_list[[v]] <- data.frame(variable = v, cl[, c("marker", "locus", "lead_marker")])
    vy <- var(nf$y_raw)
    var_list[[v]] <- dplyr::tibble(variable = v, model = mod, n_sig_markers = length(mk),
      n_independent_loci = max(cl$locus), n_lead_markers_fitted = length(lead),
      structure_PCs_percent = 100 * var(as.numeric(nf$X %*% nf$fit$beta)) / vy,
      covariates_percent = 100 * (1 - var(nf$y_an) / vy),
      polygenic_percent = 100 * nf$fit$Vu * nf$Kd / vy,
      lead_markers_R2_adjusted_percent = jr[["adjusted"]],               # % of residual variance after structure + kinship
      lead_markers_R2_crossvalidated_percent = jr[["cross_validated"]],   # out-of-sample; can be < 0 if markers do not predict
      lead_markers_R2_in_sample_percent = jr[["in_sample"]],
      sum_single_marker_PVE_percent = sum(pve_list[[v]]$PVE_percent))
  }
  sig_markers <- sig_markers %>% dplyr::left_join(dplyr::bind_rows(pve_list), by = c("variable", "marker")) %>%
    dplyr::left_join(dplyr::bind_rows(clump_list), by = c("variable", "marker")) %>%
    dplyr::mutate(effect_direction = ifelse(beta > 0, "ALT increases variable", "ALT decreases variable"))
  variance_df <- dplyr::bind_rows(var_list)
  print(variance_df %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2))), n = Inf)
}

# ==============================================================================
# 9  MANHATTAN PLOTS
# ==============================================================================
chrom_order <- function(ch) {                       # natural order: 1,2,...,10; 1A,1B...
  u <- unique(ch); num <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", u)))
  u[order(is.na(num), num, u)]
}
mf <- setNames(vapply(env_vars, model_for, ""), env_vars)
man <- gwas_analyzed %>% dplyr::filter(model == mf[variable], !is.na(p_value)) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order(chrom)))
offs <- man %>% dplyr::distinct(chrom, pos) %>% dplyr::group_by(chrom) %>% dplyr::summarise(len = max(pos), .groups = "drop") %>%
  dplyr::mutate(offset = cumsum(as.numeric(dplyr::lag(len, default = 0))), mid = offset + len / 2)
man <- man %>% dplyr::left_join(offs, by = "chrom") %>% dplyr::mutate(bp_cum = pos + offset,
  col = factor(as.integer(chrom) %% 2), sig = paste(variable, marker) %in% paste(sig_markers$variable, sig_markers$marker))
thr <- man %>% dplyr::group_by(variable) %>% dplyr::summarise(bonf = -log10(0.05 / dplyr::n()), .groups = "drop")
if (nrow(sig_markers)) thr <- thr %>% dplyr::left_join(sig_markers %>% dplyr::group_by(variable) %>%
                                                  dplyr::summarise(fdr_line = min(score), .groups = "drop"), by = "variable")
p_man <- ggplot(man, aes(bp_cum, score)) +
  geom_point(aes(colour = col), size = 0.5, show.legend = FALSE) +
  geom_point(data = dplyr::filter(man, sig), colour = "#B2182B", size = 1.2) +
  geom_hline(data = thr, aes(yintercept = bonf), linetype = 2, colour = "grey30") +
  { if ("fdr_line" %in% names(thr)) geom_hline(data = dplyr::filter(thr, !is.na(fdr_line)), aes(yintercept = fdr_line),
                                               linetype = 3, colour = "#B2182B") } +
  scale_colour_manual(values = c("0" = "blue4", "1" = "orange3")) +
  scale_x_continuous(breaks = offs$mid, labels = offs$chrom, expand = c(0.01, 0)) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  labs(x = "Chromosome", y = expression(-log[10](italic(p))),
       title = paste0("Manhattan plots (", if (sig_model == "best") "model with \u03bb closest to 1" else model_labels[sig_model], ")"),
       subtitle = paste0("Dashed = Bonferroni 5%; dotted red = FDR ", fdr_threshold, "; red points = significant")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 7, angle = if (nrow(offs) > 15) 90 else 0))
save_fig("egwas_manhattan", p_man, width = 12, height = 1.9 * n_vars + 1.5)

# ==============================================================================
# 10  SAVE RESULTS
# ==============================================================================
save(gwas_combined, gwas_analyzed, p_lambda, sig_markers, variance_df, egwas_models, K,
     snp_info_final, geno_mat_final, env_data, file = file.path(output_dir, "egwas_results.RData"))
if (nrow(sig_markers)) {
  out_cols <- c("marker", "chrom", "pos", "variable", "model", "score", "p_value", "q_value", "fdr_method",
                "alt_allele_freq", "beta", "beta_se", "effect_direction", "PVE_percent", "partial_R2_percent", "R2_raw_percent",
                "locus", "lead_marker")
  write.csv(sig_markers[, intersect(out_cols, names(sig_markers))], file.path(output_dir, "significant_markers.csv"), row.names = FALSE)
  write.csv(variance_df, file.path(output_dir, "variance_explained.csv"), row.names = FALSE)
}
write.csv(p_lambda, file.path(output_dir, "genomic_inflation_factors.csv"), row.names = FALSE)
if (!nrow(sig_markers)) { msg("No significant markers at FDR", fdr_threshold, "- visualizations of hits skipped"); quit(save = "no") }

# ==============================================================================
# 11  VISUALIZATION OF SIGNIFICANT MARKERS
# ==============================================================================
msg("Visualizations")
sig_per_variable <- sig_markers %>% dplyr::group_by(variable) %>%
  dplyr::summarise(n_markers = dplyr::n_distinct(marker), n_associations = dplyr::n(), max_score = max(score),
            max_PVE_percent = max(PVE_percent, na.rm = TRUE), .groups = "drop") %>% dplyr::arrange(desc(n_markers))
sig_per_marker <- sig_markers %>% dplyr::group_by(marker, chrom, pos) %>%
  dplyr::summarise(n_variables = dplyr::n_distinct(variable), variables = paste(variable, collapse = ", "),
            max_score = max(score), max_PVE_percent = max(PVE_percent, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(n_variables), desc(max_score))
markers_per_chrom <- sig_markers %>% dplyr::distinct(marker, chrom) %>% dplyr::count(chrom, name = "n_markers") %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order(chrom))) %>% dplyr::arrange(chrom)
write.csv(sig_per_variable, file.path(figures_dir, "table_summary_by_variable.csv"), row.names = FALSE)
write.csv(sig_per_marker, file.path(figures_dir, "table_summary_by_marker.csv"), row.names = FALSE)
write.csv(markers_per_chrom, file.path(figures_dir, "table_chromosome_distribution.csv"), row.names = FALSE)

p1 <- ggplot(sig_per_variable, aes(reorder(variable, n_markers), n_markers)) + geom_col(fill = "steelblue") +
  geom_text(aes(label = n_markers), hjust = -0.3, size = 3) + coord_flip() +
  labs(title = "Significant markers per variable", subtitle = paste("FDR", fdr_threshold), x = NULL, y = "Markers")
save_fig("01_markers_per_variable", p1, 9, max(5, 0.3 * nrow(sig_per_variable) + 2))
p2 <- ggplot(head(sig_per_marker, 20), aes(reorder(marker, n_variables), n_variables)) + geom_col(fill = "coral") +
  geom_text(aes(label = n_variables), hjust = -0.3, size = 3) + coord_flip() +
  labs(title = "Top 20 markers associated with multiple variables", x = NULL, y = "Variables") +
  theme(axis.text.y = element_text(size = 7, family = "mono"))
save_fig("02_pleiotropic_markers", p2, 9, 7)
p3 <- ggplot(sig_markers, aes(score)) + geom_histogram(bins = 30, fill = "steelblue", colour = "black") +
  geom_vline(xintercept = median(sig_markers$score), linetype = 2, colour = "red") +
  labs(title = "Association scores of significant markers", subtitle = "Red = median", x = expression(-log[10](italic(p))), y = "Count")
save_fig("03_score_distribution", p3, 7, 5)
p4 <- ggplot(sig_markers, aes(reorder(variable, PVE_percent, median), PVE_percent)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA) + geom_jitter(width = 0.2, size = 1, alpha = 0.6) + coord_flip() +
  labs(title = "Percent variance explained by significant markers", x = NULL, y = "PVE (%) per marker")
save_fig("04_PVE_by_variable", p4, 9, max(5, 0.3 * nrow(sig_per_variable) + 2))
p5 <- ggplot(markers_per_chrom, aes(chrom, n_markers)) + geom_col(fill = "darkgreen") +
  geom_text(aes(label = n_markers), vjust = -0.4, size = 3) + labs(title = "Significant markers per chromosome", x = "Chromosome", y = "Markers")
save_fig("05_markers_per_chromosome", p5, 9, 5)

env_type <- function(v) {
  n <- suppressWarnings(as.integer(sub("^bio[_ ]?0*([0-9]+)$", "\\1", tolower(v))))
  ifelse(!is.na(n) & n <= 11, "Temperature", ifelse(!is.na(n) & n <= 19, "Precipitation", "Soil/Other"))
}
d6 <- sig_markers %>% dplyr::mutate(env_type = env_type(variable), chrom = factor(chrom, levels = chrom_order(chrom)))
p6 <- ggplot(d6, aes(pos / 1e6, chrom, colour = env_type, size = PVE_percent)) + geom_point(alpha = 0.7) +
  scale_colour_manual(values = c(Temperature = "#d62728", Precipitation = "#1f77b4", "Soil/Other" = "#2ca02c"), name = "Variable type") +
  scale_size_continuous(range = c(1, 5), name = "PVE (%)") +
  labs(title = "Genomic positions of significant markers", x = "Position (Mb)", y = "Chromosome")
save_fig("06_marker_positions", p6, 12, 7)

top_markers <- head(sig_per_marker$marker, 30)
p7 <- ggplot(dplyr::filter(sig_markers, marker %in% top_markers), aes(variable, marker, fill = PVE_percent)) +
  geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.1f", PVE_percent)), size = 2.3) +
  scale_fill_viridis_c(option = "plasma", name = "PVE (%)") +
  labs(title = "Percent variance explained: top 30 markers", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7, family = "mono"), panel.grid = element_blank())
save_fig("07_PVE_heatmap", p7, 11, 9)

# ---- allele-environment boxplots ----
top_m <- head(sig_per_marker$marker, n_top_markers); top_v <- head(sig_per_variable$variable, 10)
bx <- sig_markers %>% dplyr::filter(marker %in% top_m, variable %in% top_v) %>%
  dplyr::select(marker, variable, PVE_percent) %>%
  dplyr::left_join(map_dfr(top_m, function(m) data.frame(marker = m, sample_id = rownames(geno_mat_final), genotype = geno_mat_final[, m])),
            by = "marker", relationship = "many-to-many") %>%
  dplyr::left_join(env_data %>% dplyr::select(sample_id, all_of(top_v)) %>% pivot_longer(-sample_id, names_to = "variable", values_to = "env_value"),
            by = c("sample_id", "variable")) %>%
  dplyr::mutate(genotype_class = genotype_class(genotype))
if (nrow(bx)) {
  lab9 <- bx %>% dplyr::distinct(marker, variable, PVE_percent) %>% dplyr::mutate(label = sprintf("PVE %.1f%%", PVE_percent))
  p9 <- ggplot(bx, aes(genotype_class, env_value, fill = genotype_class)) + geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.4) +
    geom_text(data = lab9, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 2.3, fontface = "bold") +
    facet_grid(variable ~ marker, scales = "free_y") + scale_fill_manual(values = geno_cols, name = "Genotype") +
    labs(title = "Environmental values by genotype at significant markers", x = NULL, y = "Environmental value",
         subtitle = "PVE = % of the variable's variance explained by the marker (mixed-model estimate)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), strip.text.x = element_text(size = 6, angle = 90),
          strip.text.y = element_text(size = 7), legend.position = "bottom")
  save_fig("09_allele_environment_boxplots", p9, 2 + 1.1 * dplyr::n_distinct(bx$marker), 2 + 1.5 * dplyr::n_distinct(bx$variable))
  write.csv(bx %>% dplyr::group_by(marker, variable, genotype_class) %>%
              dplyr::summarise(n = dplyr::n(), mean_env = mean(env_value), sd_env = sd(env_value), PVE_percent = dplyr::first(PVE_percent), .groups = "drop"),
            file.path(figures_dir, "table_allele_environment_stats.csv"), row.names = FALSE)
}

# ---- geographic allele distribution on a world map ----
if (HAS_COORDS) {
  map_m <- head(sig_per_marker$marker, n_map_markers)
  md <- map_dfr(map_m, function(m) data.frame(marker = m, sample_id = rownames(geno_mat_final), genotype = geno_mat_final[, m])) %>%
    dplyr::left_join(env_data %>% dplyr::select(sample_id, latitude, longitude), by = "sample_id") %>%
    dplyr::filter(!is.na(latitude), !is.na(longitude)) %>%
    dplyr::mutate(genotype_class = genotype_class(genotype),
           marker_lab = paste0(marker, "\n", sig_per_marker$variables[match(marker, sig_per_marker$marker)]))
  world <- ggplot2::map_data("world")
  lon_r <- range(md$longitude); lat_r <- range(md$latitude)
  pad_lon <- max(3, diff(lon_r) * 0.1); pad_lat <- max(3, diff(lat_r) * 0.1)
  xl <- c(max(-180, lon_r[1] - pad_lon), min(180, lon_r[2] + pad_lon))
  yl <- c(max(-90, lat_r[1] - pad_lat), min(90, lat_r[2] + pad_lat))
  # continent labels: first candidate anchor (on land) that lies well inside the map extent
  cand <- data.frame(label = rep(c("AFRICA", "EUROPE", "ASIA", "NORTH\nAMERICA", "SOUTH\nAMERICA", "OCEANIA"), c(4, 2, 5, 2, 2, 2)),
                     lon = c(20, 30, 25, 35, 15, 30, 90, 100, 78, 105, 60, -100, -95, -60, -55, 135, 145),
                     lat = c(5, -5, 15, -15, 50, 55, 45, 30, 22, 15, 45, 45, 35, -15, -5, -25, -30))
  mx <- 0.08 * diff(xl); my <- 0.08 * diff(yl)
  cand <- cand[cand$lon > xl[1] + mx & cand$lon < xl[2] - mx & cand$lat > yl[1] + my & cand$lat < yl[2] - my, ]
  continents <- cand[!duplicated(cand$label), ]
  p_map <- ggplot() +
    geom_polygon(data = world, aes(long, lat, group = group), fill = "#F2EFE6", colour = "grey55", linewidth = 0.15) +
    geom_text(data = continents, aes(lon, lat, label = label), colour = "grey45", size = 3.2, fontface = "bold", alpha = 0.6, lineheight = 0.85) +
    geom_point(data = md, aes(longitude, latitude, colour = genotype_class, shape = genotype_class), size = 1.8, alpha = 0.85) +
    facet_wrap(~ marker_lab, ncol = min(3, length(map_m))) +
    scale_colour_manual(values = geno_cols, name = "Genotype") + scale_shape_manual(values = geno_shape, name = "Genotype") +
    coord_quickmap(xlim = xl, ylim = yl, expand = FALSE) +
    labs(title = "Geographic distribution of genotypes at the top pleiotropic markers",
         subtitle = "Associated variables listed under each marker", x = "Longitude", y = "Latitude") +
    theme(panel.background = element_rect(fill = "#D6EAF8"), panel.grid = element_line(colour = "white", linewidth = 0.2),
          legend.position = "bottom", strip.text = element_text(size = 8, face = "bold"))
  asp <- diff(yl) / diff(xl); nc <- min(3, length(map_m)); nr <- ceiling(length(map_m) / nc)
  save_fig("11_geographic_allele_distribution_map", p_map, width = 4.5 * nc + 1, height = 4.5 * nc * asp / nc * nr + 2)
  write.csv(md %>% dplyr::group_by(marker, genotype_class) %>%
              dplyr::summarise(n = dplyr::n(), mean_lat = mean(latitude), sd_lat = sd(latitude), mean_lon = mean(longitude), sd_lon = sd(longitude), .groups = "drop"),
            file.path(figures_dir, "table_geographic_allele_stats.csv"), row.names = FALSE)
} else msg("No latitude/longitude columns: geographic map skipped")

summary_plot <- (p1 | p5) / (p3 | p4)
save_fig("10_combined_summary", summary_plot, 15, 11)

# ---- HTML report (uses the PNG copies) ----
img <- function(f, cap) if (file.exists(file.path(figures_dir, f)))
  sprintf('<h3>%s</h3><img src="%s"><p><a href="%s">PDF version</a></p>', cap, f, sub("png$", "pdf", f)) else ""
rows <- paste(sprintf("<tr><td>%s</td><td>%d</td><td>%d</td><td>%.2f</td><td>%.1f</td></tr>", sig_per_variable$variable,
                      sig_per_variable$n_markers, sig_per_variable$n_associations, sig_per_variable$max_score, sig_per_variable$max_PVE_percent), collapse = "\n")
html <- paste0('<!DOCTYPE html><html><head><meta charset="utf-8"><title>Environmental GWAS summary</title><style>
body{font-family:Arial,sans-serif;margin:40px;background:#f5f5f5}h1{color:#2c3e50;border-bottom:3px solid #3498db}
.box{background:#fff;padding:20px;margin:20px 0;border-radius:8px;box-shadow:0 2px 4px rgba(0,0,0,.1)}
img{max-width:100%;margin:10px 0}table{border-collapse:collapse;width:100%}th,td{border:1px solid #ddd;padding:8px}
th{background:#3498db;color:#fff}</style></head><body><h1>Environmental GWAS summary</h1>',
  sprintf('<div class="box"><p><b>%d</b> associations, <b>%d</b> markers, <b>%d</b> variables, <b>%d</b> chromosomes.
  FDR %.2f (%s); significance model: %s.</p></div>', nrow(sig_markers), dplyr::n_distinct(sig_markers$marker),
  dplyr::n_distinct(sig_markers$variable), dplyr::n_distinct(sig_markers$chrom), fdr_threshold, sig_markers$fdr_method[1],
  if (sig_model == "best") "per-variable model with lambda closest to 1" else model_labels[sig_model]),
  '<div class="box"><h2>Variables</h2><table><tr><th>Variable</th><th>Markers</th><th>Associations</th><th>Max -log10(p)</th><th>Max PVE (%)</th></tr>',
  rows, '</table></div><div class="box">',
  img("egwas_qqplots_all_models.png", "QQ plots (all models)"), img("egwas_manhattan.png", "Manhattan plots"),
  img("01_markers_per_variable.png", "Markers per variable"), img("04_PVE_by_variable.png", "Percent variance explained"),
  img("07_PVE_heatmap.png", "PVE heatmap"), img("09_allele_environment_boxplots.png", "Allele-environment boxplots"),
  img("11_geographic_allele_distribution_map.png", "Geographic allele distribution"),
  img("10_combined_summary.png", "Combined summary"),
  sprintf('</div><p style="color:#7f8c8d;font-size:12px">Generated %s</p></body></html>', format(Sys.time())))
writeLines(html, file.path(figures_dir, "summary_report.html"))

msg("Done. Tables:", normalizePath(output_dir), "| figures and report:", normalizePath(figures_dir))
