################################################################################
# finger_millet_genomic_offset_revised.R
#
# Genomic (genetic) offset for finger millet, using the same method as the yam
# script 11_genomic_offset_revised.R. Replaces finger_millet_Genetic_Offset_Analysis.Rmd.
#
# Why the Rmd maps looked odd, and what this script does instead (matching the yam fixes):
#  1. GLOBAL projection on a coarse 1.8-degree grid: offsets were computed for the whole
#     world (Sahara, Arctic, Amazon ...), climates finger millet was never sampled in, and
#     the colour scale was set by those extremes; the "regional" map was a crop of that
#     coarse grid (blocky). -> Here the map covers the sampled range + buffer at native
#     WorldClim resolution, like the yam script.
#  2. No extrapolation control. -> Cells whose CURRENT climate lies outside the climate of
#     the sampled accessions (ExDet: NT1 < 0 or NT2 > 1; Mesgaran et al. 2014) are masked
#     (grey); the share of cells with non-analog FUTURE climate is reported, not masked.
#  3. All 19 bioclim variables entered LFMM2 together. They are strongly collinear, which
#     makes the LFMM2 effect sizes unstable and the offset dominated by arbitrary variable
#     combinations. -> Variables are thinned by VIF (< CFG$vif_max), as in the yam analysis.
#  4. LFMM2 used only SNPs with NO missing call: 1,497 of 8,778 SNPs (17%) in
#     finger_millet_404.vcf. -> SNPs are filtered on missingness/MAF and the remaining
#     missing calls imputed, so ~all SNPs are used.
#  5. The offset used only candidate SNPs (p < 1e-5). -> LEA::genetic.gap on all SNPs, as in
#     the yam script (CFG$candidate_p can restrict it to candidates if wanted).
#  6. The main map was a GRAND MEAN over SSP2-4.5 + SSP5-8.5 and 2041-2060 + 2081-2100.
#     -> One map per SSP x period (each the mean over GCMs); scenarios are never averaged.
#  7. Rainbow palette. -> Single sequential palette (viridis), low -> high in the legend.
#
# Offset at a location = genetic gap between current and future climate AT THAT SAME
# location (LEA::genetic.gap; Gain et al. 2023), averaged over GCMs, with the SD among
# GCMs mapped as uncertainty.
#
# Inputs : VCF; coordinate CSV (SampleID, Latitude, Longitude, [Population]).
# Climate: WorldClim 2.1 (current) and CMIP6 (future) via geodata, cached in CFG$clim_dir.
# Outputs (CFG$out_dir): offset maps (PDF + PNG) per scenario incl. regional panels,
#   GeoTIFFs (mean, SD among GCMs, extrapolation mask), per-accession offsets, group
#   summaries, extrapolation summary, variables used.
################################################################################

# ==============================================================================
# 0  CONFIGURATION
# ==============================================================================
CFG <- list(
  vcf          = "finger_millet_404.vcf",
  coords_csv   = "millet_passport_test.csv",     # SampleID, Latitude, Longitude, [Population or Region]
  clim_dir     = "cmip6_cache",                   # WorldClim/CMIP6 download cache
  clim_res     = 5,                               # arc-minutes: 2.5, 5 or 10 (yam used CFG$clim_res)
  gcms         = c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CNRM-CM6-1", "IPSL-CM6A-LR", "MIROC6"),
  ssps         = c("585"),                        # one map per SSP x period; e.g. c("245", "585")
  periods      = c("2081-2100"),                  #   e.g. c("2041-2060", "2081-2100")
  map_buffer   = 3,                               # degrees around the sampled range
  K            = 5,                               # LFMM2 latent factors
  vif_max      = 10,                              # stepwise VIF threshold for bioclim variables
  vars         = NULL,                            # NULL = VIF selection; or e.g. c("bio1","bio4","bio12","bio15")
  maf          = 0.05,
  max_snp_missing = 0.20,                         # SNPs above this missingness are dropped; the rest imputed
  candidate_p  = NA,                              # NA = all SNPs (as yam); e.g. 1e-5 = LFMM2 candidates only
  # regional panels cut from the same map (lon_min, lon_max, lat_min, lat_max); NULL = none
  regions      = list(Africa = c(5, 52, -36, 20), `South Asia` = c(65, 92, 4, 32)),
  common_scale = TRUE,                            # same colour limits across scenario maps
  out_dir      = "genomic_offset_results",
  seed         = 2026
)

# ==============================================================================
# 1  PACKAGES AND HELPERS
# ==============================================================================
need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
}
need(c("terra", "geodata", "ggplot2", "dplyr", "maps", "patchwork"))
if (!requireNamespace("LEA", quietly = TRUE)) {
  need("BiocManager"); BiocManager::install("LEA", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({ library(terra); library(ggplot2); library(dplyr) })
set.seed(CFG$seed)
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$clim_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(CFG$out_dir, f)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
save_fig <- function(name, p, width, height) {
  ggsave(out(paste0(name, ".pdf")), p, width = width, height = height,
         device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf)
  ggsave(out(paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
}
bio_names <- function(n) {                         # "wc2.1_5m_bio_1" / "bio01" -> "bio1"
  d <- suppressWarnings(as.integer(sub("^.*?([0-9]+)$", "\\1", tolower(n), perl = TRUE)))
  ifelse(is.na(d), n, paste0("bio", d))
}
# climate sources (geodata); kept as functions so another source can be swapped in
get_current <- function() {
  r <- geodata::worldclim_global("bio", res = CFG$clim_res, path = CFG$clim_dir)
  names(r) <- bio_names(names(r)); r
}
get_future <- function(gcm, ssp, period) {
  r <- geodata::cmip6_world(gcm, ssp, period, var = "bioc", res = CFG$clim_res, path = CFG$clim_dir)
  names(r) <- bio_names(names(r)); r
}
# ExDet (Mesgaran et al. 2014): NT1 = univariate novelty (< 0 outside the range of ref),
# NT2 = Mahalanobis distance relative to the most extreme reference point (> 1 = novel combination)
exdet <- function(ref, new) {
  mn <- apply(ref, 2, min); mx <- apply(ref, 2, max); rg <- mx - mn
  lo <- sweep(sweep(new, 2, mn, "-"), 2, rg, "/"); hi <- sweep(sweep(-new, 2, -mx, "-"), 2, rg, "/")
  nt1 <- rowSums(pmin(lo, hi, 0))
  mu <- colMeans(ref); S <- cov(ref)
  Si <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
  d_ref <- max(mahalanobis(ref, mu, Si, inverted = TRUE))
  nt2 <- mahalanobis(new, mu, Si, inverted = TRUE) / d_ref
  list(nt1 = nt1, nt2 = nt2, analog = nt1 >= 0 & nt2 <= 1)
}
vif_select <- function(E, thr) {                   # drop the highest-VIF variable until all < thr
  v <- colnames(E)
  repeat {
    if (length(v) < 3) break
    vif <- sapply(v, function(x) { r2 <- summary(lm(E[, x] ~ E[, setdiff(v, x)]))$r.squared; 1 / (1 - r2) })
    if (max(vif) < thr) break
    v <- setdiff(v, names(which.max(vif)))
  }
  v
}

# ==============================================================================
# 2  GENOTYPES AND SAMPLES
# ==============================================================================
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
  G <- t(d); rownames(G) <- trimws(smp); colnames(G) <- paste(v[[1]], v[[2]], sep = "_"); G
}
msg("Reading", CFG$vcf)
G <- read_vcf_dosage(CFG$vcf)
co <- read.csv(CFG$coords_csv, stringsAsFactors = FALSE, check.names = FALSE)
nm <- tolower(names(co))
meta <- data.frame(SampleID = trimws(co[[1]]),
                   Latitude = as.numeric(co[[grep("^lat", nm)[1]]]), Longitude = as.numeric(co[[grep("^lon", nm)[1]]]),
                   Population = if (length(grep("^pop|^region|^group", nm))) as.character(co[[grep("^pop|^region|^group", nm)[1]]]) else NA)
meta <- meta[!duplicated(meta$SampleID) & !is.na(meta$Latitude) & !is.na(meta$Longitude), ]
ids <- intersect(rownames(G), meta$SampleID)
msg("Samples in VCF / coordinates / both:", nrow(G), "/", nrow(meta), "/", length(ids))
G <- G[ids, , drop = FALSE]; meta <- meta[match(ids, meta$SampleID), ]

# SNP filters and imputation (LFMM2 / genetic.gap need complete genotypes)
p <- colMeans(G, na.rm = TRUE) / 2; miss <- colMeans(is.na(G))
keep <- !is.na(p) & miss <= CFG$max_snp_missing & pmin(p, 1 - p) >= CFG$maf
G <- G[, keep, drop = FALSE]
n_imp <- sum(is.na(G))
for (j in which(colSums(is.na(G)) > 0)) { x <- G[, j]; x[is.na(x)] <- round(mean(x, na.rm = TRUE)); G[, j] <- x }
msg("SNPs used:", ncol(G), "(missingness <=", CFG$max_snp_missing, ", MAF >=", CFG$maf, ");",
    n_imp, "missing calls imputed (", round(100 * n_imp / length(G), 2), "% )")

# ==============================================================================
# 3  CURRENT CLIMATE, VARIABLES, EXTENT
# ==============================================================================
ext0 <- ext(min(meta$Longitude) - CFG$map_buffer, max(meta$Longitude) + CFG$map_buffer,
            min(meta$Latitude) - CFG$map_buffer, max(meta$Latitude) + CFG$map_buffer)
EV <- as.vector(ext0)
now_all <- crop(get_current(), ext0)
pts <- vect(meta, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
E_all <- as.matrix(terra::extract(now_all, pts)[, names(now_all)])
no_clim <- !complete.cases(E_all)
if (any(no_clim)) {                                 # coastal sites in a sea cell: use the nearest land cells
  filled <- terra::focal(now_all, w = 5, fun = "mean", na.policy = "only", na.rm = TRUE)
  E_fill <- as.matrix(terra::extract(filled, pts[which(no_clim)])[, names(now_all)])
  E_all[no_clim, ] <- E_fill
  msg(sum(no_clim), "sample(s) fell in a sea/lake cell; climate taken from adjacent land cells (<= 2 cells away) for",
      sum(complete.cases(E_fill)))
  no_clim <- !complete.cases(E_all)
}
if (any(no_clim)) {
  msg("WARNING:", sum(no_clim), "sample(s) have no climate data nearby (check coordinates) and are left out:",
      paste(head(meta$SampleID[no_clim], 10), collapse = ", "))
  G <- G[!no_clim, , drop = FALSE]; meta <- meta[!no_clim, ]; E_all <- E_all[!no_clim, , drop = FALSE]
  pts <- vect(meta, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
}
E_all <- E_all[, apply(E_all, 2, sd) > 0, drop = FALSE]
vars <- if (is.null(CFG$vars)) vif_select(E_all, CFG$vif_max) else intersect(CFG$vars, colnames(E_all))
vars <- vars[order(as.integer(sub("bio", "", vars)))]
msg("Offset predictors (VIF <", CFG$vif_max, "):", paste(vars, collapse = ", "))
now <- now_all[[vars]]
E_train <- E_all[, vars, drop = FALSE]; rownames(E_train) <- meta$SampleID
cells <- which(complete.cases(values(now)))         # land cells only - no imputation
E_now_grid <- values(now)[cells, , drop = FALSE]
novel <- !exdet(E_train, E_now_grid)$analog         # mask: CURRENT climate outside the sampled envelope
msg("Map:", paste(round(EV, 1), collapse = " / "), "| land cells:", length(cells),
    "| masked (current climate non-analog):", round(100 * mean(novel), 1), "%")
write.csv(data.frame(variable = vars, min_at_samples = apply(E_train, 2, min), max_at_samples = apply(E_train, 2, max)),
          out("offset_variables_used.csv"), row.names = FALSE)

# ---- LFMM2 fitted ONCE; offset computed exactly as LEA::genetic.gap ----
# genetic.gap(scale = TRUE) scales env with the training mean/SD, fits lfmm2, and returns
# offset = rowSums(((X_now - X_fut) %*% t(B))^2) / L. That equals rowSums((D %*% C) * D) with
# D = X_now - X_fut and C = t(B) %*% B / L, which avoids a locations x SNPs matrix (GBs for a map)
# and refitting LFMM2 for every call. Checked against LEA::genetic.gap below.
mx <- colMeans(E_train); sx <- apply(E_train, 2, sd)
sc <- function(E) sweep(sweep(E, 2, mx, "-"), 2, sx, "/")
lf <- LEA::lfmm2(input = G, env = sc(E_train), K = CFG$K, effect.sizes = TRUE)
cand <- NULL
if (!is.na(CFG$candidate_p)) {
  pv <- LEA::lfmm2.test(lf, input = G, env = sc(E_train), full = TRUE)$pvalues
  cand <- which(pv < CFG$candidate_p)
  msg("Candidate loci (p <", CFG$candidate_p, "):", length(cand))
  if (length(cand) < 2) stop("Fewer than 2 candidate loci; use candidate_p = NA (all SNPs)")
}
B <- as.matrix(lf@B); if (!is.null(cand)) B <- B[cand, , drop = FALSE]
Cb <- crossprod(B) / nrow(B)
ggap <- function(new_env, pred_env) { D <- sc(new_env) - sc(pred_env); rowSums((D %*% Cb) * D) }
checked <- FALSE

# ==============================================================================
# 4  FUTURE CLIMATE AND OFFSET, one SSP x period at a time
# ==============================================================================
world <- ggplot2::map_data("world")
pretty_ssp <- function(s) sprintf("SSP%s-%s.%s", substr(s, 1, 1), substr(s, 2, 2), substr(s, 3, 3))
scen <- expand.grid(ssp = CFG$ssps, period = CFG$periods, stringsAsFactors = FALSE)
results <- list(); extrap <- list(); acc_all <- list()
for (k in seq_len(nrow(scen))) {
  ssp <- scen$ssp[k]; period <- scen$period[k]; sid <- paste0("ssp", ssp, "_", period)
  msg("==", pretty_ssp(ssp), period, "==")
  grid_off <- acc_off <- list(); fut_novel <- rep(FALSE, length(cells))
  for (gcm in CFG$gcms) {
    fut <- tryCatch(crop(get_future(gcm, ssp, period), ext0)[[vars]],
                    error = function(e) { msg("  ", gcm, "failed:", conditionMessage(e)); NULL })
    if (is.null(fut)) next
    E_fut_grid <- values(fut)[cells, , drop = FALSE]; ok <- complete.cases(E_fut_grid)
    E_fut_acc <- as.matrix(terra::extract(fut, pts)[, vars])
    na_acc <- !complete.cases(E_fut_acc)              # coastal sites: same neighbour fill as current climate
    if (any(na_acc)) E_fut_acc[na_acc, ] <- as.matrix(terra::extract(
      terra::focal(fut, w = 5, fun = "mean", na.policy = "only", na.rm = TRUE), pts[which(na_acc)])[, vars])
    go <- rep(NA_real_, length(cells))
    go[ok] <- ggap(E_now_grid[ok, , drop = FALSE], E_fut_grid[ok, , drop = FALSE])
    grid_off[[gcm]] <- go
    a <- rep(NA_real_, nrow(E_train)); oka <- complete.cases(E_fut_acc)
    a[oka] <- ggap(E_train[oka, , drop = FALSE], E_fut_acc[oka, , drop = FALSE]); acc_off[[gcm]] <- a
    if (!checked) {                                 # identical to LEA::genetic.gap (as used in the yam script)?
      ref <- LEA::genetic.gap(input = G, env = E_train, new.env = E_train[oka, , drop = FALSE],
                              pred.env = E_fut_acc[oka, , drop = FALSE], K = CFG$K, scale = TRUE, candidate.loci = cand)$offset
      rel <- max(abs(ref - a[oka])) / max(abs(ref))
      msg("  check vs LEA::genetic.gap at the accessions: max relative difference", signif(rel, 3))
      if (rel > 1e-6) warning("Offset differs from LEA::genetic.gap (relative difference ", signif(rel, 3), ")")
      checked <- TRUE
    }
    xd <- exdet(E_train, E_fut_grid[ok, , drop = FALSE]); fn <- rep(TRUE, length(cells)); fn[ok] <- !xd$analog
    fut_novel <- fut_novel | fn
    msg("  ", gcm, ": median offset over unmasked land", signif(median(go[!novel], na.rm = TRUE), 3))
  }
  if (!length(grid_off)) { msg("  no GCM succeeded - scenario skipped"); next }
  GM <- do.call(cbind, grid_off)
  off_mean <- rowMeans(GM, na.rm = TRUE); off_sd <- if (ncol(GM) > 1) apply(GM, 1, sd, na.rm = TRUE) else rep(NA, nrow(GM))
  r_mean <- r_sd <- r_mask <- now[[1]]; values(r_mean) <- values(r_sd) <- values(r_mask) <- NA
  r_mean[cells] <- ifelse(novel, NA, off_mean); r_sd[cells] <- ifelse(novel, NA, off_sd); r_mask[cells] <- as.numeric(novel)
  rs <- c(r_mean, r_sd, r_mask); names(rs) <- c("offset_mean", "offset_sd_among_GCMs", "extrapolation_masked")
  writeRaster(rs, out(paste0("offset_", sid, ".tif")), overwrite = TRUE)
  results[[sid]] <- list(ssp = ssp, period = period, r = rs, n_gcm = length(grid_off))
  extrap[[sid]] <- data.frame(scenario = sid, gcms = length(grid_off), land_cells = length(cells),
                              current_nonanalog_pct = 100 * mean(novel), future_nonanalog_pct = 100 * mean(fut_novel),
                              grid_min = min(off_mean[!novel], na.rm = TRUE), grid_median = median(off_mean[!novel], na.rm = TRUE),
                              grid_max = max(off_mean[!novel], na.rm = TRUE), mean_sd_among_GCMs = mean(off_sd[!novel], na.rm = TRUE))
  AM <- do.call(cbind, acc_off)
  acc_all[[sid]] <- data.frame(SampleID = meta$SampleID, Population = meta$Population, Latitude = meta$Latitude,
                               Longitude = meta$Longitude, scenario = sid, offset_mean = rowMeans(AM, na.rm = TRUE),
                               offset_sd = if (ncol(AM) > 1) apply(AM, 1, sd, na.rm = TRUE) else NA)
  msg("  Land cells with non-analog FUTURE climate (any GCM):", round(100 * mean(fut_novel), 1), "% (reported, not masked)")
}
if (!length(results)) stop("No scenario produced an offset map")
write.csv(bind_rows(extrap), out("offset_extrapolation_and_grid_summary.csv"), row.names = FALSE)
print(bind_rows(extrap))
acc <- bind_rows(acc_all)
write.csv(acc, out("offset_per_accession.csv"), row.names = FALSE)

# ==============================================================================
# 5  MAPS
# ==============================================================================
lims <- if (CFG$common_scale) range(unlist(lapply(results, function(x) values(x$r[["offset_mean"]]))), na.rm = TRUE) else NULL
offset_map <- function(r, xlim, ylim, title, fill_name, limits = NULL, pt_size = 1.4) {
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE); names(df)[3:4] <- c("value", "masked")
  df <- df[!is.na(df$masked) & df$x >= xlim[1] & df$x <= xlim[2] & df$y >= ylim[1] & df$y <= ylim[2], ]
  sm <- meta[meta$Longitude >= xlim[1] & meta$Longitude <= xlim[2] & meta$Latitude >= ylim[1] & meta$Latitude <= ylim[2], ]
  ggplot() + geom_raster(data = df, aes(x, y, fill = value)) +
    scale_fill_viridis_c(name = fill_name, na.value = "grey80", limits = limits) +
    geom_path(data = world, aes(long, lat, group = group), colour = "grey25", linewidth = 0.25) +
    geom_point(data = sm, aes(Longitude, Latitude), shape = 21, fill = "white", size = pt_size, stroke = 0.3) +
    coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = "Longitude", y = "Latitude", title = title) +
    theme_bw(base_size = 10) + theme(panel.background = element_rect(fill = "#DCEBF5"), panel.grid = element_blank())
}
cap <- "Grey = land whose current climate lies outside the sampled climate (ExDet), masked; blue = sea; points = accessions"
for (sid in names(results)) {
  R <- results[[sid]]; lab <- paste0(pretty_ssp(R$ssp), ", ", R$period, " (mean of ", R$n_gcm, " GCMs)")
  m <- offset_map(c(R$r[["offset_mean"]], R$r[["extrapolation_masked"]]), EV[1:2], EV[3:4],
                  paste("Genomic offset:", lab), "Genomic offset\n(low \u2192 high)", lims) + labs(caption = cap)
  asp <- diff(EV[3:4]) / diff(EV[1:2])
  save_fig(paste0("offset_map_", sid), m, width = 9, height = max(4.5, 9 * asp * 0.95))
  s <- offset_map(c(R$r[["offset_sd_among_GCMs"]], R$r[["extrapolation_masked"]]), EV[1:2], EV[3:4],
                  paste("Uncertainty (SD among GCMs):", lab), "SD among\nGCMs") + labs(caption = cap)
  save_fig(paste0("offset_sd_map_", sid), s, width = 9, height = max(4.5, 9 * asp * 0.95))
  # regional panels cut from the same map (same colour scale)
  if (length(CFG$regions)) {
    rp <- list()
    for (rn in names(CFG$regions)) {
      b <- CFG$regions[[rn]]
      xl <- c(max(b[1], EV[1]), min(b[2], EV[2])); yl <- c(max(b[3], EV[3]), min(b[4], EV[4]))
      if (diff(xl) <= 0 || diff(yl) <= 0) next
      if (!any(meta$Longitude >= xl[1] & meta$Longitude <= xl[2] & meta$Latitude >= yl[1] & meta$Latitude <= yl[2])) next
      rp[[rn]] <- offset_map(c(R$r[["offset_mean"]], R$r[["extrapolation_masked"]]), xl, yl, rn,
                             "Genomic offset\n(low \u2192 high)", if (is.null(lims)) range(values(R$r[["offset_mean"]]), na.rm = TRUE) else lims, 1.8)
    }
    if (length(rp)) {
      pr <- patchwork::wrap_plots(rp, nrow = 1, guides = "collect") +
        patchwork::plot_annotation(title = paste("Genomic offset by region:", lab), caption = cap)
      save_fig(paste0("offset_map_regions_", sid), pr, width = 6 * length(rp), height = 6.5)
    }
  }
}

# ==============================================================================
# 6  PER-ACCESSION OFFSETS
# ==============================================================================
if (!all(is.na(acc$Population))) {
  gs <- acc %>% group_by(scenario, Population) %>%
    summarise(n = n(), median_offset = median(offset_mean, na.rm = TRUE), mean_offset = mean(offset_mean, na.rm = TRUE),
              max_offset = max(offset_mean, na.rm = TRUE), .groups = "drop")
  write.csv(gs, out("offset_by_population.csv"), row.names = FALSE); print(gs)
  pb <- ggplot(acc, aes(reorder(Population, offset_mean, median), offset_mean)) +
    geom_boxplot(fill = "grey90", outlier.shape = NA) + geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
    facet_wrap(~ scenario, scales = "free_y") + labs(x = NULL, y = "Genomic offset at the accession's site",
                                                    title = "Genomic offset by population") + theme_bw()
  save_fig("offset_by_population", pb, width = 4 + 2.5 * n_distinct(acc$scenario), height = 5)
}
for (sid in unique(acc$scenario)) {                 # direction check, as in the yam script
  d <- acc[acc$scenario == sid, ]
  msg(sid, ": offset ~ latitude slope", signif(coef(lm(offset_mean ~ Latitude, d))[2], 3),
      "(p =", signif(summary(lm(offset_mean ~ Latitude, d))$coefficients[2, 4], 2), ")")
}
writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
msg("Done ->", normalizePath(CFG$out_dir))
