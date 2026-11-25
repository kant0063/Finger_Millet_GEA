# ==== 03_TESS3_ancestry.R - CENTROID-FILTERED PIPELINE (one accession per putative country centroid) ====
# TESS3 K = 1-8 with cross-validation, Q matrix K = 4, kriged maps, outliers.
# (Nucleotide diversity is now computed and plotted by 03.5_diversity_rTASSEL.R.)
# Inputs are read from data_filtered/ (created by 00_filter_centroid_accessions.R).

################################################################################
# finger_millet_interpolation.R
#
# TESS3 spatial ancestry for finger millet: K = 1-8, cross-validation, ancestry (Q)
# matrix for K = 4, kriged ancestry maps, TESS3 outlier loci, and extraction of
# WorldClim bioclim values at the collection sites.
#
# Changes from the previous version
#  - Cross-validation is now actually computed: tess3() is run with mask = 0.05 (5% of
#    genotypes held out) and CFG$rep repetitions. With the default mask = 0 TESS3 computes
#    no cross-validation, and plot(tess3.obj) shows the training error (RMSE), not a
#    cross-validation score.
#  - New tables: TESS3_cross_validation_K1-8.csv and TESS3_Qmatrix_K4.csv (+ p-values).
#  - Accessions without coordinates are found automatically (no hard-coded row numbers).
#  - The Manhattan plot now plots -log10(p) (previously raw p under a -log10 label).
#  - All figures are saved as PDF (vector) and 600-dpi PNG.
#
# Inputs : genotype CSV (first column = accession name, then SNPs coded 0/1/2, NA = missing)
#          passport CSV (first column = accession name, then Longitude, Latitude)
# Outputs: CFG$out_dir/tables and CFG$out_dir/figures
################################################################################

# ---- settings ------------------------------------------------------------------
CFG <- list(
  genotype_csv = "data_filtered/test_for_numeric.csv",
  passport_csv = "data_filtered/millet_passport_test.csv",
  K            = 1:8,
  K_main       = 4,          # K used for the Q matrix, maps and outlier test
  rep          = 5,          # TESS3 repetitions per K (best run kept for Q; all used for CV)
  mask         = 0.05,       # proportion of genotypes held out for cross-validation
  cores        = 4,
  fdr          = 0.2,        # BH threshold for TESS3 outlier loci
  colors       = c("tomato", "orange", "lightblue", "olivedrab"),
  map_xlim     = c(-20, 100),
  map_ylim     = c(-25, 40),
  bioclim_dir  = "/Users/mikey/Downloads/wc2.1_30s_bio",   # WorldClim 2.1 30 s GeoTIFFs; NA = skip
  out_dir      = "results/03_TESS3",
  dpi          = 600,
  seed         = 2026
)

# install.packages("devtools"); devtools::install_github("bcm-uga/TESS3_encho_sen")
suppressPackageStartupMessages({ library(tess3r); library(maps); library(ggplot2) })
set.seed(CFG$seed)
for (d in c("tables", "figures")) dir.create(file.path(CFG$out_dir, d), recursive = TRUE, showWarnings = FALSE)
tab <- function(f) file.path(CFG$out_dir, "tables", f)
fig <- function(f) file.path(CFG$out_dir, "figures", f)

# base-graphics figure -> PDF + high-resolution PNG
save_base <- function(name, draw, width = 7, height = 5) {
  pdf(fig(paste0(name, ".pdf")), width = width, height = height); draw(); invisible(dev.off())
  png(fig(paste0(name, ".png")), width = width, height = height, units = "in", res = CFG$dpi); draw(); invisible(dev.off())
  draw()                                             # also show on screen
}
save_gg <- function(name, p, width = 9, height = 6) {
  ggsave(fig(paste0(name, ".pdf")), p, width = width, height = height)
  ggsave(fig(paste0(name, ".png")), p, width = width, height = height, dpi = CFG$dpi, bg = "white")
  print(p)
}

# ---- data ----------------------------------------------------------------------
genotype    <- read.csv(CFG$genotype_csv, check.names = FALSE)
coordinates <- read.csv(CFG$passport_csv, check.names = FALSE)
stopifnot("Accession order differs between genotype and passport files" =
            all(genotype[[1]] == coordinates[[1]]))
accession_names <- coordinates[[1]]
genotype    <- as.matrix(genotype[, -1])
coordinates <- as.matrix(coordinates[, c("Longitude", "Latitude")])

# accessions without a georeference are removed (found from the data, not hard-coded)
no_geo <- !complete.cases(coordinates)
cat("Accessions without coordinates (removed):", sum(no_geo), "\n"); print(accession_names[no_geo])
genotype <- genotype[!no_geo, , drop = FALSE]; coordinates <- coordinates[!no_geo, , drop = FALSE]
accession_names <- accession_names[!no_geo]
cat("Accessions analysed:", nrow(genotype), "| SNPs:", ncol(genotype), "\n")

save_base("sampling_locations", function() {
  plot(coordinates, pch = 19, cex = .5, xlab = "Longitude (°E)", ylab = "Latitude (°N)")
  map(add = TRUE, interior = FALSE)
})

# ---- TESS3 ---------------------------------------------------------------------
tess3.obj <- tess3(X = genotype, coord = coordinates, K = CFG$K, method = "projected.ls", ploidy = 2,
                   rep = CFG$rep, mask = CFG$mask, keep = "best", openMP.core.num = CFG$cores)

# cross-validation scores for every K (summarised over repetitions)
cv <- do.call(rbind, lapply(seq_along(tess3.obj), function(i) {
  r <- tess3.obj[[i]]
  data.frame(K = r$K,
             crossvalid_rmse_median = median(r$crossvalid.rmse), crossvalid_rmse_min = min(r$crossvalid.rmse),
             crossvalid_rmse_max = max(r$crossvalid.rmse),
             crossvalid_crossentropy_median = median(r$crossvalid.crossentropy),
             crossvalid_crossentropy_min = min(r$crossvalid.crossentropy),
             training_rmse_median = median(r$rmse), training_crossentropy_median = median(r$crossentropy),
             n_rep = length(r$rmse))
}))
write.csv(cv, tab(paste0("TESS3_cross_validation_K", min(CFG$K), "-", max(CFG$K), ".csv")), row.names = FALSE)
print(cv)

save_base("TESS3_cross_validation", function() {
  plot(tess3.obj, crossvalid = TRUE, crossentropy = FALSE, pch = 19, col = "blue",
       xlab = "Number of ancestral populations (K)", ylab = "Cross-validation score (RMSE)")
})

# ---- Q matrix for K_main -----------------------------------------------------------
q.matrix <- qmatrix(tess3.obj, K = CFG$K_main)
Q <- as.data.frame(unclass(q.matrix)); colnames(Q) <- paste0("Q", seq_len(ncol(Q)))
q_out <- data.frame(Accession.name = accession_names, Longitude = coordinates[, "Longitude"],
                    Latitude = coordinates[, "Latitude"], Q,
                    main_cluster = paste0("Q", max.col(as.matrix(Q))), max_ancestry = apply(Q, 1, max))
write.csv(q_out, tab(paste0("TESS3_Qmatrix_K", CFG$K_main, ".csv")), row.names = FALSE)
cat("Main cluster sizes:\n"); print(table(q_out$main_cluster))

my.palette <- CreatePalette(CFG$colors[seq_len(CFG$K_main)], 9)
save_base(paste0("TESS3_ancestry_barplot_K", CFG$K_main), function() {
  bp <- barplot(q.matrix, border = NA, space = 0, main = paste0("Ancestry matrix (K = ", CFG$K_main, ")"),
                xlab = "Individuals", ylab = "Ancestry proportions", col.palette = my.palette)
  axis(1, at = 1:nrow(q.matrix), labels = accession_names[bp$order], las = 3, cex.axis = .25)
}, width = 12, height = 4.5)

save_base(paste0("TESS3_ancestry_map_kriged_K", CFG$K_main), function() {
  plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10),
       main = paste0("Ancestry coefficients (K = ", CFG$K_main, ")"), xlab = "Longitude", ylab = "Latitude",
       resolution = c(300, 300), cex = .4, col.palette = my.palette)
}, width = 9, height = 6)

if (requireNamespace("rworldmap", quietly = TRUE)) {
  map.polygon <- rworldmap::getMap(resolution = "low")
  pl <- ggtess3Q(q.matrix, coordinates, map.polygon = map.polygon, col.palette = my.palette) +
    geom_path(data = map.polygon, aes(x = long, y = lat, group = group), linewidth = 0.2) +
    coord_equal(xlim = CFG$map_xlim, ylim = CFG$map_ylim) +   # zoom without cutting country outlines
    geom_point(data = as.data.frame(coordinates), aes(x = Longitude, y = Latitude), size = 0.3) +
    labs(x = "Longitude", y = "Latitude", title = paste0("Interpolated ancestry (K = ", CFG$K_main, ")")) + theme_bw()
  save_gg(paste0("TESS3_ancestry_map_K", CFG$K_main), pl, width = 10, height = 6)
}

# ---- TESS3 outlier loci for K_main --------------------------------------------------
p.values <- pvalue(tess3.obj, K = CFG$K_main)
adj.p <- p.adjust(p.values, method = "BH")
candidates <- which(adj.p <= CFG$fdr)
cat("Outlier loci at BH-adjusted p <=", CFG$fdr, ":", length(candidates), "\n")
write.csv(data.frame(locus = colnames(genotype), p_value = p.values, p_BH = adj.p, outlier = adj.p <= CFG$fdr),
          tab(paste0("TESS3_pvalues_K", CFG$K_main, ".csv")), row.names = FALSE)
save_base(paste0("TESS3_pvalue_histogram_K", CFG$K_main), function() hist(p.values, col = "lightblue", main = "TESS3 p-values", xlab = "p-value"))
save_base(paste0("TESS3_manhattan_K", CFG$K_main), function() {
  plot(-log10(p.values), main = paste0("TESS3 outlier test (K = ", CFG$K_main, ")"), xlab = "Locus index",
       ylab = expression(-log[10](italic(p))), cex = .3, col = "grey", pch = 19)
  points(candidates, -log10(p.values)[candidates], pch = 19, cex = .3, col = "blue")
}, width = 10, height = 4.5)

# ---- WorldClim bioclim values at the collection sites -------------------------------
if (!is.na(CFG$bioclim_dir) && dir.exists(CFG$bioclim_dir)) {
  library(terra)
  pp <- read.csv(CFG$passport_csv, check.names = FALSE)
  r <- rast(list.files(CFG$bioclim_dir, pattern = "\\.tif$", full.names = TRUE))
  pts <- vect(pp[, c("Longitude", "Latitude")][complete.cases(pp[, c("Longitude", "Latitude")]), ],
              geom = c("Longitude", "Latitude"), crs = crs(r))
  vals <- terra::extract(r, pts)
  result <- cbind(pp[complete.cases(pp[, c("Longitude", "Latitude")]), ], vals[, -1])
  write.csv(result, tab("finger_millet_bioclim.csv"), row.names = FALSE)
}
cat("Done ->", normalizePath(CFG$out_dir), "\n")