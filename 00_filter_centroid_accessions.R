################################################################################
# 00_filter_centroid_accessions.R
#
# Many accessions share coordinates that are almost certainly geocoded country names
# (country centroids), not collection sites. Environmental values extracted at those
# points do not describe where the accessions were collected. This script keeps ONE
# accession per putative centroid and writes filtered copies of every input file to
# data_filtered/, which scripts 01-09 read.
#
# Putative centroids
#   - CFG$centroids: coordinates treated as country centroids (always collapsed).
#     Defaults: Nepal (57 accessions), India (52) and Kenya (48) in the full panel.
#   - Any other coordinate shared by >= CFG$flag_shared accessions is REPORTED for review
#     (data_filtered/shared_coordinates_report.csv); set CFG$collapse_flagged = TRUE to
#     collapse those as well.
# Which accession is kept: the one with the fewest missing genotype calls in the VCF
#   (ties broken alphabetically), so the most informative genotype is retained.
#
# Output (data_filtered/): filtered copies of the VCF and every CSV/TXT listed in
#   CFG$tables (rows for removed accessions dropped, row order kept), plus
#   kept_accessions.csv, removed_accessions.csv, centroid_summary.csv and
#   shared_coordinates_report.csv.
################################################################################

CFG <- list(
  vcf          = "finger_millet_404.vcf",
  passport_csv = "millet_passport_test.csv",        # used to find the coordinates
  # every per-accession table used by scripts 01-09 (missing files are skipped with a message)
  tables       = c("millet_passport_test.csv",           # passport: ID, Longitude, Latitude, Population ...
                   "sample_coordinates.csv",             # coordinates for genomic prediction (07)
                   "millet_extracted_bioclim_values.csv",# environment at collection sites (04, 05, 07)
                   "test_for_numeric.csv",               # numeric genotypes for TESS3 (03)
                   "metadata_for_PCA.csv",               # PCA groups (02)
                   "millet_core_n100.csv"),              # optional supplied core (07); 01 recreates it
  centroids    = data.frame(label = c("Nepal centroid", "India centroid", "Kenya centroid"),
                            Latitude = c(28.108, 22.351, 1.442), Longitude = c(84.092, 78.668, 38.431)),
  coord_digits = 3,          # coordinates are matched after rounding to this many decimals
  flag_shared  = 10,         # report other coordinates shared by at least this many accessions
  collapse_flagged = FALSE,  # TRUE = also keep only one accession at those coordinates
  out_dir      = "data_filtered",
  filtered_vcf_name = "finger_millet_filtered.vcf"
)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

# ---- VCF: sample names and missingness per sample ---------------------------------------
con <- if (grepl("\\.gz$", CFG$vcf)) gzfile(CFG$vcf, "r") else file(CFG$vcf, "r")
hdr <- character(0)
repeat { l <- readLines(con, n = 1); if (!length(l)) stop("No #CHROM line in ", CFG$vcf)
         hdr <- c(hdr, l); if (startsWith(l, "#CHROM")) break }
body <- readLines(con); close(con)
cols <- strsplit(hdr[length(hdr)], "\t")[[1]]; samples <- trimws(cols[-(1:9)])
gt <- do.call(rbind, strsplit(body, "\t", fixed = TRUE))
miss <- colMeans(substr(gt[, -(1:9), drop = FALSE], 1, 1) == ".")
names(miss) <- samples
msg("VCF:", length(samples), "samples x", length(body), "variants")

# ---- passport coordinates ---------------------------------------------------------------
id_col_for <- function(d) {                      # column whose values best match the VCF sample IDs
  hit <- vapply(d, function(x) sum(trimws(as.character(x)) %in% samples), 0)
  if (max(hit) == 0) stop("No column matches the VCF sample IDs")
  names(d)[which.max(hit)]
}
pp <- read.csv(CFG$passport_csv, check.names = FALSE, stringsAsFactors = FALSE)
idc <- id_col_for(pp); nm <- tolower(names(pp))
latc <- names(pp)[grep("^lat", nm)[1]]; lonc <- names(pp)[grep("^lon", nm)[1]]
pp$.id <- trimws(pp[[idc]]); pp <- pp[pp$.id %in% samples, ]
pp$.coord <- ifelse(is.na(pp[[latc]]) | is.na(pp[[lonc]]), NA,
                    paste(round(pp[[latc]], CFG$coord_digits), round(pp[[lonc]], CFG$coord_digits), sep = ","))

cent_keys <- paste(round(CFG$centroids$Latitude, CFG$coord_digits), round(CFG$centroids$Longitude, CFG$coord_digits), sep = ",")
shared <- sort(table(pp$.coord[!is.na(pp$.coord)]), decreasing = TRUE)
flagged <- setdiff(names(shared)[shared >= CFG$flag_shared], cent_keys)
report <- data.frame(coordinate = names(shared), n_accessions = as.integer(shared))
report$status <- ifelse(report$coordinate %in% cent_keys, "putative country centroid (collapsed)",
                 ifelse(report$coordinate %in% flagged,
                        if (CFG$collapse_flagged) "shared by many (collapsed)" else "shared by many - review", "ok"))
write.csv(report[report$n_accessions > 1, ], file.path(CFG$out_dir, "shared_coordinates_report.csv"), row.names = FALSE)
if (length(flagged)) msg("Other coordinates shared by >=", CFG$flag_shared, "accessions (",
                         if (CFG$collapse_flagged) "collapsed" else "reported only", "):", paste(flagged, collapse = " | "))

collapse <- c(cent_keys, if (CFG$collapse_flagged) flagged)
keep <- rep(TRUE, nrow(pp)); summ <- list()
for (k in collapse) {
  i <- which(pp$.coord == k)
  if (!length(i)) { msg("Centroid", k, "not found in the passport data"); next }
  ids <- pp$.id[i]; o <- order(miss[ids], ids); kept <- ids[o[1]]
  keep[i] <- pp$.id[i] == kept
  lab <- if (k %in% cent_keys) CFG$centroids$label[match(k, cent_keys)] else "shared coordinate"
  summ[[k]] <- data.frame(coordinate = k, label = lab, n_accessions = length(i), kept_accession = kept,
                          kept_missing_rate = round(miss[kept], 4), n_removed = length(i) - 1)
}
summ <- do.call(rbind, summ); print(summ, row.names = FALSE)
write.csv(summ, file.path(CFG$out_dir, "centroid_summary.csv"), row.names = FALSE)

keep_ids <- pp$.id[keep]
# accessions without coordinates in the passport file are kept (they are handled by each script)
no_pass <- setdiff(samples, pp$.id); keep_ids <- union(keep_ids, no_pass)
removed <- pp[!keep, c(idc, latc, lonc)]; removed$reason <- "duplicate at putative centroid"
write.csv(data.frame(SampleID = samples[samples %in% keep_ids]), file.path(CFG$out_dir, "kept_accessions.csv"), row.names = FALSE)
write.csv(removed, file.path(CFG$out_dir, "removed_accessions.csv"), row.names = FALSE)
msg("Accessions kept:", sum(samples %in% keep_ids), "| removed:", nrow(removed))

# ---- filtered VCF -----------------------------------------------------------------------------
kc <- c(1:9, 9 + which(samples %in% keep_ids))
out_vcf <- file.path(CFG$out_dir, CFG$filtered_vcf_name)
writeLines(c(hdr[-length(hdr)], paste(cols[kc], collapse = "\t"),
             apply(gt[, kc, drop = FALSE], 1, paste, collapse = "\t")), out_vcf)
msg("Filtered VCF ->", out_vcf)

# ---- filtered tables (same file names; rows of removed accessions dropped) ----------------------
for (f in CFG$tables) {
  if (!file.exists(f)) { msg("Not found, skipped:", f); next }
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  ic <- tryCatch(id_col_for(d), error = function(e) NA)
  if (is.na(ic)) { msg("No sample-ID column found, copied unchanged:", f); file.copy(f, file.path(CFG$out_dir, basename(f)), overwrite = TRUE); next }
  d2 <- d[!trimws(d[[ic]]) %in% removed[[idc]], , drop = FALSE]
  write.csv(d2, file.path(CFG$out_dir, basename(f)), row.names = FALSE)
  msg(basename(f), ":", nrow(d), "->", nrow(d2), "rows")
}
msg("Done. Run scripts 01-09 next; they read from", CFG$out_dir)
