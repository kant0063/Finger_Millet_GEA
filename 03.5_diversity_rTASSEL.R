################################################################################
# 03.5_diversity_rTASSEL.R  -  CENTROID-FILTERED PIPELINE (run after 00)
#
# Sliding-window sequence diversity with rTASSEL (TASSEL SequenceDiversityPlugin) on the
# filtered VCF: nucleotide diversity (pi per bp), Watterson's theta per bp and Tajima's D,
# computed separately for each chromosome so that no window spans two chromosomes.
# Windows are defined in SNP sites (TASSEL's definition), not base pairs.
#
# Outputs (results/03.5_diversity/):
#   finger_millet_diversity.txt              all windows (TASSEL Diversity report columns)
#   diversity_summary_by_chromosome.csv      mean / median pi, theta and Tajima's D per chromosome
#   diversity_summary_by_subgenome.csv       A vs B sub-genome (+ Wilcoxon test)
#   figures/nucleotide_diversity_by_chromosome.(pdf|png)
#   figures/tajimas_D_by_chromosome.(pdf|png)
#
# Note: DArTseq data contain only variable sites, so pi per bp is a relative measure for
# comparing windows, chromosomes and sub-genomes, not an absolute nucleotide diversity.
################################################################################

CFG <- list(
  vcf          = "data_filtered/finger_millet_filtered.vcf",
  window_sites = 20,           # SNPs per window
  step_sites   = 10,           # SNPs between window starts
  min_maf      = 0,            # optional TASSEL site filter before the analysis (0 = keep all sites)
  java_memory  = "-Xmx4g",
  chr_map      = setNames(paste0(rep(1:9, each = 2), c("A", "B")), 1:18),   # TASSEL 1-18 -> 1A..9B; 19 = scaffolds ("Un")
  include_unplaced = FALSE,    # TRUE = also compute windows on the concatenated scaffolds (chromosome 19)
  out_dir      = "results/03.5_diversity",
  copy_to      = NA,           # optional extra copy of the window table (path); NA = none
  dpi          = 600
)

options(java.parameters = CFG$java_memory)               # must be set before rJava starts
need <- function(p) { m <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
if (length(m)) install.packages(m, repos = "https://cloud.r-project.org") }
need(c("ggplot2", "dplyr"))
if (!requireNamespace("rTASSEL", quietly = TRUE)) {
  need("remotes"); remotes::install_github("maize-genetics/rTASSEL", build_vignettes = FALSE)
}
if (!requireNamespace("rTASSEL", quietly = TRUE))
  stop("rTASSEL is not installed or cannot be loaded (it needs Java and rJava). Restart R and run the script again.")
suppressPackageStartupMessages({ library(rTASSEL); library(ggplot2); library(dplyr) })
# rTASSEL functions are called as rTASSEL::fn, and checked by name, because function names
# differ between rTASSEL versions (e.g. readGenotypeTableFromPath() is being replaced)
rt_has <- function(fn) fn %in% getNamespaceExports("rTASSEL")
msg0 <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
msg0("rTASSEL version", as.character(utils::packageVersion("rTASSEL")))
for (fn in c("seqDiversity")) if (!rt_has(fn))
  stop("This rTASSEL version has no ", fn, "(); install the GitHub version: remotes::install_github(\"maize-genetics/rTASSEL\")")
dir.create(file.path(CFG$out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
save_gg <- function(name, p, width, height) {
  ggsave(file.path(CFG$out_dir, "figures", paste0(name, ".pdf")), p, width = width, height = height)
  ggsave(file.path(CFG$out_dir, "figures", paste0(name, ".png")), p, width = width, height = height, dpi = CFG$dpi, bg = "white")
  print(p)
}

# ---- read the filtered VCF with rTASSEL ------------------------------------------------------
# optional TASSEL log file; the arguments differ between rTASSEL versions, so try both forms
if (rt_has("startLogger")) {
  lg <- names(formals(rTASSEL::startLogger))
  try(if ("filePath" %in% lg) rTASSEL::startLogger(filePath = file.path(tempdir(), "rtassel_log.txt"), verbose = FALSE)
      else rTASSEL::startLogger(fullPath = tempdir(), fileName = "rtassel_log.txt"), silent = TRUE)
}
read_geno <- function() {
  if (rt_has("readGenotypeTableFromPath"))
    return(suppressMessages(rTASSEL::readGenotypeTableFromPath(path = CFG$vcf, keepDepth = FALSE)))
  if (rt_has("readGenotype")) return(suppressMessages(rTASSEL::readGenotype(CFG$vcf)))
  stop("No genotype reader found in this rTASSEL version")
}
# TASSEL's Java libraries (JARs) must be on the Java class path. After a first-time download with
# setupTASSEL(), Java is already running without them: add the cached JARs to the running session
# and retry; if that is not enough, ask for an R restart (the JARs are then loaded automatically).
add_cached_jars <- function() {
  jars <- list.files(tools::R_user_dir("rTASSEL", "cache"), pattern = "\\.jar$", recursive = TRUE, full.names = TRUE)
  if (length(jars)) rJava::.jaddClassPath(jars)
  length(jars)
}
restart_msg <- paste("TASSEL is installed but this R session started Java without it.",
                     "Restart R (RStudio: Session > Restart R) and run 03.5_diversity_rTASSEL.R again.")
geno <- tryCatch(read_geno(), error = function(e) {
  if (add_cached_jars() > 0) {                       # JARs already downloaded earlier: add them and retry
    g <- tryCatch(read_geno(), error = function(e2) NULL)
    if (!is.null(g)) return(g)
  }
  if (!rt_has("setupTASSEL")) stop("Could not read the VCF with rTASSEL: ", conditionMessage(e))
  msg("TASSEL libraries not found; running setupTASSEL() (one-time download)")
  rTASSEL::setupTASSEL()
  add_cached_jars()
  tryCatch(read_geno(), error = function(e3) stop(restart_msg, call. = FALSE))
})
if (!methods::.hasSlot(geno, "jGenotypeTable"))
  stop("rTASSEL returned an object of class ", class(geno)[1], " without a TASSEL genotype table; ",
       "install the GitHub version of rTASSEL, which provides readGenotypeTableFromPath() and seqDiversity()")
if (CFG$min_maf > 0) {
  geno <- rTASSEL::filterGenotypeTableSites(geno, siteMinCount = 0, siteMinAlleleFreq = CFG$min_maf,
                                            siteMaxAlleleFreq = 1, siteRangeFilterType = "none")
}
jg <- geno@jGenotypeTable
msg("Genotype table:", jg$numberOfTaxa(), "accessions x", jg$numberOfSites(), "sites")

# ---- sliding-window diversity, one chromosome at a time -----------------------------------------
as_vec <- function(x) if (methods::is(x, "jarrayRef")) rJava::.jevalArray(x) else x
chroms <- as_vec(jg$chromosomes())
win <- list()
for (ch in chroms) {
  nm <- ch$getName()
  if (!CFG$include_unplaced && !nm %in% names(CFG$chr_map)) { msg("  chromosome", nm, "(unplaced scaffolds) skipped"); next }
  fl <- as.integer(as_vec(jg$firstLastSiteOfChromosome(ch)))      # first and last site index
  n_sites <- fl[2] - fl[1] + 1
  if (n_sites < CFG$window_sites) { msg("  chromosome", nm, ": only", n_sites, "sites, skipped"); next }
  d <- tryCatch(rTASSEL::seqDiversity(geno, startSite = fl[1], endSite = fl[2], slidingWindowAnalysis = TRUE,
                                      stepSize = CFG$step_sites, windowSize = CFG$window_sites)$Diversity,
                error = function(e) { msg("  chromosome", nm, "failed:", conditionMessage(e)); NULL })
  if (is.null(d)) next
  d <- as.data.frame(d)
  if ("Site_Type" %in% names(d) && length(unique(d$Site_Type)) > 1) d <- d[d$Site_Type == d$Site_Type[1], ]
  d$Chromosome <- nm                                      # TASSEL chromosome name (1-19)
  win[[nm]] <- d
  msg("  chromosome", nm, ":", n_sites, "sites,", nrow(d), "windows")
}
div <- bind_rows(win)
if (!nrow(div)) stop("No diversity windows computed")
num <- c("StartChrPosition", "EndChrPosition", "StartSite", "EndSite", "MidSite", "SiteCount", "AvgSiteCount",
         "SegSites", "PiPerBP", "ThetaPerBP", "TajimaD")
for (cc in intersect(num, names(div))) div[[cc]] <- suppressWarnings(as.numeric(gsub(",", "", div[[cc]])))
write.table(div, file.path(CFG$out_dir, "finger_millet_diversity.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.na(CFG$copy_to)) write.table(div, CFG$copy_to, sep = "\t", quote = FALSE, row.names = FALSE)
msg("Windows:", nrow(div), "->", file.path(CFG$out_dir, "finger_millet_diversity.txt"))

# ---- summaries: chromosome and sub-genome ----------------------------------------------------------
div$Chromosome <- as.character(div$Chromosome)
div$Chr_name <- ifelse(div$Chromosome %in% names(CFG$chr_map), CFG$chr_map[div$Chromosome], "Un")
div$Subgenome <- ifelse(div$Chr_name == "Un", "Un", substring(div$Chr_name, nchar(div$Chr_name)))
chr_order <- unique(c(unname(CFG$chr_map), "Un"))
by_chr <- div %>% group_by(Chr_name, Subgenome) %>%
  summarise(n_windows = n(), mean_pi_per_bp = mean(PiPerBP, na.rm = TRUE), median_pi_per_bp = median(PiPerBP, na.rm = TRUE),
            mean_theta_per_bp = mean(ThetaPerBP, na.rm = TRUE), mean_tajimas_D = mean(TajimaD, na.rm = TRUE),
            span_Mb = (max(EndChrPosition) - min(StartChrPosition)) / 1e6, .groups = "drop") %>%
  mutate(Chr_name = factor(Chr_name, levels = chr_order)) %>% arrange(Chr_name)
write.csv(by_chr, file.path(CFG$out_dir, "diversity_summary_by_chromosome.csv"), row.names = FALSE)
by_sg <- div %>% filter(Subgenome %in% c("A", "B")) %>% group_by(Subgenome) %>%
  summarise(n_windows = n(), mean_pi_per_bp = mean(PiPerBP, na.rm = TRUE), median_pi_per_bp = median(PiPerBP, na.rm = TRUE),
            mean_tajimas_D = mean(TajimaD, na.rm = TRUE), .groups = "drop")
ab <- div %>% filter(Subgenome %in% c("A", "B"))
if (length(unique(ab$Subgenome)) == 2) {
  by_sg$wilcoxon_p_pi_A_vs_B <- wilcox.test(PiPerBP ~ Subgenome, data = ab)$p.value
  by_sg$wilcoxon_p_D_A_vs_B  <- wilcox.test(TajimaD ~ Subgenome, data = ab)$p.value
}
write.csv(by_sg, file.path(CFG$out_dir, "diversity_summary_by_subgenome.csv"), row.names = FALSE)
print(by_sg)

# ---- figures (same layout as script 03) ---------------------------------------------------------------
div <- div %>%
  mutate(Chr_name = factor(Chr_name, levels = chr_order[chr_order %in% Chr_name])) %>%
  arrange(Chr_name, StartSite) %>%
  mutate(window_index = row_number())                      # windows in genome order
chrom_bands <- div %>% group_by(Chr_name) %>%
  summarise(xmin = min(window_index) - 0.5, xmax = max(window_index) + 0.5, mid = mean(range(window_index)), .groups = "drop") %>%
  mutate(ymin = -Inf, ymax = Inf)
by_chrom_plot <- function(yvar, ylab, hline = NULL) {
  p <- ggplot(div, aes(x = window_index, y = .data[[yvar]], group = Chr_name, color = Chr_name)) +
    geom_rect(data = chrom_bands, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Chr_name),
              inherit.aes = FALSE, alpha = 0.2) +
    geom_line(linewidth = 0.4) +
    scale_x_continuous(breaks = chrom_bands$mid, labels = chrom_bands$Chr_name, expand = c(0.005, 0)) +
    labs(x = "Chromosome", y = ylab, color = "Chromosome", fill = "Chromosome",
         caption = paste0(nrow(div), " windows of ", CFG$window_sites, " SNPs (step ", CFG$step_sites, "); ",
                          jg$numberOfTaxa(), " accessions")) +
    theme_bw() +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(size = 8))
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = 2, colour = "grey40")
  p
}
save_gg("nucleotide_diversity_by_chromosome", by_chrom_plot("PiPerBP", expression("Nucleotide diversity (" * pi * " per bp)")), 12, 4.5)
save_gg("tajimas_D_by_chromosome", by_chrom_plot("TajimaD", "Tajima's D", hline = 0), 12, 4.5)
msg("Done ->", normalizePath(CFG$out_dir))