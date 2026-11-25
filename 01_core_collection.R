# ==== 01_core_collection.R - CENTROID-FILTERED PIPELINE (one accession per putative country centroid) ====
# Core Hunter core of 100 accessions. Output: results/01_core_collection/millet_core_n100.csv (used by 07).
# Inputs are read from data_filtered/ (created by 00_filter_centroid_accessions.R).

########################################
# Step 1: Read VCF and build genotype matrix
########################################
library(vcfR)
library(corehunter)

# Step 1: Read VCF file
vcf_file <- "data_filtered/finger_millet_filtered.vcf"
out_dir  <- "results/01_core_collection"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
vcf <- read.vcfR(vcf_file)

# Extract genotype calls (e.g. "0/0", "0/1", "1/1") and convert to 0/1/2 alt-allele
# dosage, equivalent to allel.GenotypeArray(...).to_n_alt() in the Python version
gt_char <- extract.gt(vcf, element = "GT")

gt_to_dosage <- function(x) {
  alleles <- strsplit(x, "[/|]")
  vapply(alleles, function(a) {
    if (length(a) < 2 || any(a == ".")) return(NA_real_)
    sum(as.numeric(a))
  }, numeric(1))
}

# genotypes: rows = SNPs, columns = individuals (dosage-coded)
genotypes <- apply(gt_char, 2, gt_to_dosage)
rownames(genotypes) <- rownames(gt_char)

########################################
# Step 2: Core collection selection with corehunter
########################################
# corehunter expects individuals as rows and markers as columns, so transpose
geno_matrix <- t(genotypes)

# Wrap as Core Hunter genotype data. Values are 0/1/2 dosage counts, i.e. the
# "biparental" format (equivalent to allele-dosage / SNP data for a biparental
# or diploid population).
core_data <- genotypes(data = geno_matrix, format = "biparental")

# Set the number of individuals to select (e.g., 100)
num_core_samples <- 100

# Use the entry-to-nearest-entry (EN) objective with the Modified Rogers (MR)
# genetic distance -- this maximizes the minimum genetic distance between each
# selected accession and its nearest neighbor in the core, which is the same
# "maximum diversity" criterion implemented by the greedy farthest-point
# algorithm in the original Python script.
set.seed(42)  # for reproducibility, since Core Hunter's search is stochastic
obj <- objective(type = "EN", measure = "MR")
core <- sampleCore(core_data, obj = obj, size = num_core_samples)

# Map the result back to sample names (Core Hunter returns ids by default)
core_sample <- core$sel

# Print the core sample names
cat("\nCore Sample IDs:\n")
cat(paste(core_sample, collapse = "\n"), "\n")

# Save the core sample (SRR names) to CSV
core_sample_data <- data.frame(Sample_ID = core_sample)
write.csv(core_sample_data,
          file.path(out_dir, "millet_core_n100.csv"),
          row.names = FALSE)
