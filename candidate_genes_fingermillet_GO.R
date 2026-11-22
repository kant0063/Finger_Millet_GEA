################################################################################
# candidate_genes_fingermillet_GO.R                        Discussion (candidate genes)
# Finger millet (Eleusine coracana) candidate genes and GO enrichment from the eGWAS output.
#
# Reads significant_markers.csv written by environmental_gwas_complete.R directly
# (columns marker, chrom, pos, variable, model, q_value, PVE_percent, beta,
# effect_direction ...). The older format (SNP, CHROM, POS, trait, set, PVE, q) still works.
#
# Gene Ontology enrichment (all in R, no Bioconductor needed):
#   - GO terms per gene from the Phytozome annotation_info file (GO column).
#   - Ontology structure from go-basic.obo (found in input_dir, or downloaded from
#     http://purl.obolibrary.org/obo/go/go-basic.obo on first run). Annotations are
#     propagated to all ancestor terms (is_a and part_of), roots excluded.
#   - Background ("universe"): annotated genes within the search window of ANY genotyped
#     SNP (default), so the test is not biased by where DArTseq markers fall; or all
#     annotated genes (FM$go_universe = "all").
#   - Test sets: all candidate genes, and the candidate genes of each environmental variable.
#   - One-sided Fisher's exact (hypergeometric) test per term, BH-adjusted within each
#     ontology (BP, MF, CC). Each term also reports how many independent genomic regions
#     its genes come from: terms supported by a single region usually reflect one cluster
#     of tandem paralogs near one marker, not a pathway-level signal, and are flagged.
#
#   Terms need >= 2 candidate genes (FM$go_min_hits). robust = FDR-significant AND from >= 2 regions.
#
# Usage:  Rscript candidate_genes_fingermillet_GO.R [input_dir] [output_dir]
#   input_dir  : folder searched (recursively) for the input files below; default "."
#   output_dir : where tables are written; default "candidate_genes_out"
# Input files are found by name pattern, so they can sit anywhere under input_dir and
# keep any download name (gzipped or not). To use a specific file, set its path in FM.
#
#   required  VCF                      *.vcf / *.vcf.gz
#   required  eGWAS results            significant_markers.csv (from environmental_gwas_complete.R)
#                                       or *EGWAS*.csv in the older format
#   required  gene annotation GFF3     *gene[._]gff3*  (not gene_exons)
#   optional  functional annotation    *annotation_info*
#   optional  gene descriptions        *defline*
#   optional  primary-transcript proteins *protein_primaryTranscriptOnly*  (BLAST only)
#   optional  LD decay summary         *LD_decay*.csv  (set = A / B / genome, half_decay_bp);
#                                       overrides the published values below
#   optional  chromosome map           *chr_map*.tsv   (2 cols: VCF name, reference name)
#   optional  GO ontology              *.obo (go-basic.obo); downloaded if absent
#
# Reference: ML-238 genome (Devos et al. 2023, Nat Commun 14:3694), Phytozome
#   Ecoracana_560_v1.1 + annotation P14. AABB allotetraploid (2n = 4x = 36);
#   one transcript per locus, CDS only (no UTRs annotated).
# Chromosome names: names already like 1A / Chr01A are used as is. Numeric TASSEL names
#   are mapped 1-18 -> 1A,1B,2A,2B,...,9A,9B (verified for finger_millet_404.vcf by matching
#   SNP and gene extents); 19 = concatenated unplaced scaffolds, reported but not searched.
#   A chr_map file overrides this. Every run writes chr_name_check.tsv to confirm the fit.
#
# DArTseq panel: Bancic et al. 2024, Plant Genome 17:e20392 (doi:10.1002/tpg2.20392).
#   8,778 SNPs (8,096 on chromosomes, 682 on scaffolds); per-chromosome SNP counts in their
#   Table 1 match the 1-18 mapping above exactly (built in as a check). The VCF has 404
#   accessions = the paper's known-origin accessions (423 minus 19 "other"); groups are
#   East Africa, Southern Africa, India, Nepal (DAPC clusters 1-4).
#   LD decay (r2 = 0.2, structure-corrected): A sub-genome 88 kb, B sub-genome 168 kb,
#   whole genome 106 kb -> used as subgenome-specific search windows. The paper listed
#   candidate genes within +/-25 kb of significant MTAs; a column flags genes inside that.
#   The paper warns that most significant SNPs had low MAF and that collapsed homeologous
#   loci show excess heterozygosity (> 2pq(1-F), F = 0.5); MAF and heterozygosity are
#   computed from the VCF and flagged for every marker.
#
# Steps: (1) map chromosomes, check positions, marker MAF/heterozygosity; (2) rank markers:
#   priority = 2 x within-group significance (only if the table has a set column)
#   + 1 if significant for >1 variable + 1 if max PVE >= FM$pve_high + 1 if min q <= FM$q_strict;
#   (3) genes within the subgenome-specific LD window of every significant marker (FM$min_priority = 0);
#   (4) SNP position (CDS / intron / upstream / downstream); (5) functional annotation;
#   (6) optional BLASTP; (7) GO enrichment.
# Outputs: candidate_markers_prioritised.csv, candidate_genes.csv, chr_name_check.tsv,
#          GO_enrichment.csv, GO_enrichment_significant.csv, GO_universe_summary.csv,
#          GO_enrichment.pdf, blast_<db>.tsv (optional)
################################################################################
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[1] else "."
out_dir <- if (length(args) >= 2) args[2] else "candidate_genes_out"

## ---- Settings ---------------------------------------------------------------------
# Leave a path NULL to auto-detect it in in_dir, or set it explicitly.
FM <- list(
  vcf = NULL, egwas = NULL, gff = NULL, annot_info = NULL, defline = NULL,
  proteins_fa = NULL, ld_summary = NULL, chr_map = NULL,
  # numeric VCF chromosome -> reference chromosome (used when no chr_map file is given)
  vcf_chr_map = setNames(paste0(rep(1:9, each = 2), c("A", "B")), 1:18),
  blast_dbs   = c(TAIR10 = "blastdb/TAIR10_pep", rice = "blastdb/rice_MSU7_pep"),  # optional
  # +/- search window = LD decay distance (r2 = 0.2) per sub-genome, Bancic et al. 2024
  ld_window   = c(A = 88000, B = 168000, genome = 106000),
  min_window  = 25000,     # never search less than the paper's +/-25 kb
  max_window  = 300000,    # cap on the window
  paper_window = 25000,    # flag genes within the paper's candidate-gene interval
  low_maf     = 0.05,      # flag markers below this MAF
  inbreeding_F = 0.5,      # F for the homeolog heterozygosity filter (paper used 0.5)
  # published SNPs per chromosome (Table 1), used to check the chromosome mapping
  paper_snps  = c(`1A` = 742, `2A` = 683, `3A` = 570, `4A` = 618, `5A` = 1057, `6A` = 357, `7A` = 553,
                  `8A` = 307, `9A` = 486, `1B` = 196, `2B` = 323, `3B` = 136, `4B` = 158, `5B` = 659,
                  `6B` = 244, `7B` = 267, `8B` = 294, `9B` = 446),
  full_set    = "full",    # name of the all-accessions set (only if the table has a set column)
  pve_high    = 5,         # priority +1 if the marker explains >= this % of a variable's variance
  q_strict    = 0.05,      # priority +1 if q <= this
  min_priority = 0,        # markers searched: priority >= this (0 = every significant marker)
  threads     = 4,
  # ---- GO enrichment ----
  go_obo      = NULL,      # NULL = *.obo in input_dir, else downloaded to the output folder
  go_obo_url  = "http://purl.obolibrary.org/obo/go/go-basic.obo",
  go_universe = "near_snps",  # "near_snps" or "all"
  go_min_term = 5,         # tested GO terms need >= this many universe genes
  go_max_term = 500,       #   ... and <= this many (very broad terms are uninformative)
  go_min_genes = 5,        # test a gene set only if it has >= this many GO-annotated genes
  go_min_hits = 2,         # test a term only if >= this many study genes carry it
  go_fdr      = 0.05,      # FDR for calling a term enriched
  go_plot_top = 10)        # terms per ontology in the figure

## ---- locate input files -------------------------------------------------------------
all_files <- list.files(in_dir, recursive = TRUE, full.names = TRUE)
find_in <- function(key, patterns, exclude = NULL, required = FALSE) {
  if (!is.null(FM[[key]])) { if (!file.exists(FM[[key]])) stop(key, " not found: ", FM[[key]]); return(FM[[key]]) }
  for (p in patterns) {
    f <- all_files[grepl(p, basename(all_files), ignore.case = TRUE)]
    if (!is.null(exclude)) f <- f[!grepl(exclude, basename(f), ignore.case = TRUE)]
    if (length(f)) {
      f <- f[order(nchar(basename(f)), f)]
      if (length(f) > 1) message("  ", key, ": ", length(f), " matches, using ", f[1])
      return(f[1]) }
  }
  if (required) stop("No ", key, " file found in ", normalizePath(in_dir), " (patterns: ", paste(patterns, collapse = " | "), ")")
  NULL }
FM$vcf         <- find_in("vcf", "\\.vcf(\\.gz)?$", required = TRUE)
FM$egwas       <- find_in("egwas", c("^significant_markers.*\\.csv$", "EGWAS.*signif.*\\.csv$", "EGWAS.*\\.csv$"), required = TRUE)
FM$gff         <- find_in("gff", c("gene[._]gff3", "\\.gff3?(\\.gz)?$"), exclude = "exons", required = TRUE)
FM$annot_info  <- find_in("annot_info", "annotation_info")
FM$defline     <- find_in("defline", "defline")
FM$proteins_fa <- find_in("proteins_fa", c("protein_primaryTranscriptOnly", "protein.*\\.fa"))
FM$ld_summary  <- find_in("ld_summary", "LD_decay.*\\.csv$")
FM$chr_map     <- find_in("chr_map", "chr_map.*\\.(tsv|txt)$")
FM$go_obo      <- find_in("go_obo", "\\.obo$")

## ---- helpers --------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tab <- function(f) file.path(out_dir, f)
res <- tab
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
has <- function(f) !is.null(f) && length(f) == 1 && file.exists(f)
cat_cmd <- function(f) paste(if (grepl("gz$", f)) "zcat" else "cat", shQuote(f))
to_ref  <- function(x) {                          # VCF/EGWAS chromosome -> reference name
  x <- sub("^(chromosome|chrom|chr)[_.-]?", "", as.character(x), ignore.case = TRUE)
  x <- sub("^0+(?=[0-9])", "", toupper(x), perl = TRUE)
  ifelse(x %in% names(chr_lookup), chr_lookup[x],
         ifelse(grepl("^[1-9][AB]$|^SCAFFOLD_", x), sub("^SCAFFOLD", "scaffold", x), NA_character_)) }
subgenome <- function(chr) ifelse(grepl("^[1-9][AB]$", chr), substring(chr, 2), NA_character_)
get_att <- function(a, key) {
  m <- regexpr(paste0("(^|;)", key, "=[^;]+"), a); out <- rep(NA_character_, length(a))
  out[m > 0] <- sub(paste0("^;?", key, "="), "", regmatches(a, m)); out }
clean_def <- function(x) {                       # "(1 of 2) PTHR123:SF4 - NAME" -> "NAME"
  x <- sub("^\\([0-9]+ of [0-9]+\\) ", "", x)
  x <- sub("^(PTHR[0-9]+(:SF[0-9]+)?|[0-9]+\\.[0-9.-]+|K[0-9]{5}|KOG[0-9]+|PF[0-9]+) - ", "", x)
  x[x == ""] <- NA; x }
db_exists <- function(db) any(file.exists(paste0(db, c(".pin", ".pal", ".psq", ".pdb", ".00.pin"))))
tbl_str <- function(x) { t <- table(x, useNA = "ifany"); paste(names(t), t, sep = "=", collapse = " ") }

chr_lookup <- FM$vcf_chr_map
if (!is.null(FM$chr_map)) {
  cm <- fread(FM$chr_map, header = FALSE, col.names = c("vcf_name", "ref_name"), colClasses = "character")
  chr_lookup <- setNames(cm$ref_name, toupper(sub("^0+(?=[0-9])", "", sub("^(chromosome|chrom|chr)[_.-]?", "", cm$vcf_name, ignore.case = TRUE), perl = TRUE)))
}
fk <- c("vcf", "egwas", "gff", "annot_info", "defline", "proteins_fa", "ld_summary", "chr_map", "go_obo")
msg("Inputs:", paste(fk, sapply(FM[fk], function(x) if (is.null(x)) "-" else basename(x)), sep = "=", collapse = "; "))

## ---- 1. read GFF and VCF; check the chromosome mapping -------------------------------
gff <- fread(cmd = paste(cat_cmd(FM$gff), "| grep -v '^#'"), sep = "\t", header = FALSE, quote = "",
             col.names = c("seqid", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
gff <- gff[type %in% c("gene", "mRNA", "CDS")]
ext <- gff[type == "gene", .(gff_max = max(end)), by = seqid]

vraw <- fread(cmd = paste(cat_cmd(FM$vcf), "| grep -v '^##'"), sep = "\t", colClasses = "character")
vcf <- data.table(CHROM = vraw[[1]], POS = as.integer(vraw[[2]]), ID = vraw[[3]])
vcf[, chr_ref := to_ref(CHROM)]
# marker MAF and observed heterozygosity from GT (first field of each sample column)
gt <- as.matrix(vraw[, -(1:9)]); rm(vraw)
gt <- gsub("\\|", "/", sub(":.*", "", gt))
dos <- matrix(NA_integer_, nrow(gt), ncol(gt))
dos[gt == "0/0"] <- 0L; dos[gt %in% c("0/1", "1/0")] <- 1L; dos[gt == "1/1"] <- 2L
n_called <- rowSums(!is.na(dos)); p_alt <- rowSums(dos, na.rm = TRUE) / (2 * n_called)
vcf[, `:=`(n_called = n_called, MAF = round(pmin(p_alt, 1 - p_alt), 4),
           het_obs = round(rowSums(dos == 1L, na.rm = TRUE) / n_called, 4))]
vcf[, het_excess := het_obs > 2 * MAF * (1 - MAF) * (1 - FM$inbreeding_F)]
msg("Genotypes:", ncol(gt), "samples; median MAF", median(vcf$MAF, na.rm = TRUE), "; median het",
    median(vcf$het_obs, na.rm = TRUE), ";", sum(vcf$het_excess, na.rm = TRUE), "SNPs with excess heterozygosity")
rm(gt, dos)
chk <- vcf[, .(n_snps = .N, vcf_max = max(POS)), by = .(CHROM, chr_ref)][ext, gff_max := i.gff_max, on = c(chr_ref = "seqid")]
chk[, paper_n_snps := FM$paper_snps[chr_ref]]
fwrite(chk[order(suppressWarnings(as.integer(CHROM)), CHROM)], tab("chr_name_check.tsv"), sep = "\t")
pc <- chk[!is.na(paper_n_snps)]
if (nrow(pc)) msg("SNP counts match Bancic et al. 2024 Table 1 on", sum(pc$n_snps == pc$paper_n_snps), "of",
                  nrow(pc), "chromosomes", if (all(pc$n_snps == pc$paper_n_snps)) "(mapping confirmed)" else "- check mapping")
msg("VCF:", nrow(vcf), "SNPs;", sum(is.na(vcf$chr_ref)), "on unplaced chromosome(s)",
    paste(unique(vcf$CHROM[is.na(vcf$chr_ref)]), collapse = ","), "(not searched)")
bad <- chk[!is.na(chr_ref) & (is.na(gff_max) | vcf_max > gff_max + 2e6)]
if (nrow(bad)) warning("VCF positions exceed the mapped reference chromosome for: ",
                       paste(bad$CHROM, collapse = ", "), " - check the chromosome map / reference version")

## ---- 2. prioritise markers ---------------------------------------------------------
hits <- read.csv(FM$egwas, stringsAsFactors = FALSE, check.names = FALSE)
# map environmental_gwas_complete.R column names onto the names used here
ren_in <- c(marker = "SNP", chrom = "CHROM", pos = "POS", variable = "trait", q_value = "q", PVE_percent = "PVE")
for (o in names(ren_in)) if (o %in% names(hits) && !ren_in[[o]] %in% names(hits)) names(hits)[names(hits) == o] <- ren_in[[o]]
if (!all(c("CHROM", "POS") %in% names(hits))) {          # fill positions from VCF IDs
  hits <- merge(hits, as.data.frame(vcf[, .(SNP = ID, CHROM, POS)]), by = "SNP")
  msg("CHROM/POS taken from VCF IDs for", nrow(hits), "rows")
}
if (length(miss <- setdiff(c("SNP", "CHROM", "POS", "trait"), names(hits))))
  stop("eGWAS table lacks columns: ", paste(miss, collapse = ", "))
HAS_SET <- "set" %in% names(hits)
if (!HAS_SET) hits$set <- FM$full_set
for (cc in c("PVE", "q", "beta")) if (!cc %in% names(hits)) hits[[cc]] <- NA_real_
if (!"model" %in% names(hits)) hits$model <- NA_character_
if (!"effect_direction" %in% names(hits)) hits$effect_direction <- NA_character_
hits$CHROM <- as.character(hits$CHROM)
msg("eGWAS hits:", nrow(hits), "marker-variable associations;", n_distinct(hits$SNP), "markers;",
    n_distinct(hits$trait), "variables", if (!HAS_SET) "(no set column: within-group term not used)" else "")

ldw <- FM$ld_window
if (has(FM$ld_summary)) {                         # own LD estimates override the published ones
  ld <- read.csv(FM$ld_summary); ld <- ld[ld$set %in% names(ldw), ]
  ldw[ld$set] <- ld$half_decay_bp; msg("LD windows from", basename(FM$ld_summary))
}
ldw <- pmin(pmax(ldw, FM$min_window), FM$max_window)
win_for <- function(sg) unname(ifelse(!is.na(sg) & sg %in% names(ldw), ldw[sg], ldw["genome"]))
msg("Search windows (+/- bp):", paste(names(ldw), ldw, sep = "=", collapse = " "))

mmax <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
mmin <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
mk <- hits %>% group_by(SNP, CHROM, POS) %>%
  summarise(n_vars = n_distinct(trait[set == FM$full_set]), within_group = any(set != FM$full_set),
            sets = paste(sort(unique(set)), collapse = ";"), traits = paste(sort(unique(trait)), collapse = ";"),
            models = paste(sort(unique(na.omit(model))), collapse = ";"),
            effect = paste(unique(na.omit(paste0(trait, ": ", effect_direction)[!is.na(effect_direction)])), collapse = "; "),
            max_PVE = mmax(PVE), min_q = mmin(q), .groups = "drop") %>%
  mutate(chr_ref = to_ref(CHROM), subgenome = subgenome(chr_ref), window_bp = win_for(subgenome),
         priority = within_group * 2 + (n_vars > 1) + coalesce(max_PVE >= FM$pve_high, FALSE) +
                    coalesce(min_q <= FM$q_strict, FALSE)) %>%
  left_join(as.data.frame(vcf[!duplicated(ID), .(SNP = ID, MAF, het_obs, het_excess)]), by = "SNP") %>%
  mutate(low_MAF = MAF < FM$low_maf) %>%
  arrange(desc(priority), min_q) %>% relocate(chr_ref, subgenome, .after = POS)
if (all(is.na(mk$MAF))) { msg("EGWAS SNP names do not match VCF IDs: MAF/heterozygosity not added")
} else msg("Markers with low MAF (<", FM$low_maf, "):", sum(mk$low_MAF, na.rm = TRUE), "; with excess heterozygosity:",
         sum(mk$het_excess, na.rm = TRUE))
write.csv(mk, tab("candidate_markers_prioritised.csv"), row.names = FALSE)
msg(nrow(mk), "markers; by subgenome:", tbl_str(mk$subgenome))

## ---- 3. one row per locus with its primary transcript ----------------------------------
tx <- gff[type == "mRNA"]
tx[, `:=`(tx_ID = get_att(attr, "ID"), transcript = get_att(attr, "Name"), Parent = get_att(attr, "Parent"),
          rank = as.integer(!(get_att(attr, "longest") %in% "1")))]
setorder(tx, Parent, rank, start)
ptx <- tx[!duplicated(Parent), .(Parent, tx_ID, transcript)]
genes <- gff[type == "gene"][, `:=`(ID = get_att(attr, "ID"), locus = get_att(attr, "Name"))]
genes <- merge(genes, ptx, by.x = "ID", by.y = "Parent", all.x = TRUE)[, attr := NULL]
cds <- gff[type == "CDS", .(tx_ID = get_att(attr, "Parent"), start, end)]
setkey(cds, tx_ID)
msg("GFF:", nrow(genes), "loci on", uniqueN(genes$seqid), "sequences;", nrow(cds), "CDS segments")

# functional annotation (annotation_info: one row per transcript; defline: description)
if (has(FM$annot_info)) {
  an <- fread(cmd = cat_cmd(FM$annot_info), sep = "\t", quote = "", colClasses = "character")
  setnames(an, sub("^#", "", names(an)))
  ren  <- c(Pfam = "Pfam", Panther = "PANTHER", KOG = "KOG", KO = "KO", ec = "EC", GO = "GO",
            `Best-hit-arabi-name` = "arabidopsis_hit", `Best-hit-arabi-defline` = "arabidopsis_desc",
            `Best-hit-rice-name` = "rice_hit", `Best-hit-rice-defline` = "rice_desc")
  keep <- intersect(names(ren), names(an))
  an <- an[, c("transcriptName", keep), with = FALSE]; setnames(an, keep, ren[keep])
  for (cc in c("arabidopsis_desc", "rice_desc")) if (cc %in% names(an)) an[, (cc) := clean_def(get(cc))]
  genes <- merge(genes, an[!duplicated(transcriptName)], by.x = "transcript", by.y = "transcriptName", all.x = TRUE)
}
if (has(FM$defline)) {
  dl <- fread(cmd = cat_cmd(FM$defline), sep = "\t", header = FALSE, quote = "",
              col.names = c("transcript", "def_type", "description"))
  dl[, description := clean_def(description)]
  genes <- merge(genes, dl[!duplicated(transcript), .(transcript, description)], by = "transcript", all.x = TRUE)
}
ann_present <- intersect(c("description", "arabidopsis_hit", "arabidopsis_desc", "rice_hit", "rice_desc",
                           "PANTHER", "Pfam", "GO", "KOG", "KO", "EC"), names(genes))
msg("Annotated loci:", if ("description" %in% names(genes)) sum(!is.na(genes$description)) else 0, "with a description")

## ---- 4. genes within +/- window of priority markers --------------------------------
snp_context <- function(p, tx_id, strand, s, e) {
  if (p < s) return(if (strand == "-") "downstream" else "upstream")
  if (p > e) return(if (strand == "-") "upstream" else "downstream")
  cc <- cds[.(tx_id), nomatch = NULL]
  if (!nrow(cc)) return("genic")
  if (any(p >= cc$start & p <= cc$end)) "CDS"
  else if (p < min(cc$start) || p > max(cc$end)) "genic_noncoding" else "intron" }

top <- mk %>% filter(priority >= FM$min_priority, !is.na(chr_ref))
if (!HAS_SET) top$sets <- NULL
if (!nrow(top)) top <- head(filter(mk, !is.na(chr_ref)), 20)
cg <- rbindlist(lapply(seq_len(nrow(top)), function(i) {
  p <- top$POS[i]; win <- top$window_bp[i]
  g <- copy(genes[seqid == top$chr_ref[i] & end >= p - win & start <= p + win])
  if (!nrow(g)) return(NULL)
  g[, `:=`(SNP = top$SNP[i], CHROM_vcf = top$CHROM[i], POS = p, traits = top$traits[i],
           sets = if (HAS_SET) top$sets[i] else NA_character_, priority = top$priority[i], max_PVE = top$max_PVE[i], min_q = top$min_q[i],
           effect = top$effect[i], MAF = top$MAF[i], low_MAF = top$low_MAF[i],
           het_excess = top$het_excess[i], window_bp = win, dist_bp = pmax(0L, start - p, p - end))]
  g[, in_paper_window := dist_bp <= FM$paper_window]
  g[, snp_position := mapply(snp_context, p, tx_ID, strand, start, end)]
  g[, nearest := dist_bp == min(dist_bp)]
  g[, c("SNP", "CHROM_vcf", "POS", "traits", "sets", "priority", "max_PVE", "min_q", "effect", "seqid", "start", "end", "strand",
        "locus", "transcript", "MAF", "low_MAF", "het_excess", "window_bp", "dist_bp", "snp_position",
        "nearest", "in_paper_window", ann_present), with = FALSE] }), fill = TRUE)
if (nrow(cg)) setnames(cg, "seqid", "chr_ref")[, subgenome := subgenome(chr_ref)]
msg("Markers searched:", nrow(top), "; candidate loci:", nrow(cg), "(", sum(cg$in_paper_window), "within +/-",
    FM$paper_window, "bp)",
    "; markers with no gene in window:", nrow(top) - uniqueN(cg$SNP))

## ---- 5. optional BLASTP (primary transcripts) --------------------------------------------
if (nrow(cg) && has(FM$proteins_fa) && nzchar(Sys.which("blastp")) && any(sapply(FM$blast_dbs, db_exists))) {
  fa  <- readLines(if (grepl("gz$", FM$proteins_fa)) gzfile(FM$proteins_fa) else FM$proteins_fa)
  hdr <- grepl("^>", fa); rec <- cumsum(hdr)
  keep <- sub("^>(\\S+).*", "\\1", fa[hdr]) %in% cg$transcript
  sel <- fa[keep[rec]]; sel[!grepl("^>", sel)] <- gsub("\\*", "", sel[!grepl("^>", sel)])
  writeLines(sel, res("candidate_proteins.fa"))
  msg("Proteins extracted:", sum(keep), "of", uniqueN(cg$transcript))
  for (db in names(FM$blast_dbs)) {
    if (!db_exists(FM$blast_dbs[[db]])) { msg("BLAST db not found, skipped:", FM$blast_dbs[[db]]); next }
    out <- tab(paste0("blast_", db, ".tsv"))
    system(sprintf("blastp -query %s -db %s -evalue 1e-10 -max_target_seqs 3 -num_threads %d -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' > %s",
                   shQuote(res("candidate_proteins.fa")), shQuote(FM$blast_dbs[[db]]), FM$threads, shQuote(out)))
    if (file.exists(out) && file.size(out) > 0) {
      b <- fread(out, sep = "\t", header = FALSE, quote = "",
                 col.names = c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle"))
      b <- b[order(qseqid, -bitscore)][!duplicated(qseqid)]
      nm <- paste0("blast_", db, c("_hit", "_pident", "_evalue", "_desc"))
      setnames(b, c("sseqid", "pident", "evalue", "stitle"), nm)
      cg <- merge(cg, b[, c("qseqid", nm), with = FALSE], by.x = "transcript", by.y = "qseqid", all.x = TRUE, sort = FALSE)
      msg("BLAST", db, ":", nrow(b), "queries with hits ->", out)
    }
  }
}

if (nrow(cg)) {
  setorder(cg, -priority, SNP, dist_bp)
  setcolorder(cg, c("SNP", "CHROM_vcf", "chr_ref", "subgenome", "POS", "traits", "sets", "priority", "max_PVE", "min_q", "effect",
                    "MAF", "low_MAF", "het_excess", "locus", "transcript", "start", "end", "strand", "window_bp",
                    "dist_bp", "snp_position", "nearest", "in_paper_window"))
}
write.csv(cg, tab("candidate_genes.csv"), row.names = FALSE)
msg("Candidate loci by subgenome:", tbl_str(cg$subgenome), "| SNP position:", tbl_str(cg$snp_position))

## ---- 7. GO enrichment ----------------------------------------------------------------
if (!nrow(cg)) {
  msg("No candidate genes: GO enrichment skipped")
} else if (!"GO" %in% names(genes) || all(is.na(genes$GO) | genes$GO == "")) {
  msg("No GO annotations (annotation_info file with a GO column needed): GO enrichment skipped")
} else {
  # ontology: local *.obo, else download go-basic.obo once
  obo <- FM$go_obo
  if (is.null(obo)) {
    obo <- tab("go-basic.obo")
    if (!file.exists(obo) || file.size(obo) < 1e6) {
      msg("Downloading GO ontology:", FM$go_obo_url)
      ok <- tryCatch({ download.file(FM$go_obo_url, obo, quiet = TRUE, mode = "wb"); TRUE }, error = function(e) FALSE)
      if (!ok || !file.exists(obo) || file.size(obo) < 1e6) {
        obo <- NULL; msg("GO ontology could not be downloaded; put go-basic.obo in the input folder. GO enrichment skipped")
      }
    }
  }
}
if (nrow(cg) && exists("obo") && !is.null(obo)) {
  msg("GO ontology:", basename(obo))
  # ---- parse OBO ([Term] stanzas: id, name, namespace, alt_id, is_a, part_of, is_obsolete)
  L <- readLines(obo, warn = FALSE)
  is_hdr <- grepl("^\\[", L); blk <- cumsum(is_hdr)
  term_blk <- which(L[is_hdr] == "[Term]")
  L <- L[blk %in% term_blk & !is_hdr]; b <- blk[blk %in% term_blk & !is_hdr]
  tg <- sub(":.*$", "", L); vl <- sub("^[^:]+: ", "", L)
  od <- data.table(b, tg, vl)[tg %in% c("id", "name", "namespace", "alt_id", "is_a", "relationship", "is_obsolete")]
  terms <- dcast(od[tg %in% c("id", "name", "namespace", "is_obsolete")], b ~ tg, value.var = "vl", fun.aggregate = function(x) x[1])
  if (!"is_obsolete" %in% names(terms)) terms[, is_obsolete := NA]
  terms <- terms[is.na(is_obsolete) | is_obsolete != "true"]
  terms[, ontology := c(biological_process = "BP", molecular_function = "MF", cellular_component = "CC")[namespace]]
  ids <- od[tg == "id", .(b, id = vl)]
  alt <- merge(od[tg == "alt_id", .(b, alt = vl)], ids, by = "b")
  par <- od[tg == "is_a" | (tg == "relationship" & grepl("^part_of ", vl))]
  par[, parent := ifelse(grepl("GO:[0-9]{7}", vl), sub("^.*?(GO:[0-9]{7}).*$", "\\1", vl, perl = TRUE), NA_character_)]
  par <- merge(par[!is.na(parent), .(b, parent)], ids, by = "b")[id %in% terms$id & parent %in% terms$id]
  ns <- setNames(terms$ontology, terms$id); nm <- setNames(terms$name, terms$id)
  par <- par[ns[id] == ns[parent]]                                  # stay within one ontology
  par_list <- split(par$parent, par$id)
  anc_env <- new.env(hash = TRUE)
  get_anc <- function(t) {
    r <- anc_env[[t]]; if (!is.null(r)) return(r)
    p <- par_list[[t]]; r <- if (is.null(p)) character(0) else unique(c(p, unlist(lapply(p, get_anc))))
    assign(t, r, envir = anc_env); r }
  roots <- c("GO:0008150", "GO:0003674", "GO:0005575")
  msg("GO terms:", nrow(terms), "| parent links:", nrow(par))

  # ---- gene -> GO (direct), alt ids mapped, then propagated to ancestors
  gd <- genes[!is.na(GO) & GO != "", .(term = unlist(strsplit(GO, "[ ,;|]+"))), by = locus][grepl("^GO:", term)]
  am <- setNames(alt$id, alt$alt); gd[term %in% names(am), term := am[term]]
  n_unknown <- uniqueN(gd$term[!gd$term %in% terms$id])
  gd <- unique(gd[term %in% terms$id])
  gp <- gd[, .(term = unique(c(term, unlist(lapply(unique(term), get_anc))))), by = locus][!term %in% roots]
  gp[, ontology := ns[term]]
  msg("Genes with GO terms:", uniqueN(gd$locus), "| direct annotations:", nrow(gd), "| after propagation:", nrow(gp),
      if (n_unknown) paste("|", n_unknown, "GO IDs not in this ontology version (ignored)") else "")

  # GO names on the candidate gene table (direct annotations)
  gnames <- gd[, .(GO_names = paste(sprintf("%s %s", term, nm[term]), collapse = "; ")), by = locus]
  cg <- merge(cg, gnames, by = "locus", all.x = TRUE, sort = FALSE)
  setorder(cg, -priority, SNP, dist_bp)
  write.csv(cg, tab("candidate_genes.csv"), row.names = FALSE)

  # ---- universe
  if (FM$go_universe == "near_snps") {
    sw <- vcf[!is.na(chr_ref), .(chr = chr_ref, POS, w = win_for(subgenome(chr_ref)))]
    sw <- sw[, .(chr, s = pmax(1L, POS - w), e = POS + w)]; setkey(sw, chr, s, e)
    gi <- genes[, .(locus, chr = seqid, s = start, e = end)]
    uni <- unique(foverlaps(gi, sw, by.x = c("chr", "s", "e"), type = "any", nomatch = NULL)$locus)
    uni_lab <- "annotated genes within the search window of any genotyped SNP"
  } else { uni <- genes$locus; uni_lab <- "all annotated genes" }
  uni <- intersect(uni, gp$locus)
  msg("GO universe:", length(uni), uni_lab)

  # ---- independent regions: overlapping marker windows merged
  rg <- as.data.table(top)[, .(SNP, chr = chr_ref, s = POS - window_bp, e = POS + window_bp)][order(chr, s)]
  rg[, region := paste0(chr, "_", cumsum(c(1L, as.integer(s[-1] > head(cummax(e), -1))))), by = chr]
  gene_region <- unique(merge(cg[, .(locus, SNP)], rg[, .(SNP, region)], by = "SNP")[, .(locus, region)])

  # ---- study sets: all candidates + candidates of each variable
  cgt <- cg[, .(trait = unlist(strsplit(traits, ";"))), by = .(locus, SNP)]
  study_sets <- c(list(`all candidate genes` = unique(cg$locus)), split(cgt$locus, cgt$trait))
  study_sets <- lapply(study_sets, unique)

  U <- gp[locus %in% uni]; N <- length(uni)
  Nt <- U[, .(N_t = .N), by = .(term, ontology)][N_t >= FM$go_min_term & N_t <= FM$go_max_term]
  enrich <- function(study, label) {
    st <- intersect(study, uni); n <- length(st)
    if (n < FM$go_min_genes) { msg("  ", label, ": only", n, "GO-annotated genes - not tested"); return(NULL) }
    S <- U[locus %in% st, .(k = .N, genes = paste(sort(locus), collapse = ";")), by = term]
    r <- merge(Nt, S, by = "term")[k >= FM$go_min_hits]
    if (!nrow(r)) { msg("  ", label, ":", n, "genes; no term carried by >=", FM$go_min_hits, "genes"); return(NULL) }
    r[, `:=`(gene_set = label, n_study = n, n_universe = N, expected = n * N_t / N)]
    r[, `:=`(fold_enrichment = k / expected, p_value = phyper(k - 1, N_t, N - N_t, n, lower.tail = FALSE))]
    r[, q_value := p.adjust(p_value, "BH"), by = ontology]
    r[, n_regions := vapply(strsplit(genes, ";"), function(g) uniqueN(gene_region$region[gene_region$locus %in% g]), 0L)]
    r[, `:=`(term_name = nm[term], single_region = n_regions < 2,
             significant = q_value <= FM$go_fdr)]
    r[, robust := significant & !single_region]
    msg("  ", label, ":", n, "genes;", nrow(r), "terms tested;", sum(r$significant), "enriched at FDR", FM$go_fdr,
        if (sum(r$significant)) paste0("(", sum(r$significant & !r$single_region), " from >1 region)") else "")
    r }
  msg("GO enrichment (one-sided Fisher, BH within ontology):")
  go_res <- rbindlist(Map(enrich, study_sets, names(study_sets)), fill = TRUE)
  if (nrow(go_res)) {
    setcolorder(go_res, c("gene_set", "ontology", "term", "term_name", "k", "n_study", "N_t", "n_universe",
                          "expected", "fold_enrichment", "p_value", "q_value", "significant", "n_regions",
                          "single_region", "robust", "genes"))
    setnames(go_res, c("term", "k", "N_t"), c("GO_ID", "n_genes_in_set", "n_genes_in_universe"))
    go_res <- go_res[order(gene_set != "all candidate genes", gene_set, ontology, p_value)]
    fwrite(go_res, tab("GO_enrichment.csv"))
    fwrite(go_res[significant == TRUE], tab("GO_enrichment_significant.csv"))
    fwrite(data.table(universe = uni_lab, n_universe_genes = N, n_terms_tested_size_range = nrow(Nt),
                      min_term_size = FM$go_min_term, max_term_size = FM$go_max_term, ontology_file = basename(obo)),
           tab("GO_universe_summary.csv"))

    # ---- figure: top terms per ontology for each gene set (one page per set)
    dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    dev(tab("GO_enrichment.pdf"), width = 10, height = 8, onefile = TRUE)
    for (gs in unique(go_res$gene_set)) {
      d <- go_res[gene_set == gs][order(p_value)][, head(.SD, FM$go_plot_top), by = ontology]
      if (!nrow(d)) next
      d[, label := paste0(substr(term_name, 1, 55), ifelse(nchar(term_name) > 55, "...", ""), " (", GO_ID, ")")]
      d[, label := factor(label, levels = rev(unique(label)))]
      p <- ggplot(d, aes(fold_enrichment, label)) +
        geom_segment(aes(x = 0, xend = fold_enrichment, yend = label), colour = "grey80") +
        geom_point(aes(size = n_genes_in_set, colour = -log10(p_value), shape = single_region)) +
        scale_colour_viridis_c(option = "plasma", name = expression(-log[10](italic(p)))) +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), name = "Single region") +
        scale_size_continuous(range = c(2, 7), name = "Genes") +
        geom_text(data = d[significant == TRUE], aes(label = "*"), nudge_x = 0, vjust = -0.6, size = 5) +
        facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
        labs(x = "Fold enrichment", y = NULL, title = paste("GO enrichment:", gs),
             subtitle = paste0("Top ", FM$go_plot_top, " terms per ontology by p; * = FDR <= ", FM$go_fdr,
                               "; open circle = genes from a single genomic region\nUniverse: ", N, " ", uni_lab)) +
        theme_bw(base_size = 9) + theme(strip.text.y = element_text(angle = 0, face = "bold"))
      print(p)
    }
    invisible(dev.off())
    msg("GO results ->", tab("GO_enrichment.csv"), "and GO_enrichment.pdf")
  } else msg("No GO terms could be tested")
}
