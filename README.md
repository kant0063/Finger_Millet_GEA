# Finger millet GEA — centroid-filtered analysis

This is a re-analysis of the finger millet genome–environment association workflow, keeping **one accession per putative country-centroid coordinate**.

In the full panel, 157 of 403 georeferenced accessions sit on three points that match geocoded country names rather than collection sites:

| Point | Coordinates | Accessions |
|---|---|---|
| Nepal | 28.108°N, 84.092°E | 57 |
| India | 22.351°N, 78.668°E | 52 |
| Kenya | 1.442°N, 38.431°E | 48 |

Environmental values extracted at these points do not describe where the accessions were collected. Script `00` keeps the accession with the fewest missing genotype calls at each point, removes the other 154, and writes filtered copies of every input file to `data_filtered/`. Scripts `00.5`–`09` then repeat the analyses on the filtered set (250 accessions).

## Run order

| Step | Script | Analysis | Reads | Main outputs (in `results/`) |
|---|---|---|---|---|
| 00 | `00_filter_centroid_accessions.R` | Collapse centroid coordinates to one accession each; filter all inputs | original VCF and tables | `data_filtered/` (filtered VCF and tables, kept/removed lists, centroid summary, shared-coordinate report) |
| 00.5 | `00.5_LD_decay.R` | LD decay with the four LDcorSV measures (r², r²<sub>S</sub>, r²<sub>V</sub>, r²<sub>VS</sub>) and the Hill & Weir fit; distance at r² = 0.2 by chromosome, sub-genome and genome | filtered VCF | `reference_data/LD_decay_summary.csv` (read by 06) and `00.5_LD_decay/` (summary, binned means, PDF/PNG figures) |
| 01 | `01_core_collection.R` | Core collection of 100 accessions (Core Hunter, entry-to-nearest-entry, Modified Rogers' distance) | filtered VCF | `01_core_collection/millet_core_n100.csv` |
| 02 | `02_PCA_analysis.Rmd` | PCA on LD-pruned SNPs (SNPRelate) | filtered VCF, `metadata_for_PCA.csv` | `02_PCA/` scores, variance explained, PDF/PNG figures |
| 03 | `03_TESS3_ancestry.R` | TESS3 K = 1–8 with cross-validation, Q matrix (K = 4), kriged maps, outlier loci | `test_for_numeric.csv`, `millet_passport_test.csv` | `03_TESS3/tables`, `03_TESS3/figures` |
| 03.5 | `03.5_diversity_rTASSEL.R` | Sliding-window nucleotide diversity (π per bp), Watterson's θ and Tajima's D with rTASSEL (TASSEL `seqDiversity`), one chromosome at a time (windows of 20 SNPs, step 10; unplaced scaffolds skipped) | filtered VCF | `03.5_diversity/` (window table `finger_millet_diversity.txt`, chromosome and sub-genome summaries, PDF/PNG figures of π and Tajima's D) |
| 04 | `04_environmental_GWAS.R` | Environmental GWAS, six mixed models (significance from the model with λ closest to 1), PVE per marker, significant markers clumped into independent loci, joint variance of lead markers (adjusted and cross-validated R²) | filtered VCF, environment, passport | `04_environmental_GWAS/output/significant_markers.csv` (with `locus`, `lead_marker`), `variance_explained.csv`, figures, HTML report |
| 05 | `05_pRDA_FST_LFMM2.R` | Partial RDA and variance partitioning, Weir & Cockerham F<sub>ST</sub> outliers, LFMM2 (full set, Africa, Asia) | filtered VCF, passport, environment, step-04 hits | `05_pRDA_FST_LFMM2/tables`, `figures` |
| 06 | `06_candidate_genes_GO.R` | Candidate genes in sub-genome LD windows around eGWAS hits; GO enrichment | filtered VCF, step-04 hits, `reference_data/` | `06_candidate_genes_GO/` |
| 07 | `07_genomic_prediction.R` | Environmental genomic prediction (whole genome, A/B sub-genomes; 8 models), core training sets, within-population GEAVs, candidate parents for every variable (ranked by GEAV, or by observed environment of origin where not predictable) | filtered VCF, coordinates, environment, step-01 core | `07_genomic_prediction/tables`, `figures` |
| 08 | `08_genomic_offset.R` | Genomic offset (genetic gap), SSP5-8.5, 2081–2100, mean of six CMIP6 models, ExDet masking | filtered VCF, passport, WorldClim/CMIP6 | `08_genomic_offset/` maps, GeoTIFFs, per-accession offsets |
| 09 | `09_koppen_geiger.R` | Köppen–Geiger class change at collection sites, 1991–2020 → 2041–2070 | passport, Köppen–Geiger rasters | `09_koppen_geiger/` tables and maps |

Run `00` first, then `03` (TESS3), then `00.5`: with the TESS3 ancestry matrix available, `00.5` uses it for its structure correction instead of running sNMF. Scripts `01`, `02`, `04`, `08` and `09` can then run in any order. Run `05` and `06` after `04`, `06` also after `00.5`, and `07` after `01`.

```bash
Rscript 00_filter_centroid_accessions.R
Rscript 03.5_diversity_rTASSEL.R          # nucleotide diversity and Tajima's D (rTASSEL)
Rscript 03_TESS3_ancestry.R               # optional here: lets 00.5 reuse its ancestry matrix
Rscript 00.5_LD_decay.R                   # progress log: results/00.5_LD_decay/LD_decay_progress.log
                                          # (03.5 can run any time after 00)
Rscript 01_core_collection.R
Rscript -e 'rmarkdown::render("02_PCA_analysis.Rmd")'
Rscript 04_environmental_GWAS.R
Rscript 05_pRDA_FST_LFMM2.R
Rscript 06_candidate_genes_GO.R            # optional args: <reference_data folder> <output folder>
Rscript 07_genomic_prediction.R
Rscript 08_genomic_offset.R
Rscript 09_koppen_geiger.R
```

All scripts are run from the project folder. Paths and settings are in the configuration block at the top of each script.

## Files needed

### 1. Genotype and accession files (project folder; filtered by step 00)

| File | Contents | Used by |
|---|---|---|
| `finger_millet_404.vcf` | DArTseq SNPs, 404 accessions × 8,778 SNPs (TASSEL export; chromosomes 1–19) | 00 → filtered VCF for 01, 02, 04–08 |
| `millet_passport_test.csv` | `SampleID, Longitude, Latitude, Population` (genetic cluster 1–4), optionally `Region, Country` | 00 → 03, 04, 05, 08, 09 |
| `sample_coordinates.csv` | `SampleID, Latitude, Longitude, Population`, optionally `Region, Country` | 00 → 07 |
| `millet_extracted_bioclim_values.csv` | `SampleID` + environmental variables (WorldClim bioclim and monthly temperatures, SoilGrids) at each site | 00 → 04, 05, 07 |
| `test_for_numeric.csv` | Accession name + SNPs coded 0/1/2 (same accession order as the passport file) | 00 → 03 |
| `metadata_for_PCA.csv` | `SampleID, Cluster` | 00 → 02 |
| `millet_core_n100.csv` | *Optional.* Previous core; replaced by the output of step 01 | 00 |

### 2. Reference genome and annotation (`reference_data/`; used by 06)

| File | Source |
|---|---|
| `Ecoracana_560_v1.1.gene.gff3.gz` (or `…_gene_gff3.gz`) | Phytozome, *E. coracana* v1.1 |
| `Ecoracana_560_v1.1.P14.annotation_info.txt.gz` | Phytozome, annotation release P14 |
| `Ecoracana_560_v1.1.P14.defline.txt.gz` | Phytozome, annotation release P14 |
| `Ecoracana_560_v1.1.protein_primaryTranscriptOnly.fa.gz` | *Optional.* Only needed for BLAST |
| `LD_decay_summary.csv` | **Created by step 00.5** from the filtered VCF; no need to supply it |
| `go-basic.obo` | *Optional.* Downloaded automatically from geneontology.org on the first run |

Script 06 finds these files by name pattern anywhere under `reference_data/`.

### 3. Climate data

| File / folder | Used by | Notes |
|---|---|---|
| `koppen_geiger/1991_2020/koppen_geiger_0p00833333.tif` | 09 | Beck et al. (2023), 1-km Köppen–Geiger maps |
| `koppen_geiger/2041_2070/ssp585/koppen_geiger_0p00833333.tif` | 09 | Same source |
| WorldClim 2.1 bioclim, 30 s (`wc2.1_30s_bio/*.tif`) | 03 | *Optional.* Only for re-extracting bioclim values; set `bioclim_dir` |
| `cmip6_cache/` | 08 | Created automatically: WorldClim 2.1 and CMIP6 layers are downloaded with geodata on the first run (needs internet) |

### 4. Files produced and passed between scripts

| File | From | To |
|---|---|---|
| `data_filtered/finger_millet_filtered.vcf` and filtered tables | 00 | 00.5, 01–09 |
| `reference_data/LD_decay_summary.csv` | 00.5 | 06 |
| `results/01_core_collection/millet_core_n100.csv` | 01 | 07 |
| `results/04_environmental_GWAS/output/significant_markers.csv` | 04 | 05, 06 |

## Expected folder layout

```
project/
├── 00_filter_centroid_accessions.R, 00.5_LD_decay.R, 01_… 03_TESS3_ancestry.R, 03.5_diversity_rTASSEL.R, 04_… 09_koppen_geiger.R, README.md
├── finger_millet_404.vcf
├── millet_passport_test.csv, sample_coordinates.csv, millet_extracted_bioclim_values.csv
├── test_for_numeric.csv, metadata_for_PCA.csv, [millet_core_n100.csv]
├── reference_data/   (Phytozome files, [go-basic.obo]; LD_decay_summary.csv written by 00.5)
├── koppen_geiger/    (1991_2020/…, 2041_2070/ssp585/…)
├── data_filtered/    (created by 00)
└── results/          (created by 00.5–09)
```

## R packages

R ≥ 4.2. Install the CRAN, Bioconductor and GitHub packages with:

```r
install.packages(c("data.table", "dplyr", "tidyr", "purrr", "tibble", "readr", "ggplot2", "patchwork",
                   "maps", "terra", "geodata", "rrBLUP", "BGLR", "ranger", "vegan", "corehunter",
                   "vcfR", "fields", "rworldmap", "rmarkdown"))
install.packages("BiocManager")
BiocManager::install(c("LEA", "SNPRelate", "qvalue"))
remotes::install_github("maize-genetics/rTASSEL")   # step 03.5
remotes::install_github("bcm-uga/TESS3_encho_sen")   # tess3r
```

Other requirements:
- `corehunter` and `rTASSEL` (step 03.5) need Java (a JDK; check with `java -version`).
- **rTASSEL first run:** script 03.5 calls `rTASSEL::setupTASSEL()` to download TASSEL (about 70 MB) from Maven Central, once. TASSEL is cached (on macOS in `~/Library/Caches/org.R-project.R/R/rTASSEL/`); **restart R after this first download** and run 03.5 again.
- Script 06 uses the Unix tools `gzip`, `awk`, `grep` and `cut`. These are standard on Linux and macOS; on Windows, run it under WSL.

## Things to know about the filtered set

- **Nepal is reduced to one accession.** All 57 Nepal accessions shared one coordinate. Nepal can therefore no longer be analysed as a region in any environment-based analysis.
- **Genetic cluster 2 is reduced to 16 accessions,** most of them Indian. Expect wide confidence intervals for this cluster. It may also have too few non-admixed accessions for pairwise F<sub>ST</sub> (minimum 5).
- **Other shared coordinates are reported, not removed.** Two other points are shared by ≥ 10 accessions: 12 at −19.189°S, 29.939°E (Zimbabwe) and 11 at 15.924°N, 80.186°E (India). They are listed in `data_filtered/shared_coordinates_report.csv`. They may be genuine sites, such as a market or research station; set `collapse_flagged = TRUE` in script 00 to collapse them too.
- **LD is re-estimated on the filtered set** (step 00.5) and the candidate-gene windows use those values. Run on the 250 filtered accessions, r²<sub>VS</sub> = 0.2 is reached at about 121 kb genome-wide (A ≈ 100 kb, B ≈ 190 kb), longer than the full-panel 78.7 / 64.7 / 127.8 kb. Run on the full 404-accession VCF, the same script gives 82.4 / 67.6 / 134.1 kb, so most of the increase comes from the change in panel rather than the method. To keep the full-panel windows instead, put the original `LD_decay_summary.csv` in `reference_data/` and skip 00.5.
- **Joint variance of eGWAS markers (step 04).** Fitting all significant markers of a variable jointly overfits (in-sample values above 100%). Significant markers are therefore clumped into independent loci (a marker joins a stronger marker's locus if it lies within the sub-genome LD window, set in `clump_window`, or has r² ≥ 0.5 with it), and one lead marker per locus (at most one per ten accessions) is fitted. `variance_explained.csv` reports the adjusted and 10-fold cross-validated R² of the lead markers, as a percentage of the variance remaining after structure and kinship.
- **Candidate parents (step 07).** `unpredictable_rule = "observed"` ranks variables that are not predictable within populations by the observed environment of origin (centred within population) instead of skipping them; `rank_all_variables = TRUE` adds both extremes of every variable. Each row states `rank_basis` (GEAV or observed environment) and `confidence`; report the two kinds separately.
- **Genomic offset (step 08).** Offsets are the mean across the climate models that downloaded. If the geodata server is temporarily unavailable, a model is skipped with a "failed" line in the console, and `offset_extrapolation_and_grid_summary.csv` shows the number used (`gcms`). In the current run BCC-CSM2-MR failed, so five models were averaged; rerun step 08 later to include all six (the other layers are cached).
- **Nucleotide diversity (step 03.5).** DArTseq data contain only variable sites, so π per bp is a relative measure for comparing windows, chromosomes and sub-genomes, not an absolute nucleotide diversity. Windows are counted in SNPs, not base pairs.
- **Filtered and full-panel scripts are the same analyses.** Each numbered script matches its full-panel counterpart line for line, except for input and output paths and:
  - script 04: LD windows used for clumping (filtered panel 100.2 / 190.2 / 121.2 kb; full panel 64.7 / 127.8 / 78.7 kb) and `sig_model = "best"`;
  - script 05: Africa/Asia subsets defined by genetic cluster (`subset_col = "Population"`) and the eGWAS hits read from step 04;
  - script 07: training core taken from step 01.

