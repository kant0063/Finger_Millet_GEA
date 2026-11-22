# install.packages("devtools")
devtools::install_github("bcm-uga/TESS3_encho_sen")


library(tess3r)

genotype = read.csv("test_for_numeric.csv")
coordinates = read.csv("millet_passport_test.csv")

# Check if they align now
all(genotype$Accession.name == coordinates$Accession.name)

head(genotype[,1:3])
genotype<-genotype[,-1]

head(coordinates[,1:3])
accession_names <- coordinates$Accession.name
coordinates <-coordinates[,-1]


library(maps)

plot(coordinates, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
map(add = T, interior = F)


which(!complete.cases(coordinates))

# Names of the accessions with no georeference (missing Latitude/Longitude)
missing_geo_names <- accession_names[!complete.cases(coordinates)]
missing_geo_names

# Example: remove specific rows from 'genotype'
rows_to_remove <- c(1, 2, 3, 4, 71, 72, 73, 74, 75, 76, 77,
                    115, 166, 167, 168, 183, 316, 417, 423)

genotype <- genotype[-rows_to_remove, ]
coordinates <- coordinates[-rows_to_remove, ]




tess3.obj <- tess3(X = genotype, coord = as.matrix(coordinates), K = 1:8, 
                   method = "projected.ls", ploidy = 2, openMP.core.num = 4) 


plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")


# retrieve tess3 Q matrix for K = 4 clusters 
q.matrix <- qmatrix(tess3.obj, K = 4)
# STRUCTURE-like barplot for the Q-matrix 
barplot(q.matrix, border = NA, space = 0, 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp


## Use CreatePalette() to define color palettes.
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 

my.colors <- c("tomato", "orange", "lightblue", "olivedrab")
my.palette <- CreatePalette(my.colors, 9)
barplot(q.matrix, border = NA, space = 0, 
        main = "Ancestry matrix", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 

plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10),  
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)

library(ggplot2)
library(rworldmap)

map.polygon <- getMap(resolution = "low")

pl <- ggtess3Q(q.matrix, coordinates, map.polygon = map.polygon)
pl +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  xlim(-20, 100) +     # Covers Africa (-20) to far East Asia (150)
  ylim(-25, 40) +      # From South Africa (-40) to Central Asia (60)
  coord_equal() + 
  geom_point(data = as.data.frame(coordinates), aes(x = Longitude, y = Latitude), size = 0.2) + 
  xlab("Longitute") +
  ylab("Latitude") + 
  theme_bw()

#####Get outlier loci
# retrieve tess3 results for K = 4 
p.values <- pvalue(tess3.obj, K = 4)
hist(p.values, col = "lightblue") 

# Benjamini-Hochberg algorithm
# BH correction using R's built-in method
adj.p <- p.adjust(p.values, method = "BH")

# Select candidates below FDR level
candidates <- which(adj.p <= 0.2)
length(candidates)

# manhattan plot 
plot(p.values, main = "Manhattan plot", 
     xlab = "Locus id", 
     ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates, -log10(p.values)[candidates], 
       pch = 19, cex = .2, col = "blue")


#######Diveristy by chromsome

div<-read.table("finger_millet_diversity.txt", header = T)
head(div)

div$new_column <- paste(div$Chromosome, div$StartSite, sep = "_")
div$Chromosome<-as.factor(div$Chromosome)


library(ggplot2)
library(dplyr)

# Example: assuming div has Chromosome, new_column, PiPerBP
# Make a chromosome shading data frame

chrom_bands <- div %>%
  group_by(Chromosome) %>%
  summarise(
    xmin = min(new_column, na.rm = TRUE),
    xmax = max(new_column, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(ymin = -Inf, ymax = Inf)

print(chrom_bands)

# Plot
ggplot(div, aes(x = new_column, y = PiPerBP, group = Chromosome, color = Chromosome)) +
  geom_rect(
    data = chrom_bands,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Chromosome),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),   # remove x-axis title
    axis.text.x  = element_blank(),   # remove x-axis tick labels
    axis.ticks.x = element_blank()    # remove x-axis ticks
  )


###########pull bioclim data
library(terra)

coordinates = read.csv("millet_passport_test.csv")
coords<-coordinates[,-1]

# Path to your directory of GeoTIFFs
tif_dir <- "/Users/mikey/Downloads/wc2.1_30s_bio"

# List all .tif files
tif_files <- list.files(tif_dir, pattern = "\\.tif$", full.names = TRUE)

# Load all rasters as a single multi-layer SpatRaster
r <- rast(tif_files)


# Convert to SpatVector
points <- vect(coords, geom = c("Longitude", "Latitude"), crs = crs(r))

# Extract values
values <- extract(r, points)

# Join back with coords
result <- cbind(coordinates, values)
head(result)

#write.csv(result,"finger_millet_bioclim.csv",row.names=F)

