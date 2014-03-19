# This script creates the "Categories.txt" file which is used to add group
# membership information to the epigenetic pipeline output.

# Set working directory.
setwd("C:/Dev/Projects/Epigenetics/EMAP")

# Load the appropriate annotation file.
VERSION <- "v2"
annotationFolder <- file.path("Annotations", VERSION)
annotations <- read.table(file.path(annotationFolder, "EDMA.Annotation.txt"), header=TRUE, sep="\t")

# Load the relevant helper functions.
source("CategoryEnrichment.R")

# Categorize lengths of CpG islands.
thresholds <- quantile(annotations$CpG_Length[annotations$CpG_Length!=0 & !is.na(annotations$CpG_Length)], c(0.20, 0.80))
lengthCategory <- vector(length=nrow(annotations))
lengthCategory[annotations$CpG_Length<thresholds[1]] <- "Small"
lengthCategory[annotations$CpG_Length>thresholds[2]] <- "Long"
lengthCategory[annotations$CpG_Length==0] <- "No-CpG-Island"
lengthCategory[annotations$CpG_Length>=thresholds[1] & annotations$CpG_Length<=thresholds[2]] <- "Intermediate"
lengthCategory[is.na(annotations$CpG_Length)] <- "Unknown"

# Categorize CpG densities.
thresholds <- quantile(annotations$CpG_Density[annotations$CpG_Length!=0 & !is.na(annotations$CpG_Length)], c(0.20, 0.80))
densityCategory <- vector(length=nrow(annotations))
densityCategory[annotations$CpG_Density<thresholds[1]] <- "Low-Density"
densityCategory[annotations$CpG_Density>thresholds[2]] <- "High-Density"
densityCategory[annotations$CpG_Density==0] <- "No-CpG-Island"
densityCategory[annotations$CpG_Density>=thresholds[1] & annotations$CpG_Density<=thresholds[2]] <- "Intermediate-Density"
densityCategory[is.na(annotations$CpG_Density)] <- "Unknown"

# Categorize by type of genic region.
geneRegionTypeCategory <- rep("No-nearby-gene", length=nrow(annotations))
geneRegionTypeCategory[annotations$Distal_Promoter != ""] <- "Distal-Promoter"
geneRegionTypeCategory[annotations$Promoter != ""] <- "Promoter"
geneRegionTypeCategory[annotations$Proximal_Promoter != ""] <- "Proximal-Promoter"
geneRegionTypeCategory[annotations$Intron != ""] <- "Intronic"
geneRegionTypeCategory[annotations$Exon != ""] <- "Exonic"
geneRegionTypeCategory[annotations$Chromosome==""] <- "Unknown"

# Categorize by proximity to CpG Islands
proximityCategory <- as.character(annotations$UCSC_CpG_Proximity)
proximityCategory[proximityCategory=="Shore"] <- "CpG Shore"
proximityCategory[proximityCategory=="Shelf"] <- "CpG Shelf"
proximityCategory[proximityCategory=="Island"] <- "CpG Islands"
proximityCategory[proximityCategory==""] <- "Unknown"
proximityCategory <- sub(" ", "-", proximityCategory, fixed=TRUE)

# Categorize repeat classes.
repeatClasses <- sub("/.*$", "", gsub("/.*? ", " ", annotations$Fragment_RepeatClass))
repeatClasses[repeatClasses==""] <- "No_Repeats"
repeatMatrix <- buildCategoryMatrix(repeatClasses)

# Build a matrix containing all categories.
categories <- cbind(lengthCategory, densityCategory, geneRegionTypeCategory, proximityCategory, 
                    repeatMatrix[,colnames(repeatMatrix) %in% c("No_Repeats", "SINE", "LINE", "Simple_repeat", "LTR", "Low_complexity", "DNA")])
row.names(categories) <- annotations$Probe

# Save the categories.
write.table(categories, file.path(annotationFolder, "Categories.txt"), sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
