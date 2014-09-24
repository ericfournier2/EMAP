# This script creates the "GeneMap.txt" and "ProbeGeneIDs.txt" files which are
# used to assign numeric IDs to all gene symbols, then precompute the set of
# IDs which are associated with each probe's Intron, Exon, Proximal Promoter
# and Distal Promoter fields.

###### To run as stand-alone, uncomment here #####
# Set working directory.
# setwd("C:/Dev/Projects/Epigenetics/cow/EMAP")
#
# Load the appropriate annotation file.
# VERSION <- "v2"
# annotationFolder <- file.path("Annotations", VERSION)
# annotation <- read.table(file.path(annotationFolder, "EDMA.Annotation.txt"), header=TRUE, sep="\t")

buildGeneAnnotation <- function() {
    # Get all gene symbols of interest together in one character string.
    geneNames <- paste(annotation$Exon, annotation$Intron, annotation$Proximal_Promoter, annotation$Promoter)

    # Remove exon/intron number markers,
    geneNames <- gsub("-[0-9]*", "", geneNames)

    # Get rid of extra whitespace.
    geneNames <- gsub(" +", " ", trimWhiteSpace(geneNames))

    # Turn vector of combined IDs into a list of vectors with separate IDs
    geneNamesList <- lapply(strsplit(geneNames, " ", fixed=TRUE), unique)

    # Get a list of all unique gene symbols in the annotation.
    uniqueGenes <- sort(unique(unlist(geneNamesList)))

    # Build a reverse-lookup vector of gene symbols to gene numeric IDs
    geneMap <- 1:length(uniqueGenes)
    names(geneMap) <- uniqueGenes

    # Convert list's gene symbols to list of numeric IDs.
    mapList <- lapply(geneNamesList, function(x) { return(geneMap[x]) })

    # Count the number of probes associated with each gene
    probeCount <- rep(0, length(geneMap))
    for(idVec in mapList) {
        if(length(idVec)!=0) {
            probeCount[idVec] <- probeCount[idVec] + 1
        }
    }

    # Convert back to string of concatenated IDs for saving.
    idStrings <- unlist(lapply(mapList, function(x) { return(ifelse(length(x)==0, "", paste(x, collapse=" "))) }))

    # Save results.
    write.table(data.frame(Gene=names(geneMap), ID=geneMap, Count=probeCount),
                file=geneMapPath,
                quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    write.table(data.frame(Probe=annotation$Probe, GeneIDs=idStrings),
                file=probeGeneIDPath,
                quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}


