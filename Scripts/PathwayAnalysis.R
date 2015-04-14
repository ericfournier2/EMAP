# Performs pathway analysis, determining which pathways are
# enriched in a dataset and producing graphical representations
# of these pathways with integrated score-data, when provided.

# Load required libraries.
library(pathview)
library(Heatplus)

# Load up saved list of pathways.
load(file.path(speciesFolder, "gsets.RData"))
speciesAbbreviation <- c(cow="bta", pig="ssc")

# Determines if members of certain KEGG pathways are significantly overrepresented
# in a given set of probes, and generates a graphical representation of the pathway
# were genes are colored according to user-specified scores.
# Parameters:
#   chosenProbes:
#       The set of probes which are part of the "chosen" subset.
#       For example, these could be the probes which are part of a given module,
#       or the differentially expressed genes.
#   universeSubset:
#       The set of all probes that should be considered in the enrichment analysis.
#       For examples, if modules were generated from all probes with an ANOVA
#       p-value below 0.05, the "chosen" set could be the probes in the module,
#       and the "universe" set would be all probes with an ANOVA p-value below 0.05.
#   scores:
#       Scores associated with all genes in the universe set. Used for generating
#       graphical representations of the enriched pathways.
# Returns:
#       The enrichment results for the chosen probes.
performPathwayEnrichment <- function(chosenProbes, universeSubset, annotation, scores=NULL) {
    # Make sure the inputs are correct.
    if(!all(chosenProbes %in% universeSubset)) {
        stop("Error: Chosen probes are not all part of the given subset")
    }
    
    if(!all(as.character(universeSubset) %in% as.character(annotation$Probe))) {
        stop("Error: Given probe subset is invalid.")
    }

    if(!is.null(scores) && length(scores) != length(universeSubset)) {
        stop("Error: Score must be NULL or of the same length as universeSubset.")
    }
    
    # Prepare the list of valid IDs.
    if(species=="cow") {
        allGeneIDs <- gsub("GeneID:", "", annotation$GeneID)
        allGeneIDs[!grepl("^\\d+$", allGeneIDs)] <- ""

        allValidIDs <- allGeneIDs[allGeneIDs != ""]
    } else {
        allGeneIDs <- annotation$GeneID
        allValidIDs <- allGeneIDs[allGeneIDs != ""]
    }

    # Which Entrez IDs are part of the chosen subset?
    chosenGeneIDs <- allGeneIDs[annotation$Probe %in% chosenProbes]
    chosenGeneIDs <- chosenGeneIDs[chosenGeneIDs != ""]
    
    # Which Entrez IDs are part of the universe subset?
    universeGeneIDs <- allGeneIDs[annotation$Probe %in% universeSubset]
    universeGeneIDs <- universeGeneIDs[universeGeneIDs != ""]   
    
    # Data-frame for output.
    results <- data.frame(Pathway=character(0),  # The pathway of interest.
                          PVal=numeric(0),       # The p-value of enrichment within the chosen subset.
                          Chosen=numeric(0),     # The number of Entrez IDs within the "Chosen" subset which are part of this pathway.
                          Universe=numeric(0),   # The number of Entrez IDs within the "Universe" subset.
                          Possible=numeric(0),   # The number of Entrez IDs which could potentially have been part of the "Chosen" subset.
                          Drawn=numeric(0),      # The number of Entrez IDs within the "Chosen" subset.
                          Expected=numeric(0))   # The expected number of Entrez IDs which should be in the "Chosen" subset if the distribution was random.

    # Loop over all non-disease pathway.                          
    for(pathway in names(gsets$kg.sets[-gsets$dise.idx])) {
        # Determine the set of pathway IDs which could possibly be drawn.
        allPathwayIDs <- gsets$kg.sets[[pathway]]
        allPathwayIDs <- allPathwayIDs[allPathwayIDs %in% universeGeneIDs]
    
        # Calculate metrics for the hypergeometric tests.
        chosen <- sum(unique(chosenGeneIDs) %in% unique(allPathwayIDs))
        universe <- length(unique(universeGeneIDs))
        possible <- length(unique(allPathwayIDs))
        drawn <- length(unique(chosenGeneIDs))
        expected <- possible*(drawn/universe)
    
        # Perform the hypergeometric test.
        pVal <- phyper(chosen, possible, universe - possible, drawn, lower.tail=FALSE)
        
        # Add this pathway to the results data.frame.
        results <- rbind(results, data.frame(Pathway=pathway, PVal=pVal, Chosen=chosen, Universe=universe, 
                                             Possible=possible, Drawn=drawn, Expected=expected))
    }

    # write out the results.
    write.table(results, file="KEGG Pathway enrichment.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    # If scores were provided, generate graphical representation of the significantly enriched pathways.
    if(!is.null(scores)) {
        for(pathway in which(results$PVal < 0.05 & results$Chosen > 0)) {
            # pathview expects the data to have row/element names corresponding to Entrez Gene IDs.
            names(scores) <- universeSubset
            
            # Only keep scores from constitutive probes to avoid duplicates.
            if(species=="cow") {
                scoresSubset <- scores[annotation$Target_Location[match(as.character(names(scores)), as.character(annotation$Probe))]=="Constitutive"]
                names(scoresSubset) <- allGeneIDs[match(as.character(names(scoresSubset)), as.character(annotation$Probe))]            
            } else {
                scoresSubset <- scores
                names(scoresSubset) <- allGeneIDs[match(as.character(names(scoresSubset)), as.character(annotation$Probe))]
            }
            
            # Generate the graphical representation for this pathway.
            pathwayID <- gsub(paste("(^", speciesAbbreviation[species], "\\d+).*", sep=""), "\\1", results$Pathway[pathway])
            pathview(scoresSubset, pathway.id=pathwayID, , species=speciesAbbreviation[species], kegg.dir=annotationKEGG)
        }
    }
    
    return(results)
}

# Generate a heatmap which compares and clusters the GO enrichment results
# from a set of orthogonal modules.
# Parameters:
#   enrichmentList:
#       A named list containing the results of multiple calls to performTopGOEnrichment.
#       All calls must have the same "universe" set of probes.
#   inclusionThreshold:
#       The p-value threshold for inclusion into the analysis. A GO term
#       must have a p-value under the threshold in at least one module
#       to be included in the clustering analysis.
compareTopKEGGSets <- function(enrichmentList, inclusionThreshold=0.001) {
    # Loop over all modules in the input list:
    firstLoop <- TRUE
    for(module in names(enrichmentList)) {
        # Fetch the enrichment results from the input list structure.
        resTable <- enrichmentList[[module]]
        resTable$PVal[resTable$PVal==0] <- 10e-10
        testResults <- -log10(as.numeric(resTable$PVal))

        if(firstLoop) {
            # On the first iteration of the loop, create a new data frame in the enrichmentDFList.
            enrichmentDF <- data.frame(KEGG=resTable$Pathway)
            tableOrdering <- seq(1:length(testResults))
        } else {
            # On subsequent iterations, reorder to results to match the initial GO ordering.
            tableOrdering <- match(as.character(enrichmentDF[,1]), as.character(resTable$Pathway))
        }
        
        # Add the results to the data.frame
        enrichmentDF[,module] <- as.numeric(testResults[tableOrdering])

        firstLoop <- FALSE
    }

    # Remove the GO ID column, and keep only the numerical data.
    rawData <- enrichmentDF[,-1]
    rownames(rawData) <- gsub("^bta\\d+\\s(.*)", "\\1", enrichmentDF[,1])
    
    # Filter out GO term which are not significant in any module.
    dataSubset <- rawData[apply(rawData>-log10(inclusionThreshold), 1, any),]
    
    # Generate the annotated heat map and plot it.
    ahm <- annHeatmap2(as.matrix(dataSubset), legend=TRUE, scale="row", 
                       labels=list(nrow=10, 
                                   Row=list(labels=substr(rownames(dataSubset), 1, 30))))
    tiff(filename="KEGG Comparison.tiff", width=9, height=9, units="in", res=600, compression="lzw")
    plot(ahm)
    dev.off()

}

