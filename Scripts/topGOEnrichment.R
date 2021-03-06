# Functions used to perform GO enrichment using the topGO library.

library(GO.db)
library(topGO)

probe2goEpi <- readMappings(file.path(annotationFolder, "GO_Filtered.map"))
probe2goTrans <- readMappings(file.path(speciesFolder, "probe2go.map"))
epiSymbols <- read.table(file.path(annotationFolder, "ProbeSymbol.mapping"), sep="\t", header=TRUE)

conditionalPastePrefix = function(x, y) {
    val <- sub("^\\s+", "", x)
    val <- val[val!=""]
    if(length(val) > 0) {
        return(paste(y, "-", val, sep="", collapse=","))
    } else {
        return("")
    }
}

conditionalPaste = function(x) {
    x <- x[x!=""]
    if(length(x) > 0) {
        return(paste(x, sep="", collapse=","))
    } else {
        return("")
    }
}

# Generates a direct mapping between probes and all symbols they are associated with.
# This can be saved, and used later when associating the significant DMR probes within
# a GO term to the genes they target.
generateEpiProbeSymbolMapping <- function() {
    # List all relevant annotation columns.
    columnsOfInterest <- c("Proximal_Promoter", "Promoter", "Exon", "Intron")
    
    # Prepare a matrix which can contain the processed content of each column,
    newMat <- matrix("", nrow=nrow(annotation), ncol=length(columnsOfInterest))
    colnames(newMat) <- columnsOfInterest
    
    # Add an appropriate prefix to all gene names.
    for(i in columnsOfInterest) {
        splitGenes = strsplit(as.character(annotation[,i]), " ")
        newMat[,i] = unlist(lapply(splitGenes, conditionalPastePrefix, i))
    }
    
    # Concatenate all symbols into a single string.
    finalSymbols <- apply(newMat, 1, conditionalPaste)
    finalDF <- data.frame(Probe=annotation$Probe, Gene_Symbol=finalSymbols, stringsAsFactors=FALSE)
    
    # Output the generate data-frame.
    write.table(finalDF, file.path(annotationFolder, "ProbeSymbol.mapping"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}

# Performs the actual calls to topGO to perform the enrichment analysis.
# Parameters:
#   ontology:
#       The ontology to use for the enrichment analysis. One of "BP", "CC" or "MF".
#   allGenes:
#       A named vector containing either scores or a 0/1 two-level factor for each gene/probe.
#       Used to determined ordering/inclusion/exclusion from the selected group.
#   goSubset:
#       The GO annotations for all probes in allGenes.
#   statistic:
#       The test statistic to use. See the documentation for runTest in the topGO package
#       for a list of valid options.
#   outputFilename:
#       A prefix/suffix for the output files generated by this function. This should not
#       contain any directories.
#   ...:
#       All remaining parameters are passed on to the constructor of the topGOdata object.
# Returns:
#   A list containing three elements:
#       Data: The topGOdata object.
#       Test: The result of the runTest function.
#       Table: The summary table provided by the GenTable method.
innerTopGo <- function(ontology, allGenes, goSubset, statistic, outputFilename, symbolAnnotation, ...) {
    # Create the topGO data object, which contains the selected subset, the universe subset, gene scores, etc.
    topGODataObj <- new("topGOdata", ontology=ontology, allGenes = allGenes, annotationFun=annFUN.gene2GO, gene2GO=goSubset, nodeSize=5, ...)
    
    # Perform the enrichment analysis.
    testResultWeight <- runTest(topGODataObj, algorithm="weight01", statistic=statistic)
    testResultClassic <- runTest(topGODataObj, algorithm="classic", statistic=statistic)
    
    # Generate graph of top enriched terms.
    printGraph(topGODataObj, testResultClassic, firstSigNodes = 5, fn.prefix = paste(ontology, outputFilename), useInfo = "all", pdfSW = TRUE)
    
    # Generate a summary table of the enriched terms and write it out.
    summaryTable <- GenTable(topGODataObj, Weight=testResultWeight, Classic=testResultClassic, orderBy="Classic", ranksOf="Classic", topNodes=length(score(testResultClassic)), numChar=2000)
    write.table(summaryTable, file=paste(ontology, " ", outputFilename, ".txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    
    # Generate lists of "significant" probes for each GO term
    allGenesInTerms <- genesInTerm(topGODataObj)
    allSigGenes <- sigGenes(topGODataObj)
    
    sigGenesInTerms <- lapply(allGenesInTerms, function(x) { x[x %in% allSigGenes] })
    
    sigGenesInTermsDF <- data.frame(GO=names(sigGenesInTerms), Probes="", Genes="", stringsAsFactors=FALSE)
    for(i in 1:length(sigGenesInTerms)) {
        goTerm = names(sigGenesInTerms)[i]
        probes = paste(sigGenesInTerms[[i]], collapse=",")
        genes = paste(symbolAnnotation$Gene_Symbol[match(sigGenesInTerms[[i]], symbolAnnotation$Probe)], collapse=",")
        sigGenesInTermsDF[i,] <- data.frame(GO=goTerm, Probes= probes, Genes=genes, stringsAsFactors=FALSE)
    }
    write.table(sigGenesInTermsDF, file=paste(ontology, " - significant genes per term - ", outputFilename, ".txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

    return(list(Data=topGODataObj, Test=testResultClassic, Table=summaryTable))
}

# For use with performTopGOEnrichment as a probeSelectionFun function.
# Returns the set of probes which have a p-value above the
# pValueSelectThreshold global variable in a given module.
pValueSelectThreshold <- 0.05
pValueSelect <- function(x) {
    return(x<pValueSelectThreshold)
}

# For use with performTopGOEnrichment as a probeSelectionFun function.
# Returns the set of probes which have an absolute membership above the
# moduleMembershipSelectThreshold global variable in a given module.
moduleMembershipSelectThreshold <- 0.8
moduleMembershipSelect <- function(x) {
    return(abs(x)>moduleMembershipSelectThreshold)
}

# Perform GO enrichment using the topGO library.
# Parameters:
#   chosenProbes:
#       The name of the probes which are part of the subset whose enrichment should
#       be assessed.
#   subsetOfProbes:
#       The "universe" set, IE all of the probes which could potentially have been
#       part of the chosenProbes set.
#   outputFilename:
#       A name to be used as a suffix/prefix when naming output files. Should not contain
#       directories.
#   probeScores:
#       If a relevant score can be attributed to all probes, this parameter can be used
#       as an alternative to the chosenProbes parameter, and these scores will be used
#       to perform a Kolmogorov-Smirnov statistical test.
#   probeSelectionFun:
#       If the probeScores argument is provided, this should be a function which
#       uses the probe scores to determine if a probe should be included in the
#       "chosen" subset.
#   scoreOrder:
#       Describes probe scores:
#           "increasing" if the lowest score is the most significant (p-values, ranks0
#           "decreasing" if the highest score is the most significant (membership percentage, etc.)
# Returns:
#   A list of three elements, one for each ontology (BP, CC, MF). See innerTopGO for the structure
#   of each individual element.
performTopGOEnrichment <- function(chosenProbes, subsetOfProbes, outputFilename, platform="Epigenetic", 
                                   probeScores=NULL, probeSelectionFun=NULL, scoreOrder="increasing") {
    # Make sure the inputs are correct.
    if(!(platform %in% c("Epigenetic", "Transcriptomic"))) {
        stop("Error: Platform must be either 'Epigenetic' or 'Transcriptomic'")
    }

    if(!all(chosenProbes %in% subsetOfProbes)) {
        stop("Error: Chosen probes are not all part of the given subset")
    }
    if(platform=="Epigenetic") {
        if(!all(as.character(subsetOfProbes) %in% as.character(annotation$Probe))) {
            stop("Error: Given probe subset is invalid.")
        }
    } else {
        if(!all(as.character(subsetOfProbes) %in% as.character(annotationTrans$Probe))) {
            stop("Error: Given probe subset is invalid.")
        }    
    }
    
    if(length(unique(chosenProbes)) != length(chosenProbes)) {
        stop("Error: duplicated probes in chosenProbes")
    }

    if(length(unique(subsetOfProbes)) != length(subsetOfProbes)) {
        stop("Error: duplicated probes in subsetOfProbes")
    }
    
    cat(paste("Selected ", length(chosenProbes), " out of ", length(subsetOfProbes), ".\n", sep=""))

    # Switch annotations depending on the platform.
    if(platform=="Epigenetic") {
        probe2go <- probe2goEpi
        symbolAnnotation <- epiSymbols
    } else {
        probe2go <- probe2goTrans
        symbolAnnotation <- annotationTrans
    }
    
    # Subset probes so that only those with GO annotations are used for the analysis.
    goSubset <- probe2go[names(probe2go) %in% subsetOfProbes]
    
    if(is.null(probeScores)) {
        # When choosing a subset of genes, we must convert them to a 2-level (0,1) factor vector.
        probeScores <- factor(as.integer(names(goSubset) %in% chosenProbes))
        names(probeScores) <- names(goSubset)

        # Without scores, we must use the count-based fisher statistical test for enrichment.
        stat <- "fisher"
        
        cat(paste("After removal of unannotated probes, ", sum(probeScores==1), " out of ", length(goSubset), " remains.\n", sep=""))
    } else {
        # If scores are available, use the Kolmogorov-Smirnov statistical test.
        stat <- "ks"
        
        cat(paste("After removal of unannotated probes, ", length(goSubset), " probes remain.\n", sep=""))
    }

    # Perform enrichment on all ontologies.
    results <- list()
    for(ontology in c("BP", "MF", "CC")) {        
        results[[ontology]] <- innerTopGo(ontology, probeScores, goSubset, stat, paste(ontology, outputFilename),  symbolAnnotation, geneSelectionFun=probeSelectionFun)
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
compareTopGOSets <- function(enrichmentList, inclusionThreshold=0.001) {
    # First, we need to fetch the enrichment values for all GO terms
    # and concatenate them into three separate data.frames, one for each
    # ontology.
    enrichmentDFList <- list()

    # Loop over all modules in the input list:
    firstLoop <- TRUE
    for(module in names(enrichmentList)) {
        # Loop over all ontologies, which are the top-level elements of the module lists.
        for(ontology in c("BP", "CC", "MF")) {
            # Fetch the enrichment results from the input list structure.
            resTable <- enrichmentList[[module]][[ontology]]$Table
            testResults <- -log10(as.numeric(resTable$Weight))

            if(firstLoop) {
                # On the first iteration of the loop, create a new data frame in the enrichmentDFList.
                enrichmentDFList[[ontology]] <- data.frame(GO=resTable$GO.ID)
                tableOrdering <- seq(1:length(testResults))
            } else {
                # On subsequent iterations, reorder to results to match the initial GO ordering.
                tableOrdering <- match(as.character(enrichmentDFList[[ontology]][,1]), as.character(resTable$GO.ID))
            }
            
            # Add the results to the data.frame
            enrichmentDFList[[ontology]][,module] <- as.numeric(testResults[tableOrdering])
        }
        firstLoop <- FALSE
    }

    # Now, for each ontology, generate a heatmap comparing the most relevant GO terms.
    for(ontology in c("BP", "CC", "MF")) {
        # Remove the GO ID column, and keep only the numerical data.
        rawData <- enrichmentDFList[[ontology]][,-1]
        
        # Associate terms to the GO IDs.
        goTermMap <- toTable(GOTERM[as.character(enrichmentDFList[[ontology]][,1])])
        goTermMap <- goTermMap[goTermMap$Ontology==ontology,]
        goTerms <- goTermMap$Term[match(enrichmentDFList[[ontology]][,1], goTermMap$go_id)]
        
        rownames(rawData) <- goTerms
        
        # Filter out GO term which are not significant in any module.
        dataSubset <- rawData[apply(rawData>-log10(inclusionThreshold), 1, any),]
        
        # Generate the annotated heat map and plot it.
        ahm <- annHeatmap2(as.matrix(dataSubset), legend=TRUE, scale="row", 
                           labels=list(nrow=10, 
                                       Row=list(labels=substr(rownames(dataSubset), 1, 30))))
        tiff(filename=paste(ontology, "GO Comparison.tiff"), width=9, height=9, units="in", res=600, compression="lzw")
        plot(ahm)
        dev.off()
    }
}
