generateBedGraph <- function (dataframe, label) {
	# Determine the names of the output file and the track.
	filename <- paste(label, ".bedgraph", sep="")
	
	# Write the track header, required so that bedgraph file will not be parsed as bed files.
	bedGraphFile <- file(filename)
	write(paste("track type=bedGraph name=\"", label, "\" description=\"", label, "\"", sep=""), file=bedGraphFile)
	close(bedGraphFile)
	
	# Write the full data.
	write.table(file=filename, dataframe, append=TRUE, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}


generateBedGraphs <- function(elmaData, annotations, fcThreshold, pvalueThreshold) {
	# Sort annotations according to Chromosome, position.
	annotations <- annotations[order(annotations$Chromosome, annotations$Fragment_Start),]
	annotations  <- annotations[annotations$Chromosome!="",]
	
	# Match annotations and Fit probes, then get them together in a single matrix.
	matchIDs <- match(annotations$Probe, elmaData$Probe)
	
    pValues <- -log10(elmaData$PValue[matchIDs])
    FC <- elmaData$FoldChange[matchIDs]
    
#    normedIntensities <- RG.MA(LimmaResults$Norm)
#    condA <- targets$Cy5[1]
#    condB <- targets$Cy3[1]
#    condMeanA <- log2((apply(cbind(normedIntensities$R[,targets$Cy5==condA], normedIntensities$G[,targets$Cy3==condA]), 1, mean, na.rm=TRUE))[matchIDs])
#    condMeanB <- log2((apply(cbind(normedIntensities$R[,targets$Cy5==condB], normedIntensities$G[,targets$Cy3==condB]), 1, mean, na.rm=TRUE))[matchIDs])
    
	probePValue <- data.frame(annotations$Chromosome, annotations$Probe_Start, annotations$Probe_End, pValues)
	probeFC <- data.frame(annotations$Chromosome, annotations$Probe_Start, annotations$Probe_End, FC)
#    probeCondMeanA <- data.frame(annotations$Chromosome, annotations$Probe_Start, annotations$Probe_End, condMeanA)
#    probeCondMeanB <- data.frame(annotations$Chromosome, annotations$Probe_Start, annotations$Probe_End, condMeanB)    

	generateBedGraph(probePValue, "Probe-P-Value")
	generateBedGraph(probeFC, "Probe-Fold-Change")
	generateBedGraph(probeCondMeanA, "Probe-Cond-Mean-A")
	generateBedGraph(probeCondMeanA, "Probe-Cond-Mean-B")

	fragPValue <- data.frame(annotations$Chromosome, annotations$Fragment_Start, annotations$Fragment_End, pValues)
	fragFC <- data.frame(annotations$Chromosome, annotations$Fragment_Start, annotations$Fragment_End, FC)
#    fragCondMeanA <- data.frame(annotations$Chromosome, annotations$Fragment_Start, annotations$Fragment_End, condMeanA)
#    fragCondMeanB <- data.frame(annotations$Chromosome, annotations$Fragment_Start, annotations$Fragment_End, condMeanB)

	generateBedGraph(fragPValue, "Fragment-P-Value")
	generateBedGraph(fragFC, "Fragment-Fold-Change")
	generateBedGraph(fragCondMeanA, "Fragment-Cond-Mean-A")
	generateBedGraph(fragCondMeanB, "Fragment-Cond-Mean-B")    

    diffExprSubset <- annotations$Probe %in% elmaData$Probe[(abs(elmaData$FoldChange) > fcThreshold) & (elmaData$PValue < pvalueThreshold)]
    diffExprFrag <- data.frame(annotations$Chromosome[diffExprSubset],
                               annotations$Fragment_Start[diffExprSubset],
                               annotations$Fragment_End[diffExprSubset],
                               pValues[diffExprSubset],
                               ifelse(FC[diffExprSubset]<0, "color=red,angle_shift=180", "color=yellow"))
    generateBedGraph(diffExprFrag, "DiffExpr-P-Value")
    
    # Generate imprinted bedgraph.
    # Concatenate all gene symbols for easier parsing
    allText <- paste(annotations$Exon, annotations$Intron, annotations$Proximal_Promoter, annotations$Promoter, annotations$Distal_Promoter, " ")

    # Find all probes which target an imprinted gene.
    annotatedImprints <- grepl("\\bPEG10\\b|MEST\\b|\\bNAP1L5\\b|\\bIGF2R\\b|\\bGNAS\\b|\\bNNAT\\b|\\bDGAT1\\b|\\bMIMT1\\b|\\bUSP29\\b|\\bPEG3\\b|\\bLOC100335527\\b|\\bRTL1\\b|\\bMAGEL2\\b|\\bSNRPN\\b|\\bH19\\b|\\bIGF2\\b|\\bTSSC4\\b|\\bMAOA\\b|\\bXIST\\b", allText)

    # A couple of genes are not in the annotation, so we need to add their intervals manually
    # MIM1: chr18	64291395	64369918
    mim1 <- (annotations$Chromosome=="chr18" & annotations$Fragment_End>=(64291395-50000) & annotations$Fragment_Start<=64291395)

    # USP29: chr18	64535015	64536532
    usp29 <- (annotations$Chromosome=="chr18" & annotations$Fragment_End>=(64535015-50000) & annotations$Fragment_Start<=64536532)

    allImprints <- annotatedImprints | mim1 | usp29
    generateBedGraph(fragFC[allImprints,], "Imprint-Fold-Change")
    
}
