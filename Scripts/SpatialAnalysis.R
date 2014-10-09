library(ggplot2)

# Determine which genomic regions are "hot spots" of methylation activity
# by averaging the p-values of differential expressions of the probes within
# a certain window.
# Arguments:
#   Fit:         The fit object from the limma analysis.
#   annotations: Annotations for all probes.
# Returns:
#   A data frame giving the following information for each probe:
#      - The probe's position.
#      - The average p-value of differential methylation of all probes within a 200k window.
#      - The number of probes within a 200k window.
doHotSpotDetection <- function(Fit, annotations) {
	# Sort annotations according to Chromosome, position.
	annotations <- annotations[order(annotations$Chromosome, annotations$Fragment_Start),]
	
	# Match annotations and Fit probes, then get them together in a single matrix.
	matchIDs <- match(annotations$Probe, Fit$genes$ID)
	annotFit <- cbind(annotations$Fragment_Start, annotations$Fragment_End, Fit$p.value[matchIDs])
	annotChr <- annotations$Chromosome[!is.na(annotFit[,1])]
	annotFit <- annotFit[!is.na(annotFit[,1]),]
	
	# Loop over probes, doing the mean of the p-value of nearby probes.
	windowSize <- 100000
	results <- matrix(nrow=nrow(annotFit), ncol=5)
	for(probe in 1:nrow(annotFit)) {
		closeProbes <- annotFit[max(0,(probe-100)):min(nrow(annotFit), (probe+100)),]
		closeEnough <- closeProbes[(closeProbes[,2] > (annotFit[probe, 1] - windowSize)) & (closeProbes[,1] < (annotFit[probe,2] + windowSize)),]
		if(length(closeEnough) > 3) {
			results[probe,] <- c(annotFit[probe,], nrow(closeEnough), mean(closeEnough[,3], na.rm=TRUE))
		} else {
			results[probe,] <- c(annotFit[probe,], 0, 1)
		}
	}
	
	return(data.frame(Chr=annotChr, Start=results[,1], End=results[,2], PValue=results[,3], CloseProbes=results[,4], MeanPValue=results[,5]))
}


# Plot the average fold-change over variable sized windows for individual chromosomes.
# Parameters:
#   fitData:
#       The limma fit that should be plotted.
#   annotationData:
#       The annotation for the microarray platform.
#   windowSize:
#       The size of the window over which fold-changes will be averaged.
#   stepSize:
#       The step size for the iterative averaging procedure. For example,
#       with a windowSize of 100 and a step size of 10, the averaged 
#       intervals would be [0-100], [10-110], [20-120], etc.
#   pValThreshold:
#       Only probes with a p-value below this threshold will be included
#       in the averaging procedure.
#   fcThreshold:
#       Only probes with an absolute threshold above this value will be
#       included in the averaging procedure.
#   filename:
#       Name to be given to the output files.       
#   chromosomes:
#       Which chromosomes should be plotted. "All" is a special value:
#       if it is present, all chromosomes are plotted.
generateSpatialFCProfile <- function(fitData, annotationData,
                                     windowSize=10000000, stepSize=500000,
                                     pValThreshold=1, fcThreshold=0,
                                     filename="",
                                     chromosomes=c("All")) {
    # Reorder annotations by their chromosome/start position.
    annotations <- annotationData[order(annotationData$Chromosome, annotationData$Fragment_Start),]
    
    # Match annotations and Fit probes, then get them together in a single matrix.
	matchIDs <- match(annotations$Probe, fitData$genes$ID)
	annotFit <- data.frame(Chr=annotations$Chromosome,
                           Start=annotations$Fragment_Start,
                           End=annotations$Fragment_End,
                           FC=fitData$coefficients[matchIDs],
                           PVal=fitData$p.value[matchIDs])
                           
    # Remove probes with undefined coordinates.
    annotFit <- annotFit[!is.na(annotFit$Start),]
    
    # We'll accumulate results in a data-frame.
    modelDF <- data.frame(Chr=character(0),     # The chromosome of the window 
                          Start=numeric(0),     # The start position of the window (in bp)
                          End=numeric(0),       # The end position of the window (in bp)
                          MeanPos=numeric(0),   # The middle position of the window (in bp)
                          MeanPerc=numeric(0),  # The relative position (in percentage) of the middle of the window 
                          FC=numeric(0))        # The average fold-change over the window

    # Process special chromosomes value.
    if("All" %in% chromosomes) {
        chromosomes <- unique(annotFit$Chr)
    }
    
    # Loop over all chromosomes:
    for(chr in unique(chromosomes)) {
        # Keep only the relevant probes.
        chrFit <- annotFit[annotFit$Chr==chr & abs(annotFit$FC) > fcThreshold & annotFit$PVal < pValThreshold,]
        
        # Reset window position, and define maximum window bound.
        lowerBound = 0
        maxBound = max(chrFit$End, na.rm=TRUE)
        
        # Iterate over the full chromosome length. Stop when the upper bound reaches the chromosome's end.
        while(lowerBound + windowSize < maxBound) {
            upperBound = lowerBound + windowSize
        
            # Determine which probes fall within our window, and average the fold-changes..
            probesInRange <- chrFit$Start >= lowerBound & chrFit$End <= upperBound
            meanFC = mean(chrFit$FC[probesInRange], na.rm=TRUE)
            middlePos = (lowerBound+upperBound)/2
            
            # Append the new entry to the data frame.
            newRow = data.frame(Chr=chr,
                                Start=lowerBound,
                                End=upperBound,
                                Mean=middlePos,
                                MeanPerc=middlePos/maxBound,
                                FC=meanFC)
            modelDF <- rbind(modelDF, newRow)
        
            # Move our window downstream.
            lowerBound = lowerBound + stepSize
        }
    }

    # Reorder/rename chromosomes.
    modelDF$Chr <- gsub("chr", "", modelDF$Chr)
    modelDF$Chr <- factor(modelDF$Chr, levels=c(1:29, "X"))

    # Generate a loess curve fitting the data average.
    loessObj <- loess(modelDF$FC ~ modelDF$MeanPerc, degree=2)
    loessPred <- predict(loessObj, seq(0, 1, by=0.01))  
    loessPredDF <- data.frame(Loess=loessPred, X=seq(0, 1, by=0.01))
    
    # Plot it.
    ggplot() +
        geom_hline(y=0) +
        geom_line(data=modelDF, mapping=aes(x=MeanPerc, y=FC, color=Chr)) +
        geom_line(data=loessPredDF, mapping=aes(y=Loess, x=X), color="black", linetype="dashed", width=2) +
        ylab(paste("Average log2(fold-change) over", windowSize/1000000, "Mbp window")) +
        xlab("Relative position on chromosome") +
        scale_colour_discrete(name="Chromosome")

    ggsave(paste(filename, " - All.tiff", sep=""), width=14, height=7, dpi=600, compression="lzw")
        
    for(chr in unique(modelDF$Chr)) {
        ggplot() +
            geom_hline(y=0) +
            geom_line(data=modelDF[modelDF$Chr==chr,], mapping=aes(x=MeanPerc, y=FC, color=Chr)) +
            geom_line(data=loessPredDF, mapping=aes(y=Loess, x=X), color="black", linetype="dashed", width=2) +
            ylab(paste("Average log2(fold-change) over", windowSize/1000000, "Mbp window")) +
            xlab("Relative position on chromosome") +
            scale_colour_discrete(name="Chromosome")    
        ggsave(paste(filename, " - ", chr, ".tiff", sep=""), width=14, height=7, dpi=600, compression="lzw")            
    }

    # Write down the results in a file.
    write.table(modelDF, file=paste(filename, ".txt", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}

# Generate a manhattan-plot like representation of p-values.
generateSpatialPValProfile <- function(fitData, annotationData,
                                     pValThreshold=1, fcThreshold=0,
                                     filename="",
                                     chromosomes=c("All")) {
    # Reorder annotations by their chromosome/start position.
    annotations <- annotationData[order(annotationData$Chromosome, annotationData$Fragment_Start),]
    
    # Match annotations and Fit probes, then get them together in a single matrix.
	matchIDs <- match(annotations$Probe, fitData$genes$ID)
	annotFit <- data.frame(Chr=annotations$Chromosome,
                           Start=annotations$Fragment_Start,
                           End=annotations$Fragment_End,
                           FC=fitData$coefficients[matchIDs],
                           PVal=fitData$p.value[matchIDs])
                           
    # Remove probes with undefined coordinates.
    annotFit <- annotFit[!is.na(annotFit$Start),]

    # Process special chromosomes value.
    if("All" %in% chromosomes) {
        chromosomes <- unique(annotFit$Chr)
    }
    
    # Loop over all chromosomes:
    for(chr in unique(chromosomes)) {
        # Keep only the relevant probes.
        chrFit <- annotFit[annotFit$Chr==chr & abs(annotFit$FC) > fcThreshold & annotFit$PVal < pValThreshold,]
        
        ggplot(data=chrFit) +
            geom_point(mapping=aes(x=Start, y=-log10(PVal))) +
            geom_hline(y=-log10(0.05), color="red")
            
        ggsave(paste("PVal profile for ", chr, ".tiff", sep=""), dpi=600, compression="lzw", width=14, height=7)
    }
}

# Calculates the ratio of DM probes to total probes for each chromosome.
getChromosomeDMRatio <- function(diffExpr, annotationData) {

    # Calculate overall DM and hyper-methylation rates.
    overallDMRate <- nrow(diffExpr) / nrow(annotation)
    overallHyperRate <- sum(diffExpr$Coef > 0)/nrow(diffExpr)
    
    # Loop over all chromosomes and accumulate chromosome-specific DM/hypermethylation rates.
    results <- data.frame(Chromosome=character(0), TotalProbes=numeric(0),
                          DMProbes=numeric(0), Ratio=numeric(0), DMEnrichment=numeric(0),
                          HyperRate=numeric(0), HyperEnrichment=numeric(0))
    for(chr in unique(annotationData$Chromosome)) {
        if(chr != "") {
            # Determine which probes lie on this chromosome, and how many there are.
            thisChr = annotationData$Chromosome==chr
            total = sum(thisChr, na.rm=TRUE)
            
            # Determine how many of these are DMRs.
            dmrOnChr <- diffExpr$ID %in% annotationData$Probe[thisChr]
            dm = sum(dmrOnChr, na.rm=TRUE)
            ratio = dm/total
            dmEnrichment = log2(ratio/overallDMRate)
            
            # Determine how many of these are hypermethylated.
            hyperProbes = sum(diffExpr$Coef[dmrOnChr] > 0, na.rm=TRUE)
            hyperRate = hyperProbes / dm
            hyperEnrichment = log2(hyperRate/overallHyperRate)
            
            results <- rbind(results, data.frame(Chromosome=chr, TotalProbes=total,
                                                 DMProbes=dm, DMRatio=ratio, DMEnrichment=dmEnrichment,
                                                 HyperProbes=hyperProbes, HyperRatio=hyperRate, HyperEnrichment=hyperEnrichment))
        }
    }
    
    # Plot DM ratios.
    ggplot(data=results) +
        geom_bar(mapping=aes(x=Chromosome, y=DMRatio), stat="identity") +
        geom_hline(y=overallDMRate, linetype="dashed", color="red") +
        theme(axis.text.x = element_text(angle = 45, hjust=1))
    ggsave("Chromosome DM Ratio.tiff", compression="lzw", dpi=600)
    
    # Plot Hyper-methylation ratios.
    ggplot(data=results) +
        geom_bar(mapping=aes(x=Chromosome, y=HyperRatio), stat="identity") +
        geom_hline(y=overallHyperRate, linetype="dashed", color="red") +
        theme(axis.text.x = element_text(angle = 45, hjust=1))
    ggsave("Chromosome Hyper Ratio.tiff", compression="lzw", dpi=600)
    
    write.table(results, file="Chromosome DM Ratio.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    return(results)
}
    
