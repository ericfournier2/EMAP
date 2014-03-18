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