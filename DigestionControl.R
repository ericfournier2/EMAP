library(ggplot2)

# Plots quality control graphs for the digestion control probes of the array.
generateDigestionPlots <- function(intensityData, annotation, filenames) {
	# Get the digestion controls' intensity data and related annotations,
	if(VERSION=="v1") {
        digData <- intensityData[substr(intensityData$genes$ID, 0, 11)=="GT_MET_CTRL",]
        digAnnot <- annotation[substr(annotation$Probe, 0, 11)=="GT_MET_CTRL",]
    } else {
        digData <- intensityData[substr(intensityData$genes$ID, 0, 8)=="EDMA_DIG",]
        digAnnot <- annotation[substr(annotation$Probe, 0, 8)=="EDMA_DIG",]    
    }
	
	# Reorder the digestion controls based on their genomic position.
	digAnnot <- digAnnot[order(digAnnot$Chromosome, digAnnot$Fragment_Start),]
	orderedDigData <- digData[match(digAnnot$Probe, digData$genes$ID),]
	
	# Determine the background cutoff.
	cutOffs <- calculateCutoffs(intensityData[intensityData$genes$ID=="(-)3xSLv1",])

    for(array in 1:ncol(orderedDigData$R)) {	
        # For each array, generate a plot showing the signal for the digestion controls.
        allDig <- data.frame(Channel=c(rep("Red", nrow(orderedDigData)), rep("Green", nrow(orderedDigData))), 
                             Chromosome=c(as.character(digAnnot$Chromosome), as.character(digAnnot$Chromosome)),
                             Value=log2(c(orderedDigData$R[,array], orderedDigData$G[,array])))

        # Remove digestion controls without annotations.
        allDig <- allDig[allDig$Chromosome != "",]
        allDig$Chromosome <- sub( "chr", "", allDig$Chromosome)

        # Reorder chromosomes.
        allDig$Chromosome <- factor(allDig$Chromosome, levels=c(1:29, "X"))
        cutoffsDF <- data.frame(Channel=c("Red", "Green"), Value=log2(c(cutOffs[array], cutOffs[ncol(orderedDigData$R)+array])))
    
		ggplot(allDig, aes(x=Chromosome, y=Value, fill=Channel)) +
            geom_boxplot() +
            geom_hline(data=cutoffsDF, mapping=aes(yintercept=Value), linetype="dashed") +
            facet_grid(Channel~.) +
            scale_fill_manual(values=c("#55CC00", "#FF1100")) +
            scale_y_continuous(name="log2(Intensity)") +                 # Change the label of the y axis.
            labs(title="Intensity of genomic digestion controls") +      # Add a plot title.
            theme(panel.grid.major.x = element_blank(),                  # Remove vertical grid lines.
                strip.text.x = element_text(size=12, face="bold"),       # Embiggen x label
                strip.text.y = element_text(size=12, face="bold"),       # Embiggen y label
                axis.text = element_text(colour="black"),                # Set axis labels to black instead of grey.
                axis.title = element_text(face="bold"),                  # Set axis title to bold
                panel.background = element_rect(fill='#F0F0F0'),         # Make the plot background a bit lighter
                legend.position = "none")                                # Do not display the legend
        ggsave(paste("Genomic digestion for", filenames[array], ".png", sep=""))
	}
}
