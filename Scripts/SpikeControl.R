library(ggplot2)

generateSpikePlots <- function(intensityData, filenames, spikeTable=NULL) {
	# Get the spike intensity data.
	if((VERSION=="v1") || (VERSION=="pigv1")) {
        spikeData <- intensityData[substr(intensityData$genes$ID, 0, 10)=="GT_MET_SPK",]
        spikeData <- spikeData[!is.na(spikeData$genes$ID),]
        # Remove all unused spikes
        spikeData <- spikeData[spikeData$genes$ID %in% spikeTable$Probe[spikeTable$Enzyme != "None"],]
        typeCode <- c("Positive", "Negative") 
    } else {	# v2
        spikeData <- intensityData[substr(intensityData$genes$ID, 0, 8)=="EDMA_SPK",]
        spikeData <- spikeData[!is.na(spikeData$genes$ID),]
        typeCode <- c("Pos", "Neg") 
    }
	
    # The way the probes are named and the full enzyme/type combination are different.
    enzymeNames <- c("AciI", "HpaII", "HinP1I")
    enzymeCode <- c("Aci", "Hpa", "Hin")
    typeName <- c("Methylated", "Non-methylated")
    
	# For all arrays
	for(array in 1:ncol(intensityData$R)) {
        # Generate data.frame to hold data.
        allSpikes2 <- data.frame(Channel=rep("Red", nrow(spikeData) * 2), 
                                 Enzyme=rep("HpaII", nrow(spikeData) * 2),
                                 Type=rep("Pos", nrow(spikeData) * 2),
                                 Value=rep(1.5, nrow(spikeData) * 2),
                                 stringsAsFactors=FALSE)
        counter <- 1

        cutoffs <- calculateCutoffs(intensityData[,array])
        names(cutoffs) <- c("Red", "Green")

        # Loop over colours.
        for(color in c("R", "G")) {
            colorName <- "Red"
            if(color=="G") {
                colorName <- "Green"
            }
            
            # Loop over enzymes, control type.            
            for(enzyme in c(1:3)) {
                for(controlType in c(1:2)) {
                    # Get the indices of the relevant controls.
                    if(VERSION=="v1" || VERSION=="pigv1") {
                        spikeNames <- spikeTable$Probe[(spikeTable$Enzyme == enzymeNames[enzyme]) & (spikeTable$Type == typeCode[controlType])]
                        indices <- spikeData$genes$ID %in% spikeNames
                    } else {
                        indices <- (substr(spikeData$genes$ID, 10, 12) == enzymeCode[enzyme]) & (substr(spikeData$genes$ID, 14, 16) == typeCode[controlType])
                    }
                    
                    # Determine which rows of the data frame we should be filling.
                    numElements <- sum(indices, na.rm=TRUE)
                    uBound <- counter + numElements - 1
                    
                    # Set data frame members.
                    allSpikes2$Channel[counter:uBound]  <- colorName
                    allSpikes2$Enzyme[counter:uBound]  <- enzymeNames[enzyme]
                    allSpikes2$Type[counter:uBound]  <- typeName[controlType]

                    # Map the raw [cutoff, 65535] interval to [0,100] as a "percentage of undigested material".
                    allSpikes2$Value[counter:uBound]  <- ((spikeData[[color]][indices,array] - cutoffs[colorName]) / (65535 - cutoffs[colorName])) * 100
                    
                    # Alternate option #1: Standard log2 of intensity.
                    #allSpikes2$Value[counter:uBound]  <- log2(spikeData[[color]][indices,array])
                    # Alternate option #2: Log2 of intensity, recentered at cutoff (IE, cutoff level = 0).
                    #allSpikes2$Value[counter:uBound]  <- log2(spikeData[[color]][indices,array]) - log2(cutoffs[colorName])
                    # Alternat option #3: Map [cutoff, 16] to [0,100] using log2 of intensities. Another version of percentage of undigested material. 
                    #allSpikes2$Value[counter:uBound]  <- ((log2(spikeData[[color]][indices,array]) - log2(cutoffs[colorName])) / (16 - log2(cutoffs[colorName]))) * 100
                    
                    # Update row counter.
                    counter <- counter + numElements
                }
            }
        }
        
        allSpikes2 <- allSpikes2[1:(counter-1),]      
       
        allSpikes2$Type <- factor(allSpikes2$Type, levels=c("Methylated", "Non-methylated"))

        # For standard graph (Alternate option #1): Build a cutoff data.frame to add cutoff lines.
#        cutoffsDF <- data.frame(Channel=c("Red", "Red", "Green", "Green"), Type=c("Methylated", "Non-methylated", "Methylated", "Non-methylated"), Value=log2(c(cutoffs["Red"], cutoffs["Red"], cutoffs["Green"], cutoffs["Green"])))
#        cutoffsDF$Type <- factor(cutoffsDF$Type, levels=c("Methylated", "Non-methylated")

        # Draw plot.
        ggplot(allSpikes2, aes(x=Enzyme, y=Value, fill=Channel)) +            # Plot allSpikes data, enzymes on the x axis, "values" on the y -axis, color=fill according to the channel.
               geom_boxplot() +                                               # Plot as a box plot.
               facet_grid(Channel~Type) +                                     # Split the graph into facets based on color channel (Vertical) and spike type (horizontal).
               theme(panel.grid.major.x = element_blank(),                    # Remove vertical grid lines.
                     strip.text.x = element_text(size=12, face="bold"),       # Embiggen x label
                     strip.text.y = element_text(size=12, face="bold"),       # Embiggen y label
                     axis.text = element_text(colour="black"),                # Set axis labels to black instead of grey.
                     axis.title = element_text(face="bold"),                  # Set axis title to bold
                     panel.background = element_rect(fill='#F0F0F0'),         # Make the plot background a bit lighter
                     legend.position = "none") +                              # Do not display the legent
               scale_fill_manual(values=c("#55CC00", "#FF1100")) +            # Change the fill of the box-plots to green/red.
               scale_y_continuous(name="Percentage of undigested material") + # Change the label of the y axis.
               coord_cartesian(ylim=c(0,103)) +                               # Set the plot's limits on the y-axis.
               labs(title="Digestion of spike controls")                      # Add a plot title.
              # Alternate option 1: add lines for cutoff, center on [0,8.5].
              # geom_hline(data=cutoffsDF, mapping=aes(yintercept=Value), linetype="dashed") +               
              # coord_cartesian(ylim=c(0,8.5))
        ggsave(filename=paste("Spike digestion for ", filenames[array], ".png", sep=""))
	}
}