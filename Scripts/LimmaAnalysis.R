# Load required libraries
library(limma)
library(grid)           # Grid is necessary for adding custom labels in the plot's margins.
library(ggplot2)
library(VennDiagram)
library(scales)

listAboveBG <- function(intensityData, targetData, reference_Condition) {
    # Determine which probes are above the background for each array.
    isAboveBG <- t(apply(cbind(intensityData$R, intensityData$G), 1, ">", calculateCutoffs(intensityData)))
    colnames(isAboveBG) <- c(targetData$Cy5, targetData$Cy3)
    otherCondition <- getOtherCondition(targetData, reference_Condition)

    # Plot a diagram showing the overlap of probes above the background in ALL arrays of a given condition.
    aboveBGRef <- apply(isAboveBG[,colnames(isAboveBG)==reference_Condition], 1, sum)
    aboveBGOther <- apply(isAboveBG[,colnames(isAboveBG)==otherCondition], 1, sum)
    
    # Place results in a data frame.
    aboveBGProbes <- data.frame(ID=intensityData$genes$ID, Reference=aboveBGRef, Other=aboveBGOther,
                                AllReference=aboveBGRef==sum(colnames(isAboveBG)==reference_Condition),
                                AllOther=aboveBGOther==sum(colnames(isAboveBG)==otherCondition))
                                
    aboveBGArray <- data.frame(Tissue=colnames(isAboveBG), Count=apply(isAboveBG, 2, sum, na.rm=TRUE))
    
    return(list(Probes=aboveBGProbes, Arrays=aboveBGArray))
}

subsetFit <- function(fitObject, indices) {
    return(data.frame(ID=fitObject$genes$ID,
                      Coef=as.vector(fitObject$coefficients),
                      PVal=as.vector(fitObject$p.value),
                      RawPVal=as.vector(fitObject$raw.p.value),
                      AdjPVal=as.vector(fitObject$adj.p.value))[indices,])
}

# Perform differential methylation analysis using the Limma package.
# Arguments:
#    targetData:       Information about which samples were hybridized
#    intensityData:    Intensity data. 
#    foldChangeCutoff: Fold-change cutoff used in creating the DiffExpr object.
#    pValueCutoff:     P-Value cutoff used in creating the DiffExpr object.
# Returns:
#    A list with the folowwing elements:
#      Norm:     Normalized intensity data.
#      Fit:      Linear fit data.
#      DiffExpr: List of differentially expressed genes.
doLimmaAnalysis <- function(targetData, intensityData, foldChangeCutoff, pValueCutoff, refCondition, useAdjustedPValue=FALSE) {
    # If no reference condition was provided, pick one randomly (Red channel of first array).
    refCond <- refCondition
    if(is.na(refCondition) || refCondition=="") {
        refCond <- targetData$Cy5[1]
    }
    
    # Determine which probes are above the background.
    aboveBG <- listAboveBG(intensityData, targetData, refCond)
    
    # Perform normalization
	Std_MA_Within <- normalizeWithinArrays(intensityData, method="loess", bc.method="none")
	Std_MA_Between <- normalizeBetweenArrays(Std_MA_Within, method="quantile")

    # Perform linear fit and bayesian correction.
	fitDesign <- modelMatrix(targetData, ref=refCond)	
	fit <- lmFit(Std_MA_Between, design=fitDesign)
	ebayes_fit <- eBayes(fit)
    ebayes_fit$raw.p.value <- ebayes_fit$p.value
    ebayes_fit$adj.p.value <- p.adjust(ebayes_fit$p.value, method="fdr")
	
    if(useAdjustedPValue) {
        ebayes_fit$p.value <- ebayes_fit$adj.p.value
    }
    
    # Determine which probes are differentially expressed/methylated.
	coef <- abs(ebayes_fit$coefficients) > foldChangeCutoff
	pval <- ebayes_fit$p.value < pValueCutoff
	indices <- coef & pval
	indices[is.na(indices)] <- FALSE
    diffExpr <- subsetFit(ebayes_fit, indices)
    
    # Build the return object.
    return(list(Norm=Std_MA_Between, Fit=ebayes_fit, DiffExpr=diffExpr, AboveBG=aboveBG))
}

# Outputs the content of a limmaresults object to text files.
# Arguments:
#    limmaResults: The result returned by the doLimmaAnalysis function.
outputLimmaResults <- function(limmaResults, annotations, reference, otherCondition, categories=NULL) {
    # Build fit results data-frame and output it.
    fitData <- cbind(limmaResults$Fit$genes$ID, limmaResults$Fit$coefficients, limmaResults$Fit$p.value)
    colnames(fitData) <- c("Probe", "Fold-change", "P-value")
    write.table(fitData, file="LimmaAnalysis.txt", quote=FALSE, col.names=, row.names=FALSE, sep="\t")
    
    write.table(limmaResults$DiffExpr, file="DiffExpr.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    
    # Build normalized result data frame and output it.
    rg <- RG.MA(limmaResults$Norm)
    normData <- cbind(limmaResults$Norm$gene$ID, limmaResults$Norm$M, limmaResults$Norm$A, rg$R, rg$G)
    colnames(normData) <- c("Probe", paste(colnames(limmaResults$Norm$M), "M", sep="."),
                                     paste(colnames(limmaResults$Norm$A), "A", sep="."),
                                     paste(colnames(rg$R), "R", sep="."),
                                     paste(colnames(rg$G), "G", sep="."))
    write.table(normData, file="NormData.txt", quote=FALSE, row.names=FALSE, sep="\t")
                
    aboveBG <- limmaResults$AboveBG$Probes
    colnames(aboveBG) <- c("ID",
                           paste("AboveBackground", reference, sep=""),
                           paste("AboveBackground", otherCondition, sep=""),
                           paste("AboveBackgroundAll", reference, sep=""),
                           paste("AboveBackgroundAll", otherCondition, sep=""))
    write.table(aboveBG, file="AboveBG.txt", quote=FALSE, row.names=FALSE, sep="\t")
    
    # Create master output file
    annMatch <- match(fitData[,1], annotations$Probe)
    isDiffExpr <- fitData[,1] %in% limmaResults$DiffExpr$ID
    if(!is.null(categories)) {
        colnames(categories) <- c("LengthCategory", "DensityCategory", "GeneRegionTypeCategory", "ProximityCategory", "DNA", "LINE", "Low_complexity", "LTR", "No_Repeats", "Simple_repeat", "SINE")
        masterResults <- cbind(fitData, DifferentiallyExpressed=isDiffExpr, normData[,-1], aboveBG[,-1], annotations[annMatch,-1], categories[annMatch,], stringsAsFactors=FALSE)
    } else {
        masterResults <- cbind(fitData, DifferentiallyExpressed=isDiffExpr, normData[,-1], aboveBG[,-1], annotations[annMatch,-1], stringsAsFactors=FALSE)
    }
    write.table(masterResults, file="MasterResults.txt", quote=FALSE, row.names=FALSE, sep="\t")
}

#  Generates a volcano-plot
generateVolcanoPlot <- function(fitData, foldchange_Threshold, pvalue_Threshold, target, reference_Condition, overLabel) {
    # Determine which points are significant when compared to threshold
    aboveFC <- (abs(fitData$coefficients) >= foldchange_Threshold)
    belowPValue <- (fitData$p.value <= pvalue_Threshold)
    sigVec <- ifelse( aboveFC & belowPValue, 0, 1)
    sigVec[!aboveFC & !belowPValue] <- 2

    # Count the number of significanthypo and hyper-methylated regions.
    numHypo <- sum(belowPValue & (fitData$coefficients <= -foldchange_Threshold), na.rm=TRUE)
    numHyper <- sum(belowPValue & (fitData$coefficients >= foldchange_Threshold) , na.rm=TRUE)
    
    # Get the name of the condition which is not the reference.
    nonRefCondition <- setdiff(unique(c(target$Cy3, target$Cy5)), reference_Condition)
    
    # Build the data frame.
    volcanoDF <- data.frame(FoldChange=as.vector(fitData$coefficients),
                            PValue=as.vector(-log10(fitData$p.value)),
                            Significant=as.vector(sigVec))
                            
    # Infer limits to create a symmetrical plot along the x-axis
    xWidth <- max(abs(min(volcanoDF$FoldChange, na.rm=TRUE)), max(volcanoDF$FoldChange, na.rm=TRUE))
    xWidth <- xWidth * 1.1

    # Generate the volcano-plot proper.
    p <- ggplot(volcanoDF, aes(x=FoldChange, y=PValue, color=Significant)) +
        geom_point() +
        scale_x_continuous(limits=c(-xWidth, xWidth)) +
        geom_hline(yintercept=-log10(pvalue_Threshold), linetype="dashed") +
        geom_vline(xintercept=-foldchange_Threshold, linetype="dashed") +
        geom_vline(xintercept=foldchange_Threshold, linetype="dashed") +
        labs(x="log2(Fold-Change)",
             y="-log10(p-value)",
             title="Volcano plot of fold-changes and p-values") +                              # Set axis labels
        theme( axis.text = element_text(colour="black"),                                # Set axis text to black
               legend.position="none",                                                  # Remove legend
               plot.margin = unit(c(2.5,1,1,1), "lines"),
               plot.title= element_text(vjust=5.5))

    # Figure out where the plot overLabels should be positioned.
    yTop <- max(volcanoDF$PValue, na.rm=TRUE)
    
    leftLabelX <- (-xWidth - foldchange_Threshold) / 2
    leftLabelY <- yTop + (yTop*0.1)
    rightLabelX <- (xWidth + foldchange_Threshold) / 2
    rightLabelY <- yTop + (yTop*0.1)

    leftCaption <- paste(overLabel, "\nin ", reference_Condition, ": ", numHypo, " regions", sep="")
    rightCaption <- paste(overLabel, "\nin ", nonRefCondition, ": ", numHyper, " regions", sep="")
    
    # Generate the custom labels.
    p <- p +
        annotation_custom(
            grob = textGrob(label = leftCaption, hjust = 0.5, gp = gpar(cex = 1.5, fontsize=8)),
            ymin = leftLabelY,      # Vertical position of the textGrob
            ymax = leftLabelY,
            xmin = leftLabelX,         # Note: The grobs are positioned outside the plot area
            xmax = leftLabelX) +
        annotation_custom(
              grob = textGrob(label = rightCaption, hjust = 0.5, gp = gpar(cex = 1.5, fontsize=8)),
              ymin = rightLabelY,      # Vertical position of the textGrob
              ymax = rightLabelY,
              xmin = rightLabelX,         # Note: The grobs are positioned outside the plot area
              xmax = rightLabelX)

    # Disable clipping to the plot's region.
    gt <- ggplot_gtable(ggplot_build(p))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    
    # Render the plot to a png file.
    png("Volcano plot.png", width = 7, height = 7, units = "in", res=300)
    grid.draw(gt)       
    dev.off()
}

# Generates an MA plot for epigenetic data, with controls marked in a specific color and
# loess curves for all data and for spike controls.
# Arguments:
#    intensityData: An MA object to be plotted.
#    filenames:     The file names associated with the columns of the MA object, for use in the 
#                   plot's title and file name.
#    qualifier:     A qualifier for the MA plot, such as "Raw" or "Normalized".
generateMAPlots <- function(intensityData, filenames, qualifier) {
    # Split all probes into types.
    if(VERSION=="v1") {
        ctrlTypes <- substr(intensityData$genes$ID, 1, 11)
        ctrlTypes[substr(ctrlTypes, 1, 2) != "GT"] <- "Agilent\nControl"
        ctrlTypes[ctrlTypes == "GT_MET_CTRL"] <- "Digestion\nControl"
        ctrlTypes[ctrlTypes == "GT_MET_SPK_"] <- "Spike-in"
        ctrlTypes[substr(ctrlTypes, 1, 5) == "GT_HQ" | substr(ctrlTypes, 1, 5) == "GT_LQ"] <- "Standard\nProbe"
    } else if(VERSION=="v2") {
        ctrlTypes <- substr(intensityData$genes$ID, 1, 8)
        ctrlTypes[substr(ctrlTypes, 1, 4) != "EDMA"] <- "Agilent\nControl"
        ctrlTypes[ctrlTypes == "EDMA_DIG"] <- "Digestion\nControl"
        ctrlTypes[ctrlTypes == "EDMA_MET"] <- "Standard\nProbe"
        ctrlTypes[ctrlTypes == "EDMA_SPK"] <- "Spike-in"    
    } else {	#pigv1
        ctrlTypes <- substr(intensityData$genes$ID, 1, 9)
        ctrlTypes[substr(ctrlTypes, 1, 2) != "GT"] <- "Agilent\nControl"
        ctrlTypes[substr(ctrlTypes, 1, 6) == "GT_DIG"] <- "Digestion\nControl"
        ctrlTypes[ctrlTypes == "GT_pig_Hq" | ctrlTypes == "GT_pig_Lq" | ctrlTypes == "GT_MET_ch" | ctrlTypes == "GT_MET_GL"] <- "Standard\nProbe"
        ctrlTypes[ctrlTypes == "GT_MET_SP"] <- "Spike-in"    	
	}
    
    # Turn types into a factor for color ordering.
    ctrlTypes <- factor(ctrlTypes, levels=c("Standard\nProbe", "Digestion\nControl", "Agilent\nControl", "Spike-in"))

    # For each array in the intensitydata object:
    for(array in 1:ncol(intensityData$M)) {
        # Put data inside a data frame, and remove agilent controls.
        dataDF <- data.frame(M=intensityData$M[,array], A=intensityData$A[,array], Type=ctrlTypes);
        dataDF <- dataDF[order(dataDF$Type),]
        dataDF <- dataDF[dataDF$Type!="Agilent\nControl",]

        # Fit a loess curve to all data.
        completeObs <- is.finite(dataDF$A) & is.finite(dataDF$M)
        save(list=ls(), file="Dump.RData")
        allFit <- lowess(dataDF$A[completeObs], dataDF$M[completeObs], f=0.3)
        allFitDF <- data.frame(x=allFit$x, y=allFit$y)
        
        # Fit a loess curve to digestion and spike controls, giving more weight to the spike controls.
        controlSubset <- dataDF[dataDF$Type=="Digestion\nControl" | dataDF$Type=="Spike-in",]
        controlWeights <- rep(1, nrow(controlSubset))
        controlWeights[controlSubset$Type=="Digestion\nControl"] <- 0.03
        controlFit <- loess(M~A, controlSubset, weights=controlWeights, span=0.75, iterations=5, degree=2)
        
        # Generate the predicted loess curve from the loess fit.
        predictRange <- seq(min(dataDF$A, na.rm=TRUE), max(dataDF$A, na.rm=TRUE), by=0.05)
        controlDF <- data.frame(M=predict(controlFit, predictRange), A=predictRange)
        controlDF <- controlDF[!is.na(controlDF$M),]
        
        # Generate the plot, plotting the points in reverse z-order.
        ggplot() +
            geom_point(data=dataDF[dataDF$Type=="Standard\nProbe",], aes(x=A, y=M, color=Type)) +
            geom_point(data=dataDF[dataDF$Type=="Digestion\nControl",], aes(x=A, y=M, color=Type)) +
            geom_point(data=dataDF[dataDF$Type=="Spike-in",], aes(x=A, y=M, color=Type)) +
            geom_line(data=allFitDF[seq(1, nrow(allFitDF), by=100),], aes(x=x, y=y), linetype="dashed", color="yellow") +
            geom_line(data=controlDF, aes(x=A, y=M), linetype="dashed", color="orange") +
            geom_hline( mapping=aes(yintercept=0), linetype="dashed") +
            scale_colour_manual(values=c("Digestion\nControl"="#0011FF", "Standard\nProbe"="#000000", "Spike-in"="#FF1100", "#FFFFFF")) +
            scale_y_continuous(name="M-value") +                         # Change the label of the y axis.
            scale_x_continuous(name="A-value") +                         # Change the label of the x axis.
            labs(title=paste("MA plot for ", filenames[array], qualifier, sep="")) +         # Add a plot title.
            theme(panel.grid.major.x = element_blank(),                  # Remove vertical grid lines.
                strip.text.x = element_text(size=12, face="bold"),       # Embiggen x label
                strip.text.y = element_text(size=12, face="bold"),       # Embiggen y label
                axis.text = element_text(colour="black"),                # Set axis labels to black instead of grey.
                axis.title = element_text(face="bold"),                  # Set axis title to bold
                panel.background = element_rect(fill='#F0F0F0'))         # Make the plot background a bit lighter
        ggsave(paste("MA plot for ", filenames[array], qualifier, ".png", sep=""))
    }
}

# Defines lower and upper bounds for error bars, but with a minimum of zero.
# To use with bar plots where the error defines an "upper bound".
errorBarDataMinZero <- function(x) {
    results <- c(max(mean(x)-(1.5*sd(x)), 0), mean(x)+(1.5*sd(x)))
    names(results) <- c("ymin", "ymax")
    return(results)
}

# Generates a set of plots to visualize how many and which probes
# are above the background values for all arrays.
generateAboveBackgroundPlots <- function(epigeneticsData, limmaResults, reference_Condition) {
    # Determine what is the non-reference condition.
    otherCondition <- getOtherCondition(epigeneticsData$Target, reference_Condition)

    # Keep only standard probes in the diffexpr set.
    if(VERSION=="v1") {
        diffExpr <- (epigeneticsData$IntensityData$genes$ID %in% grep("^GT_[HL]Q", limmaResults$DiffExpr$ID, value=TRUE))
    } else if(VERSION=="v2") {
        diffExpr <- (epigeneticsData$IntensityData$genes$ID %in% grep("^EDMA_MET", limmaResults$DiffExpr$ID, value=TRUE))
    } else {    # pigv1
        diffExpr <- (epigeneticsData$IntensityData$genes$ID %in% grep("^GT_pig_[HL]q", limmaResults$DiffExpr$ID, value=TRUE))
    }
    
    # Plot a venn diagram of probes which are above the background in all arrays of each condition.
    png("Venn diagram of probes above the background in all arrays of a given condition.png", width = 7, height = 7, units = "in", res=300)
    grid.draw(venn.diagram(list(I   = which(limmaResults$AboveBG$Probes$AllReference),
                                II  = which(limmaResults$AboveBG$Probes$AllOther),
                                III = which(diffExpr)),
                           category=c(reference_Condition, otherCondition, "Differentially expressed"),
                           fill = c("#3399FF", "#FF1100", "yellow"),
                           margin=c(0.05),
                           filename=NULL,
                           main="Probes expressed above the background\nin all arrays of a given condition",
                           main.cex=2,
                           cat.cex=2,
                           cex=2, euler.d=FALSE, scaled=FALSE))
    dev.off()
    
    # Plot an histogram showing in how many arrays probes are above the background level.
    numProbes <- nrow(limmaResults$AboveBG$Probes)
    counts <- data.frame(Count=c(limmaResults$AboveBG$Probes$Reference, limmaResults$AboveBG$Probes$Other),
                         Tissue=c(rep(reference_Condition, numProbes),
                                  rep(otherCondition, numProbes)))
    counts$Count <- as.factor(counts$Count)
    counts <- counts[!is.na(counts$Count),]
    
    ggplot(counts, aes(x=Count, fill=Tissue)) +                         # Set data
        geom_bar(stat="bin", colour="black", position="dodge") +        # Set type (bars)
        labs(x="Number of arrays in which a probe is above the background",
             y="Number of probes") + # Set axis labels
        theme( panel.grid.major.x = element_blank(),                    # Remove x grid lines
               axis.text = element_text(colour="black", size=12)) +     # Set axis text to black
        scale_fill_manual(values=c("#3399FF", "#FF1100"))
    ggsave("Histogram of the number of arrays in which a probe is above the background.png")

    # Plot the number of probes above the background, per array, per condition.
    aboveDF <- limmaResults$AboveBG$Arrays
    pVal <- t.test(aboveDF$Count[aboveDF$Tissue==reference_Condition], aboveDF$Count[aboveDF$Tissue==otherCondition])
    
    pLabel <- paste("p-value:", sprintf("%1.3f", pVal$p.value))
    if(pVal$p.value < 0.001) {
        pLabel <- "p-value < 0.001"
    }    
    
    yMax <- max(errorBarDataMinZero(aboveDF$Count[aboveDF$Tissue==reference_Condition])["ymax"], 
                errorBarDataMinZero(aboveDF$Count[aboveDF$Tissue==otherCondition])["ymax"])
    
    ggplot(aboveDF, aes(x=Tissue, y=Count, fill=Tissue)) + 
        scale_y_continuous(labels=comma, limits=c(0, yMax*1.1)) +
        stat_summary (fun.data="errorBarDataMinZero", geom="errorbar", width=0.25, size=1, mapping = aes (group = 1)) +
        stat_summary (fun.y = mean, geom="bar", mapping = aes (group = 1), colour="black") +
        theme( panel.grid.major.x = element_blank(),                    # Remove x grid lines
               axis.text = element_text(colour="black", size=12)) +     # Set axis text to black
        scale_fill_manual(values=c("#3399FF", "#FF1100")) +
        geom_text(label=pLabel, x=2.5, y=yMax*1.1, hjust=1, vjust=0)
    ggsave("Number of probes above the background.png")
}

# For each gene, counts the number of probes which are differentially
# methylated for that gene.
getNumberOfDMProbesPerGene <- function(diffExpr) {
    # If precomputed indices do not exist, generate them.
    if(!file.exists(file.path(geneMapPath)) || !file.exists(probeGeneIDPath)) {
        print("Generating gene symbol -> numeric ID mapping... This may take a while.\n")
        buildGeneAnnotation()
        print("Done!\n")
    }

    # Load gene symbol -> numeric gene ID map
    geneMap <- read.table(geneMapPath, sep="\t", header=TRUE, as.is=TRUE)
    
    # Load Probe -> numeric gene IDs map
    mapTable <- read.table(probeGeneIDPath, sep="\t", header=TRUE, as.is=TRUE)
    mapList <- lapply(strsplit(mapTable$GeneIDs, " "), as.integer)
    
    # Count the number of differentially expressed probes associated with each
    # gene.
    diffExprInd <- annotation$Probe %in% diffExpr$ID
    results <- rep(0, nrow(geneMap))
    for(idVec in mapList[diffExprInd]) {
        if(length(idVec)!=0) {
            results[idVec] <- results[idVec] + 1
        }
    }

    return(data.frame(Gene=geneMap$Gene, ProbeCount=geneMap$Count, DiffExprProbeCount=results, Ratio=results/geneMap$Count))
}    
