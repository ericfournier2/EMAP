# Main analysis file for epigenetic analysis.

# Define error codes
ERROR_UNKNOWN <- 1
ERROR_INVALID_EPIGENETIC_FORMAT <- 2
ERROR_INVALID_TRANSCRIPTOMIC_FORMAT <- 3
ERROR_INVALID_NUMBER_OF_ARGUMENTS <- 4

# Enable printing of the stack on error.
options(error= function() { traceback(2); quit(save="no", status=ERROR_UNKNOWN)})

# Parse command line arguments if values were not provided through an OPEN_ME script.
if(!exists("epigenetic_Name")) {
    parameters <- commandArgs(trailingOnly = TRUE)
    if(length(parameters) != 13) {
        quit(save="no", status = ERROR_INVALID_NUMBER_OF_ARGUMENTS);
    }
    
    epigenetic_Name <- parameters[1]
    epigenetic_Target <-  parameters[2]
    epigenetic_Folder <-  parameters[3]

    transcriptomic_Name <-  parameters[4]
    transcriptomic_Target <-  parameters[5]
    transcriptomic_Folder <-  parameters[6]

    reference_Condition <-  parameters[7]
    
    combined_Name <- parameters[8]
    
    epi_foldchange_Threshold <- log2(as.numeric(parameters[9]))
    epi_pvalue_Threshold <- as.numeric(parameters[10])
    
    trans_foldchange_Threshold <- log2(as.numeric(parameters[11]))
    trans_pvalue_Threshold <- as.numeric(parameters[12])
    
    VERSION <- parameters[13]
}

# There are only a handful of v1 project so by default, if those projects didn't
# explicitly declare themselves as "v1", we suppose we are working with v2.
if(!exists("VERSION")) {
    VERSION <- "v2"
}

library(limma)

# Load all utility files.
source("LoadData.R")
source("LimmaAnalysis.R")
source("LimmaAnalysisControl.R")
source("DigestionControl.R")
source("SpikeControl.R")
source("CategoryEnrichment.R")
source("SpatialAnalysis.R")
source("GenerateBedgraph.R")
source("BuildGeneAnnotation.R")

# Utility function to calculate a cutoff for an array.
calculateCutoffs <- function(intensityData) {
    negControlIndices <- intensityData$genes$ID=="(-)3xSLv1"
    probeSet <- intensityData[negControlIndices,]
    return(c( mean(probeSet$R, na.rm=TRUE) + 4*apply(probeSet$R, 2, sd, na.rm=TRUE),
              mean(probeSet$G, na.rm=TRUE) + 4*apply(probeSet$G, 2, sd, na.rm=TRUE)))
}

# Utility function to get the non-reference condition from a target object.
getOtherCondition <- function(targetObject, referenceCondition) {
    return(setdiff(unique(c(targetObject$Cy3, targetObject$Cy5)), referenceCondition))
}

writeEnrichmentData <- function(enrich, folder, relativeOnly) {
    # Create and move to the output folder.
    outputFolder <- file.path("Enrichment Analysis", folder)
    dir.create(outputFolder, showWarnings=FALSE, recursive=TRUE)
    previousWD <- getwd()
    setwd(outputFolder)
    
    # Plot the enrichment data.
    plotEnrichmentData(enrich, reference_Condition, getOtherCondition(epigeneticsData$Target, reference_Condition), relativeOnly)
    
    # Replace line-breaks in captions.
    for(i in 1:(length(enrich))) {
        enrich[[i]]$"Relative DMR Count" <- gsub("\n", " ", enrich[[i]]$"Relative DMR Count")
    }

    # Write out the raw enrichment data.
    if(relativeOnly=="") {
        write.table(enrich$GeneRegion, file="Genic region/Enrichment - Gene Regions.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
        write.table(enrich$Proximity, file="Distance from CpG Island/Enrichment - CpG Proximity.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
        write.table(enrich$Length, file="CpG Island Length/Enrichment - CpG Length.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
        write.table(enrich$Density, file="CpG Island Density/Enrichment - CpG Density.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
        write.table(enrich$RepeatClasses, file="Repeat/Enrichment - Repeat Classes.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
    } else {
        write.table(enrich$GeneRegion, file="Enrichment - Gene Regions.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
        write.table(enrich$Proximity, file="Enrichment - CpG Proximity.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
        write.table(enrich$Length, file="Enrichment - CpG Length.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
        write.table(enrich$Density, file="Enrichment - CpG Density.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
        write.table(enrich$RepeatClasses, file="Enrichment - Repeat Classes.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
    }

    # Output TFBS information if it is present.
    if(!is.null(enrich$TFBS)) {
        # Add FDR corrected p-values to the TFBS data-frame.
        enrich$TFBS <- cbind(enrich$TFBS,
                             "p-value-low-fdr"=p.adjust(enrich$TFBS$"p-value-low", method="fdr"),
                             "p-value-high-fdr"=p.adjust(enrich$TFBS$"p-value-high", method="fdr"))
        
        write.table(enrich$TFBS, file="Enrichment - TFBS.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
    }
    
    setwd(previousWD)
}

annotationFolder <- file.path("Annotations", VERSION)
divergentScalePath <- file.path(getwd(), "Annotations", "DivergenceScaleNoLabel.png")

# Determine path of pre-computed gene indices/counts
geneMapPath <- file.path(file.path(getwd(), annotationFolder, "GeneMap.txt"))
probeGeneIDPath <- file.path(getwd(), annotationFolder, "ProbeGeneIDs.txt")

###################################################################################################
#                                     Epigenetic Analysis                                         #
###################################################################################################
if(epigenetic_Name!="") {
    # Load epigenetic data and annotations.
    epigeneticsData <- loadData(epigenetic_Folder, epigenetic_Target)
    if(!performSanityCheck(epigeneticsData, epigenetic_Folder)) {
        quit(save="no", status=ERROR_INVALID_EPIGENETIC_FORMAT)
    }
    otherCondition <- getOtherCondition(epigeneticsData$Target, reference_Condition)

    annotation <- read.table(file.path(annotationFolder, "EDMA.Annotation.txt"), header=TRUE, sep="\t")
    categories <- read.table(file.path(annotationFolder, "Categories.txt"), header=TRUE, sep="\t")

    # We're done reading data: move to output directory
    dir.create(file.path("Results", epigenetic_Name), showWarnings=FALSE, recursive=TRUE)
    oldWD <- getwd()
    setwd(file.path("Results", epigenetic_Name))

    # Do the differential expression analysis                      
    limmaResults <- doLimmaAnalysis(epigeneticsData$Target, epigeneticsData$IntensityData, epi_foldchange_Threshold, epi_pvalue_Threshold, reference_Condition)
    dmProbesPerGene <- getNumberOfDMProbesPerGene(limmaResults$DiffExpr)
    
    # Generate QC plots
    dir.create("QC plots", showWarnings=FALSE, recursive=TRUE)
    setwd("QC plots")
    generateMAPlots(MA.RG(epigeneticsData$IntensityData, bc.method="none"), epigeneticsData$Target$Filename, " - Raw")
    generateMAPlots(limmaResults$Norm, epigeneticsData$Target$Filename, " - Normalized")
    generateDigestionPlots(epigeneticsData$IntensityData, annotation, epigeneticsData$Target$Filename)
    # Generate Spike plot based on version.
    if(VERSION=="v1" || VERSION=="pigv1") {
        # Get additional information about the spikes.
        spikeTable <- read.table(file.path(oldWD, annotationFolder, "Spikes.txt"), sep="\t", header=TRUE)    
        generateSpikePlots(epigeneticsData$IntensityData, epigeneticsData$Target$Filename, spikeTable)
    } else {
        generateSpikePlots(epigeneticsData$IntensityData, epigeneticsData$Target$Filename)
    }
    setwd("..")
    
    # Output results of the limma analysis
    dir.create("Limma Analysis", showWarnings=FALSE, recursive=TRUE)
    setwd("Limma Analysis")
    outputLimmaResults(limmaResults, annotation, reference_Condition, otherCondition, categories=categories)
    write.table(dmProbesPerGene, file="DMProbePerGene.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    generateVolcanoPlot(limmaResults$Fit, epi_foldchange_Threshold, epi_pvalue_Threshold, epigeneticsData$Target, reference_Condition, "Hyper-methylated")
    generateAboveBackgroundPlots(epigeneticsData, limmaResults, reference_Condition)
    setwd("..")

    # Generate bedgraph files.
    dir.create("Bedgraph", showWarnings=FALSE, recursive=TRUE)
    setwd("Bedgraph")    
    generateBedGraphs(limmaResults, epigeneticsData$Target, annotation, epi_foldchange_Threshold, epi_pvalue_Threshold)
    setwd("..")
    
    # Perform enrichment analysis.
    enrichDMRs <- enrichmentAnalysis(limmaResults$DiffExpr, annotation)
    enrichRef <- enrichmentAnalysis(subsetFit(limmaResults$Fit, limmaResults$AboveBG$Probes$AllReference), annotation)
    enrichOther <- enrichmentAnalysis(subsetFit(limmaResults$Fit, limmaResults$AboveBG$Probes$AllOther), annotation)

    writeEnrichmentData(enrichDMRs, "DMRs", "")
    writeEnrichmentData(enrichRef, reference_Condition, reference_Condition)
    writeEnrichmentData(enrichOther, otherCondition, otherCondition)
    
    # Perform hot-spot detection.
    hotSpots <- doHotSpotDetection(limmaResults$Fit, annotation)
    write.table(hotSpots, file="HotSpots.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
    generateSpatialFCProfile(limmaResults$Fit, annotation,
                             filename="Spatial chromosome FC Analysis",
                             chromosomes=paste("chr", c(1:9, "X"), sep=""))
    
    # Return to the old working directory.
    setwd(oldWD)
}


###################################################################################################
#                                  Transcriptomic analysis                                        #
###################################################################################################
if(transcriptomic_Name!="") {
    # Load transcriptomic data and probe positions. 
    transcriptomicsData <- loadData(transcriptomic_Folder, transcriptomic_Target)
    if(!performSanityCheck(transcriptomicsData, transcriptomic_Folder)) {
        quit(save="no", status=ERROR_INVALID_TRANSCRIPTOMIC_FORMAT)
    }
    bedTrans <- read.table("Annotations/Circos-BESTv1.bed", sep="\t", col.names=c("BEDChromosome", "Start", "End", "Probe"))
    otherCondition <- getOtherCondition(transcriptomicsData$Target, reference_Condition)
    annotationTrans <- read.table("Annotations/EMBV3.annotation_table_2.xls", header=TRUE, sep="\t", quote="")
    
    # Move into output directory
    dir.create(file.path("Results", transcriptomic_Name), showWarnings=FALSE, recursive=TRUE)
    oldWD <- getwd()
    setwd(file.path("Results", transcriptomic_Name))

    # Do limma analysis
    limmaResultsTrans <- doLimmaAnalysis(transcriptomicsData$Target, transcriptomicsData$IntensityData, trans_foldchange_Threshold, trans_pvalue_Threshold, reference_Condition)

    # Output limma results.
    outputLimmaResults(limmaResultsTrans, annotationTrans, reference_Condition, otherCondition)
    generateVolcanoPlot(limmaResultsTrans$Fit, trans_foldchange_Threshold, trans_pvalue_Threshold, transcriptomicsData$Target, reference_Condition, "Over-expressed")
    
    
    # Get positionned and ordered transcriptomic data.
    posTransData <- cbind(ID=limmaResultsTrans$Fit$genes$ID,
                          PVal=as.vector(-log10(limmaResultsTrans$Fit$p.value)),
                          FC=as.vector(limmaResultsTrans$Fit$coefficients),
                          bedTrans[match(limmaResultsTrans$Fit$genes$ID, bedTrans[,4]), 1:3])
    posTransData <- posTransData[order(posTransData$BEDChromosome, posTransData$Start),]
                          
    # Generate p-value, Fold-change bedgraph files.
    generateBedGraph(posTransData[!is.na(posTransData$BEDChromosome),c("BEDChromosome", "Start", "End", "PVal")], "Trans-P-Value")
    generateBedGraph(posTransData[!is.na(posTransData$BEDChromosome),c("BEDChromosome", "Start", "End", "FC")], "Trans-FC") 

    # Generate differentially expressed circos track.
    #posTransData <- posTransData[(posTransData$PVal < -log10(trans_pvalue_Threshold)) & (abs(posTransData$FC) > trans_foldchange_Threshold),]
    posTransData <- posTransData[posTransData$ID %in% limmaResultsTrans$DiffExpr$ID,]
    diffExprTrans <- cbind(posTransData[,c("BEDChromosome", "Start", "End", "PVal"),],
                           ifelse(posTransData[,3]<0, "color=red,angle_shift=180", "color=green"))
    generateBedGraph(diffExprTrans, "DiffExpr-P-Value-Trans")

    # Restore old working directory.
    setwd(oldWD)
}

###################################################################################################
#                                  Cross-platform analysis                                        #
###################################################################################################
if(epigenetic_Name!="" && transcriptomic_Name!="") {
    # Match epigenetic probes to transcriptomic probes
    epiToTransMatch <- match(annotation$EMBV3_Probe, limmaResultsTrans$Fit$genes$ID)

    # Match epigenetic probes to their annotation, and obtain numerical data which is in the same
    # order as the annotations.
    annToInt <- match(annotation$Probe,  epigeneticsData$IntensityData$genes$ID)
#    allMeans <- apply(cbind(epigeneticsData$IntensityData$R, epigeneticsData$IntensityData$G), 1, mean)
#    matchedInt <- allMeans[annToInt]
    matchedFC <- limmaResults$Fit$coefficients[annToInt]
    matchedP <- limmaResults$Fit$p.value[annToInt]

    # Find genes whose methylation/expression go up/down or down/up at the same time.
    indices <- (abs(matchedFC)> epi_foldchange_Threshold) &
      (abs(limmaResultsTrans$Fit$coefficients[epiToTransMatch]) > trans_foldchange_Threshold) &
      (matchedP < epi_pvalue_Threshold) &
      (limmaResultsTrans$Fit$p.value[epiToTransMatch] < trans_pvalue_Threshold) &
      ((matchedFC*limmaResultsTrans$Fit$coefficients[epiToTransMatch])<0)

    # Remove NAs and keep only relevant rows
    concordant <- annotation[which(indices),, drop=FALSE]

    # Get single, most relevant gene name:
    # Paste together all gene symbols in order of importance. Detect exon/intron by their 
    # exon/intron number suffix then delete everything which follows it.
    geneName <- gsub("-.*", "", paste(concordant$Exon, concordant$Intron, concordant$Proximal_Promoter, concordant$Promoter))
    # If the gene symbol came from the promoter, it won't have an exon number, and it might be doubled. Remove leading spaces
    # as a first step to only keeping the frontmost entry.
    geneName <- gsub("^ *", "", geneName)
    # Once the leading space is removed, remove everything AFTER a space, which will include any secondary/tertiary gene symbols.
    geneName <- gsub(" .*", "", geneName)    

    # Build a new data frame
    concordantLabels <- cbind(as.character(concordant$Chromosome),
                              concordant$Fragment_Start,
                              concordant$Fragment_End,
                              geneName)

    # Remove entries without gene names in the Exon/Intron/Promoter part.
    concordantLabels <- concordantLabels[concordantLabels[,4]!="",, drop=FALSE]

    # Remove duplicate gene names.
    concordantLabels <- concordantLabels[!duplicated(concordantLabels[,4]),, drop=FALSE]

    # Write out the table.                          
    dir.create(file.path("Results", combined_Name), showWarnings=FALSE, recursive=TRUE)
    write.table(concordantLabels, file=file.path("Results", combined_Name, "Concordant.txt"),
                quote=FALSE, row.names=FALSE, col.names=FALSE)             
}
