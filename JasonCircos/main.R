# Main analysis file for epigenetic analysis.

# create_circos(plot_name, output_dir, expression_file, epigenetic_file, annotation_file)
parameters <- commandArgs(trailingOnly = TRUE)
output_Directory <- parameters[1]
transcriptomic_Data <-  parameters[2]
epigenetic_Data <- parameters[3]
epigenetic_Annotation <- parameters[4]


library(limma)

# Load all utility files.
source("GenerateBedgraph.R")

readELMAFile <- function(filename) {
    elmaData <- read.table(filename, sep="\t", quote="\"", header=TRUE)
    
    # Identify contrasts
    contrasts <- gsub(":logfc" ,"", grep(":logfc", colnames(elmaData), value=TRUE))
    
    resultList <- list()
    for(contrast in contrasts) {
        fc <- elmaData[,paste(contrast, ":logfc", sep="")]
        p <- elmaData[,paste(contrast, ":p_value", sep="")]
        
        resultList[[contrast]] <- data.frame(Probe=elmaData$probe,
                                             FoldChange=fc,
                                             PValue=p)
    }
    
    return(list(Contrasts=contrasts, Results=resultList))
}

###################################################################################################
#                                     Epigenetic Analysis                                         #
###################################################################################################
if(epigenetic_Data != "") {
    # Load epigenetic data and annotations.
    annotation <- read.table(epigenetic_Annotation, header=TRUE, sep="\t")
    elmaData <- readELMAFile(epigenetic_Data)
    
    # We're done reading data: move to output directory
    dir.create(output_Directory, showWarnings=FALSE, recursive=TRUE)
    oldWD <- getwd()
    setwd(file.path("Results", epigenetic_Name))

    # For each contrast
    for(contrast in elmaData[["Contrasts"]]) {
        # Move to contrast specific directory
        dir.create(contrast, showWarnings=FALSE, recursive=TRUE)
        setwd(contrast)
        
        generateBedGraphs(elmaData[["Results"]][[contrast]], 
        
        setwd("..")
    
    }
    
    # Do all analyses: QC plots, enrichment analysis, hot-spot detection, bedgraph files.
    generateBedGraphs(elmaData, epigeneticsData$Target, annotation)

    # Return to the old working directory.
    setwd(oldWD)
}


###################################################################################################
#                                  Transcriptomic analysis                                        #
###################################################################################################
if(transcriptomic_Name!="") {
    # Load transcriptomic data and probe positions. 
    transcriptomicsData <- loadData(transcriptomic_Folder, transcriptomic_Target)
    bedTrans <- read.table("Annotations/BestEMBV3.bed", sep="\t", col.names=c("BEDChromosome", "Start", "End", "Probe"))

    # Move into output directory
    dir.create(file.path("Results", transcriptomic_Name), showWarnings=FALSE, recursive=TRUE)
    oldWD <- getwd()
    setwd(file.path("Results", transcriptomic_Name))

    # Do limma analysis
    limmaResultsTrans <- doLimmaAnalysis(transcriptomicsData$Target, transcriptomicsData$IntensityData, foldchange_Threshold, pvalue_Threshold, reference_Condition)

    # Output limma results.
    write.table(cbind(limmaResultsTrans$Fit$genes$ID, limmaResultsTrans$Fit$coefficients, limmaResultsTrans$Fit$p.value),
            file="LimmaAnalysis.txt", quote=FALSE, col.names=c("Probe", "Fold-change", "P-value"), row.names=FALSE)
    write.table(limmaResultsTrans$DiffExpr, file="DiffExpr.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)

    # Get positionned and ordered transcriptomic data.
    posTransData <- cbind(limmaResultsTrans$Fit$genes$ID,
                          limmaResultsTrans$Fit$p.value,
                          limmaResultsTrans$Fit$coefficients,
                          bedTrans[match(limmaResultsTrans$Fit$genes$ID, bedTrans[,4]), 1:3])
    posTransData <- posTransData[order(posTransData$BEDChromosome, posTransData$Start),]
                          
    # Generate p-value, Fold-change bedgraph files.
    generateBedGraph(posTransData[!is.na(posTransData$BEDChromosome),c("BEDChromosome", "Start", "End", "limmaResultsTrans$Fit$p.value")], "Trans-P-Value")
    generateBedGraph(posTransData[!is.na(posTransData$BEDChromosome),c("BEDChromosome", "Start", "End", "limmaResultsTrans$Fit$coefficients")], "Trans-FC") 

    # Generate differentially expressed circos track.
    posTransData <- posTransData[(posTransData[,2] < pvalue_Threshold) & (abs(posTransData[,3]) > foldchange_Threshold),]
    diffExprTrans <- cbind(posTransData[,c("BEDChromosome", "Start", "End", "limmaResultsTrans$Fit$p.value"),],
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
    allMeans <- apply(cbind(epigeneticsData$IntensityData$R, epigeneticsData$IntensityData$G), 1, mean)
    matchedInt <- allMeans[annToInt]
    matchedFC <- limmaResults$Fit$coefficients[annToInt]
    matchedP <- limmaResults$Fit$p.value[annToInt]

    # Find genes whose methylation/expression go up/down or down/up at the same time.
    indices <- (abs(matchedFC)> foldchange_Threshold) &
      (abs(limmaResultsTrans$Fit$coefficients[epiToTransMatch]) > foldchange_Threshold) &
      (matchedP < pvalue_Threshold) &
      (limmaResultsTrans$Fit$p.value[epiToTransMatch] < pvalue_Threshold) &
      ((matchedFC*limmaResultsTrans$Fit$coefficients[epiToTransMatch])<0)

    # Remove NAs and keep only relevant rows
    concordant <- annotation[which(indices),]

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
    concordantLabels <- concordantLabels[concordantLabels[,4]!="",]

    # Remove duplicate gene names.
    concordantLabels <- concordantLabels[!duplicated(concordantLabels[,4]),]

    # Write out the table.                          
    dir.create(file.path("Results", combined_Name), showWarnings=FALSE, recursive=TRUE)
    write.table(concordantLabels, file=file.path("Results", combined_Name, "Concordant.txt"),
                quote=FALSE, row.names=FALSE, col.names=FALSE)             
}
