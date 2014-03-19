library(ggplot2)
library(reshape2)
library(png)
library(proto)
library(gridExtra)

# Add a "horizontal bar" type to ggplots, based on geom_bar. This is needed for the
# combined graph. Briefly, you cannot combined facets with space="free" and scale="free"
# with the coord_flip() object. However, geom_bar only has vertical bars, so you can't
# have horizontal bars without coord_flip(). This makes it impossible to create the
# combined graph without having all categories show up in all facets.
# By creating a new type which is equivalent to geom_bar, but with x and y reversed, 
# we remove the need for the coord_flip(), and thus we can use free space/scale in the facets.
geom_hbar <- function (mapping = NULL, data = NULL, ...) {
  GeomHBar$new(mapping = mapping, data = data, stat = "identity", position = "identity", ...)
}

GeomHBar <- proto(ggplot2:::Geom, {
  objname <- "bar"
  
  default_stat <- function(.) StatIdentity
  default_pos <- function(.) PositionIdentity
  default_aes <- function(.) aes(colour=NA, fill="grey20", size=0.5, linetype=1, weight = 1, alpha = NA)
  
  required_aes <- c("y")
 
  reparameterise <- function(., df, params) {
    df$width <- df$width %||%
      params$width %||% (resolution(df$y, FALSE) * 0.9)
    transform(df,
      xmin = pmin(x, 0), xmax = pmax(x, 0),
      ymin = y - width / 2, ymax = y + width / 2, width = NULL
    )
  }
 
  draw_groups <- function(., data, scales, coordinates, ...) {
    GeomRect$draw_groups(data, scales, coordinates, ...)
  }
  guide_geom <- function(.) "polygon"
})

# Given a list of character vectors, builds a matrix representing which elements
# are part of which element of the list. For example, given the following list:
#  [[1]] ""
#  [[2]] "A B C"
#  [[3]] "D"
#  [[4]] "B D"
# The function will return the following matrix:
#      A B C D
#  [1] 
#  [2] 1 1 1
#  [3]       1
#  [4]   1   1
#
# Parameters:
#   listOfValues - A list of character vectors
#   fillPresent  - The value to put in the matrix when an element is part of the vector.
#   fillAbsent   - The value to put in the matrix when an element is not part of the vector.
#
# Returns:
#   The matrix indicating the presence of each element in the character vectors.
# 
# Notes
# Obtained from Stack Overflow when looking for a way to optimize my own very slow function:
# http://stackoverflow.com/questions/19594253/optimizing-a-set-in-a-string-list-to-a-set-as-a-matrix-operation/19594838#19594838
# Original author says he will eventually add it to his "splitstackshape" package on CRAN.
charBinaryMat <- function(listOfValues, fillPresent=TRUE, fillAbsent = FALSE) {
  lev <- sort(unique(unlist(listOfValues, use.names = FALSE)))
  m <- matrix(fillAbsent, nrow = length(listOfValues), ncol = length(lev))
  colnames(m) <- lev
  for (i in 1:nrow(m)) {
    m[i, listOfValues[[i]]] <- fillPresent
  }
  m
}

# Wrapper for charBinaryMat which splits the set of individual strings into
# a list of character vectors and removes leading/trailing spaces.
buildCategoryMatrix <- function(str) {
  charBinaryMat(strsplit(sub("^ +", "", str), " ", fixed=TRUE))
}

# Makes all category enrichment calculations for a single category.
# Parameters:
#   allCategoryMembership 
#       A vector of size nrow(annotations) indicating whether each probe within
#       the annotation fits within the given category.
#   annotations
#       The annotation data-frame for the whole chip.
#   diffExpr
#       The DiffExpr object returned by the doLimmaAnalysis function.
# Returns:
#   A vector containing the following information regarding the category enrichment:
#       # of successes          - Number DE/DM probes in the given category
#       # of drawings           - Total number of DE/DM probes
#       # of possible successes - Number of probes in the given category
#       # of cases              - The total number of probes
#       % of successes          - Proportion of DE/DR probes in the given category within all DE/DR probes.
#       % of possible successes - Proportion of DE/DR probes in the given category within all probes.
#       p-value-low             - Significance of the lower tail hypergeometric test
#       p-value-high            - Significance of the "high" tail hypergeometric test
#       # Hyper                 - The number of hypermethylated probes within this category.
#       % Hyper                 - The proportion of hypermethylated DE/DM probes within this category.
#       # HyperAll              - The total number of hypermethylated probes.
#       % HyperAll              - The percentage of hypermethylated DE/DM probes within all probes.
evaluateCategory <- function(allCategoryMembership, annotations, diffExpr) {
    # Determine which probes are part of the differentially expressed list.
    drawnMembership <- annotations$Probe %in% diffExpr$ID
    
    # Calculate enrichment values.
    numSuccess <- sum(allCategoryMembership & drawnMembership)
    numPossibleSuccess <- sum(allCategoryMembership)
    numTotalCases <- length(allCategoryMembership)
    numDrawings <- sum(drawnMembership)
    percSuccess <- numSuccess/numDrawings
    percPossibleSuccess <- numPossibleSuccess/numTotalCases 
    
    diffExprSubset <- diffExpr$ID %in% annotations$Probe[allCategoryMembership]
    hyper <- sum(diffExpr[diffExprSubset,2]>0)
    percHyper <- hyper/sum(diffExprSubset)
    hyperAll <- sum(diffExpr[,2]>0)
    hyperRef <- numSuccess-hyper
    
    pvalLow <- phyper(numSuccess, numPossibleSuccess, numTotalCases - numPossibleSuccess, numDrawings, lower.tail=TRUE)
    pvalHigh <- phyper(numSuccess, numPossibleSuccess, numTotalCases - numPossibleSuccess, numDrawings, lower.tail=FALSE)
    
    results <- c("# of successes"          = numSuccess,
                 "# of drawings"           = numDrawings,
                 "# of possible successes" = numPossibleSuccess,
                 "# of cases"              = numTotalCases, 
                 "% of successes"          = percSuccess, 
                 "% of possible successes" = numPossibleSuccess/numTotalCases, 
                 "p-value-low"             = pvalLow,
                 "p-value-high"            = pvalHigh,
                 "# Hyper Other"           = hyper,
                 "# Hyper Ref"             = hyperRef,
                 "% Hyper Other"           = percHyper,
                 "% Hyper Ref"             = 1-percHyper,
                 "# HyperAll"              = hyperAll,
                 "% HyperAll"              = hyperAll/nrow(diffExpr),
                 "Relative DMR"            = log2(percSuccess/(percPossibleSuccess)),
                 "Relative DMR Count"      = paste(numSuccess, "\n(", sprintf("%2.1f", percSuccess*100), "%)", sep=""),
                 "Relative Hyper"          = log2(percHyper/(1-percHyper)),
                 "Relative Hyper Count"    = paste(numSuccess - hyper, ":", hyper, sep=""),
                 "Enrich Hyper Other"      = log2((hyper/hyperAll)/percPossibleSuccess),
                 "Enrich Hyper Ref"        = log2((hyperRef/(numDrawings-hyperAll))/percPossibleSuccess)
                 )

    return(results)
}

# Given a vector of strings describing which categories all probes fit in, perform category 
# enrichment for all categories present in the vector, minus those in "invCategories", which
# are completely removed prior to the analysis.
#
# Parameters:
#    categories
#       A vector of strings, of the same size as nrow(annotations), describing which categories
#       each probe fits in. If multiple categories are present, they should be split using
#       spaces. Accordingly, category names should not contain spaces.
#    invCategories
#       A list of categories which should be removed from the analysis. Used for "Unknown"
#       categories.
#    annotations
#       The annotation data-frame for the whole chip.
#   diffExpr
#       The DiffExpr object returned by the doLimmaAnalysis function.
#
# Returns
#   A matrix containing the enrichment analysis results.
evaluateCategories <- function(categories, invCategories, annotations, diffExpr) {
    catMatrix <- buildCategoryMatrix(categories)
    
    if(length(invCategories) != 0) {
        invRows <- apply(as.matrix(catMatrix[,invCategories]), 1, any)
        catMatrix <- catMatrix[!invRows,!(colnames(catMatrix) %in% invCategories)]
        diffExpr <- diffExpr[!(diffExpr$ID %in% annotations$Probe[invRows]),]
        annotations <- annotations[!invRows,]
    }
    
    results <- t(apply(catMatrix, 2, evaluateCategory, annotations, diffExpr))
    
    # Convert back to numeric types
    charColumnsIndices <- grepl("Count", colnames(results))
    numericColumns <- results[,!charColumnsIndices]
    textColumns <- results[,charColumnsIndices]
    mode(numericColumns) <- "numeric"
    
    return(data.frame(numericColumns, textColumns, check.names=FALSE))

}

# Simple enrichment analysis for various categories of genes/probes, namely:
#   - Distance from a CpG island
#   - Length of the associated CpG island
#   - Density of the associated CpG island
#   - Type of repeated elements present in the fragment
# Parameters:
#   diffExpr: The list of differentially methylated genes in an experiment
#   annotations: The annotation of all probes
# Returns:
#   A list with the results of the enrichment analysis.
enrichmentAnalysis <- function(diffExpr, annotations) {
	# Categorize lengths of CpG islands.
	thresholds <- quantile(annotations$CpG_Length[annotations$CpG_Length!=0 & !is.na(annotations$CpG_Length)], c(0.20, 0.80))
	lengthCategory <- vector(length=nrow(annotations))
	lengthCategory[annotations$CpG_Length<thresholds[1]] <- "Small"
	lengthCategory[annotations$CpG_Length>thresholds[2]] <- "Long"
	lengthCategory[annotations$CpG_Length==0] <- "No-CpG-Island"
	lengthCategory[annotations$CpG_Length>=thresholds[1] & annotations$CpG_Length<=thresholds[2]] <- "Intermediate"
	lengthCategory[is.na(annotations$CpG_Length)] <- "Unknown"
	
	# Categorize CpG densities.
	thresholds <- quantile(annotations$CpG_Density[annotations$CpG_Length!=0 & !is.na(annotations$CpG_Length)], c(0.20, 0.80))
	densityCategory <- vector(length=nrow(annotations))
	densityCategory[annotations$CpG_Density<thresholds[1]] <- "Low-Density"
	densityCategory[annotations$CpG_Density>thresholds[2]] <- "High-Density"
	densityCategory[annotations$CpG_Density==0] <- "No-CpG-Island"
	densityCategory[annotations$CpG_Density>=thresholds[1] & annotations$CpG_Density<=thresholds[2]] <- "Intermediate-Density"
	densityCategory[is.na(annotations$CpG_Density)] <- "Unknown"
	
    # Categorize by type of genic region.
    geneRegionTypeCategory <- rep("No-nearby-gene", length=nrow(annotations))
    geneRegionTypeCategory[annotations$Distal_Promoter != ""] <- "Distal-Promoter"
    geneRegionTypeCategory[annotations$Promoter != ""] <- "Promoter"
    geneRegionTypeCategory[annotations$Proximal_Promoter != ""] <- "Proximal-Promoter"
    geneRegionTypeCategory[annotations$Intron != ""] <- "Intronic"
    geneRegionTypeCategory[annotations$Exon != ""] <- "Exonic"
    geneRegionTypeCategory[annotations$Chromosome==""] <- "Unknown"
    
    # Categorize by proximity to CpG Islands
    proximityCategory <- as.character(annotations$UCSC_CpG_Proximity)
    proximityCategory[proximityCategory=="Shore"] <- "CpG Shore"
    proximityCategory[proximityCategory=="Shelf"] <- "CpG Shelf"
    proximityCategory[proximityCategory=="Island"] <- "CpG Islands"
    proximityCategory[proximityCategory==""] <- "Unknown"
    proximityCategory <- sub(" ", "-", proximityCategory, fixed=TRUE)
    
	# Categorize repeat classes.
	repeatClasses <- sub("/.*$", "", gsub("/.*? ", " ", annotations$Fragment_RepeatClass))
    repeatClasses[repeatClasses==""] <- "No_Repeats"

    # Transcription factor binding sites (TFBS) classes.
    tfClasses <- as.character(annotations$TFBS)
    tfClasses[tfClasses==""] <- "None"
    
    # Uncategorized probes must be removed to remove the bias they induce.
	result <- list(
        GeneRegion=evaluateCategories(geneRegionTypeCategory, "Unknown", annotations, diffExpr),
		Proximity=evaluateCategories(proximityCategory, "Unknown", annotations, diffExpr),
		Length=evaluateCategories(lengthCategory, "Unknown", annotations, diffExpr),
		Density=evaluateCategories(densityCategory, "Unknown", annotations, diffExpr),
		RepeatClasses=evaluateCategories(repeatClasses, c(), annotations, diffExpr),
        TFBS=evaluateCategories(tfClasses, "None", annotations, diffExpr))
        
    # Fix names and display orders
    relativeOrders=list(Proximity=c("Open Sea", "CpG Shelf", "CpG Shore", "CpG Islands"),
                        Length=c("No CpG Island","Small","Intermediate","Long"),
                        Density=c("No CpG Island","Low Density","Intermediate Density","High Density"),
                        GeneRegion=c("No nearby gene", "Distal Promoter", "Promoter", "Proximal Promoter", "Exonic", "Intronic"),
                        RepeatClasses=c("No Repeats", "SINE", "LINE", "Simple repeat", "LTR", "Low complexity", "DNA"))
    
    # Put everything except TFBS in the correct order for graphical representations.
    for(i in names(relativeOrders)) {
        rownames(result[[i]]) <- gsub("[-_]", " ", rownames(result[[i]]))
        result[[i]] <- result[[i]][relativeOrders[[i]],]
    }
    
    return(result)
}

# Returns the color vector for reference vs other graphs.
getTwoColorVector <- function(refCondition, othercondition) {
    result <- c("#FF1100", "#3399FF")
    names(result) <- c(refCondition, othercondition)
    return(result)
}

# Returns the color vector for multiple categories graph.
# Arguments:
#   colorNames:  The vector of categories to which colors should be attributed.
#   appendWhite: If true, the white color is appended at the end of the list.
#   singleColor: If true, returns a vector containing the same color multiple times.
getColorVector <- function(colorNames, appendWhite=FALSE, singleColor=FALSE) {
    if(singleColor) {
        colorVector <- rep("#3399FF", length(colorNames))
    } else {
        colorVector <- c("#3399FF", "#55CC00", "#FFCC00","#FF1100" , "#FFFF00", "#999999", "#CC0099")
        colorVector <- colorVector[1:length(colorNames)]
        # Should we replace the last color with white?
        if(appendWhite) {
            colorVector[length(colorVector)] <- "#FFFFFF"
        }
        names(colorVector) <- colorNames
    }
    
    return(colorVector)   
}

# Takes a data frame returned by the enrichmentAnalysis function and converts it to
# a data frame appropriate for use with our ggplot plots.
#  1. Row names are converted into an ordered factor and put in the "Category" column.
#  2. Items in the columnMap vector are mapped in the data frame. For example, if one
#     element of the columnMap vector is EnrichPercent="% of possible successes", then
#     the "% of possible successes" column of the original data frame is mapped to the
#     "EnrichPercent" column of the returned data frame.
# Parameters:
#   enrichData      : The raw enrichment data to be converted.
#   columnMap       : The column mapping.
#   insertLineBreak : If true, spaces in the category names are converted to line-breaks.
#   reverseOrder    : If true, category levels are set in the reverse order of the row order.
# Returns:
#   A data frame appropriate for use with ggplot.
mapToDataFrame <- function(enrichData, columnMap, insertLineBreak=FALSE, reverseOrder=FALSE) {
    catValues <- rownames(enrichData)
    if(insertLineBreak) {
        catValues <- gsub(" ", "\n", catValues)
    }
    
    catLevels <- catValues
    if(reverseOrder) {
        catLevels <- rev(catLevels)
    }
    
    result <- data.frame(Category=factor(catValues, levels=catLevels))
    for(i in names(columnMap)) {
        result <- cbind(result, enrichData[,columnMap[i]])
    }
    colnames(result) <- c("Category", names(columnMap))
    
    rownames(result) <- result$Category
    return(result)
}

# Adds a row at the end of a data frame returned by enrichmentAnalysis which creates an
# "All" category which serves as a summary of all other categories.
appendAllRow <- function(enrichData) {
    allRow <- enrichData[1,]
    allRow["# of successes"] <- allRow["# of drawings"]
    allRow["% of successes"] <- 1
    allRow["p-value-low"] <- 1
    allRow["p-value-high"] <- 1
    allRow["# Hyper Other"] <- allRow["# HyperAll"]
    allRow["% Hyper Other"] <- allRow["% HyperAll"]
    allRow["% Hyper Ref"] <- 1 - allRow["% HyperAll"]    
    allRow["Relative DMR"] <- 0
    allRow["Relative Hyper"] <- log2(allRow["% HyperAll"]/(1-allRow["% HyperAll"]))
    allRow["Relative Hyper Count"] <- paste(allRow["# of successes"] - allRow["# Hyper Other"], ":", allRow["# Hyper Other"], sep="")
    
    return(rbind(enrichData, "All"=allRow))
}


# Produces a stacked bar plot comparing the distribution of probes on the whole array with that
# of those within the list of differentially methylated probes.
# Parameters:
#   enData: A matrix, with as many rows as there are categories, and two columns.
#           Column 1 should contain the proportion of probes in this category for the whole array,
#           column 2 should contain the proportion of probes in this category for differentially expressed probes.
#   categoryNames: The ordered display names of the categories.
#   legendName: The name to give to the plot and the categories' legend.
doStackedBarPlot <- function(enrichData, legendName) {
    dataDF <- mapToDataFrame(enrichData, 
                             c("Proportion within\nall EDMA probes"="% of possible successes",
                               "Proportion within\ndifferentially methylated probes"="% of successes"))
                              
    mData <- melt(dataDF, id.vars="Category", variable.name="Type", value.name="Value")
    
    colorVector <- getColorVector(rownames(enrichData))
    
    # Generate and save the plot.
    ggplot(mData, aes(x=Type, y=Value, fill=Category)) +                     # Set plot data.
        geom_bar(stat="identity", colour="black") + labs(x="", y="") +       # Set pot type (stacked bars) and axis labels (none).
        theme( panel.grid.major.x = element_blank(),                         # Remove X grid lines.
               axis.text = element_text(colour="black")) +                   # Set the axis text to black rather than grey.
        scale_fill_manual(name=legendName, breaks=rev(rownames(enrichData)), # Set legend order so it corresponds to the stacked block order.
                          values=colorVector)                                # Set legend colors.
    ggsave(filename=paste(legendName, "enrichment - Stacked.png"),
           width=7, height=7, units="in")                                   # Save plot to file.
}

# Generates a plot showing the percentage of DMRs which are methylated in otherCondition.
# Each bar is accompanied by a percentage showing the ratio the percentage of hypermethylaed
# DMRs in the other condition in this category and the percentage of hypermethylaed DMRs in 
# the other condition for all DMRs.
doHyperPlot <- function(enrichData, categoryNames, legendName, otherCondition) {
    # Get value of the "All" bar:
    hyperAll <- enrichData[,"% HyperAll"][1]
    
    # Replace spaces with line-breaks since category names will be written horizontally.
    rownames(enrichData) <- gsub(" ", "\n", categoryNames, fixed=TRUE)
    
    # Reorder category names so that the first will be up on top.
    categoryNames <- factor(c("All", rownames(enrichData)), levels=c(rownames(enrichData), "All"))
    
    # Build vector of values, which are the values in the input argument prepended with the value
    # of the HyperAll column
    hyperValues <- c(hyperAll, enrichData[,"% Hyper Other"])
    
    # Build labels, which are the proportions of the bar length compared to the "All" bar,
    # formatted as a percentage
    hyperLabels <- sprintf("%2.1f%%", hyperValues/hyperAll*100)
    
    # If we're above the "All" line, but there isn't enough space to put the label inside of the
    # bar without overlapping said line, switch the label to outside of the bar.
    hyperLabelPos <- ifelse((hyperValues > hyperAll) & (hyperValues - 0.15 < hyperAll), hyperAll - 0.01, hyperValues - 0.01)
    # Is the previously chosen position too close to the left edge of the graph?
    tooCloseToLeft <- hyperLabelPos < 0.20
    # If so, move the label to outside of the bar.
    hyperLabelPos <- ifelse(tooCloseToLeft, hyperValues + 0.01 , hyperLabelPos)
    # Finally, if the label is outside of the bar but that it overlaps the "All" line, move it to the right of the "All" line.
    hyperLabelPos <- ifelse(tooCloseToLeft & (hyperValues < hyperAll) & ((hyperValues - 20) < hyperAll), hyperAll + 0.01, hyperLabelPos)
    
    # Label justification: right-justified (1), unless there's nos pace tot he left, in which case it will be left-justified (0).
    hyperLabelJust <- ifelse(tooCloseToLeft, 0, 1)
    
    # Build data-frame for ggplot
    hyperDF <- data.frame(Category=categoryNames, Hyper=hyperValues, Label=hyperLabels,
                          LabelPos=hyperLabelPos, Just=hyperLabelJust)                   
    
    # Match colors.
    colorVector <- getColorVector(hyperDF$Category, appendWhite=TRUE)
    
    # Build the plot
    ggplot(hyperDF, aes(x=Category, y=Hyper, fill=Category)) +                   # Set data
        geom_bar(stat="identity", colour="black") +                              # Set type (bars)
        geom_text(aes(x=Category, y=LabelPos, label=Label, hjust=Just)) +        # Text labels
        geom_hline(yintercept=enrichData[,"% HyperAll"][1], linetype="dotted") + # Dotted line on "All" level.
        ylim(c(0,1)) +                                                           # Always go from 0% to 100%
        labs(x="", y=paste("Proportion of DMRs which are hyper-methylated in", otherCondition)) +          # Set axis labels
        theme( panel.grid.major.x = element_blank(),                             # Remove x grid lines
               axis.text = element_text(colour="black", size=14),                # Set axis text to black
               legend.position="none") +                                         # Remove legend
        scale_fill_manual(values=colorVector) +                                  # Set colors
        coord_flip()                                                             # Turn graphic sideways.

    ggsave(filename=paste(legendName, "enrichment - Hypermethylation.png"),
           width=7, height=7, units="in")
}

# Generates a plot showing the percentage of DMRs which are methylated in each conditions,
# as a stacked bar plot.
doStackedHyperPlot <- function(enrichData, legendName, refCondition, otherCondition) {
    enrichData <- appendAllRow(enrichData)
    hyperDF <- rbind(mapToDataFrame(enrichData, c(Hyper="% Hyper Other"), TRUE),
                     mapToDataFrame(enrichData, c(Hyper="% Hyper Ref"), TRUE))
    hyperDF <- cbind(hyperDF, Tissue=c(rep(otherCondition, nrow(enrichData)), rep(refCondition, nrow(enrichData))))
    
    hyperDF$Category <- factor(hyperDF$Category, levels=hyperDF$Category[1:nrow(enrichData)])
    
    # Build the plot
    ggplot(hyperDF, aes(x=Category, y=Hyper, fill=Tissue)) +                   # Set data
        geom_bar(stat="identity", colour="black") +                              # Set type (bars)
        geom_hline(yintercept=enrichData[,"% HyperAll"][1], linetype="dotted") + # Dotted line on "All" level.
        ylim(c(0,1.0000001)) +                                                   # Always go from 0% to 100%. Add a tiny bit for imprecisions due to rounding.
        labs(x="", y="Proportion of DMRs which are hyper-methylated") +          # Set axis labels
        theme( panel.grid.major.x = element_blank(),                             # Remove x grid lines
               axis.text = element_text(colour="black", size=14)) +              # Set axis text to black
        scale_fill_manual(values=getTwoColorVector(refCondition, otherCondition)) +                          # Set colors
        coord_flip()                                                             # Turn graphic sideways.

    ggsave(filename=paste(legendName, "enrichment - Hypermethylation Stacked.png"),
           width=7, height=7, units="in")
}

doDodgedRelativeHyperPlot <- function(enrichData, legendName, refCondition, otherCondition) {
    hyperDF <- rbind(mapToDataFrame(enrichData, c(Hyper="Enrich Hyper Other"), TRUE),
                     mapToDataFrame(enrichData, c(Hyper="Enrich Hyper Ref"), TRUE))
    hyperDF <- cbind(hyperDF, Tissue=c(rep(otherCondition, nrow(enrichData)), rep(refCondition, nrow(enrichData))))
    
    hyperDF$Category <- factor(hyperDF$Category, levels=hyperDF$Category[1:nrow(enrichData)])
    
    # Build the plot
    ggplot(hyperDF, aes(x=Category, y=Hyper, fill=Tissue)) +                     # Set data
        geom_bar(stat="identity", colour="black", position="dodge") +            # Set type (bars)
        geom_hline(yintercept=0, linetype="solid", size=1) +
        labs(x="", y="log2(Enrichment ratio for hypermethylated DMRs)") +        # Set axis labels
        theme( panel.grid.major.x = element_blank(),                             # Remove x grid lines
               axis.text = element_text(colour="black", size=14)) +              # Set axis text to black
        scale_fill_manual(values=getTwoColorVector(refCondition, otherCondition)) +                      # Set colors
        coord_flip()                                                             # Turn graphic sideways.

    ggsave(filename=paste(legendName, "enrichment - Hypermethylation Relative.png"),
           width=7, height=7, units="in")
}


doRelativePlot <- function(enrichPercent, topLabels, plotName, baseline=0, appendWhite=FALSE, singleColor=FALSE, combined=FALSE, colorColumn="", showCount=TRUE) {
    # Replace -Infinite enrichment scores with -5.
    enrichPercent$EnrichPercent[enrichPercent$EnrichPercent==-Inf] <- -5
     
    graphHeight <- 7
    leftMargin <- 5
    if(nrow(enrichPercent) > 10) {
        graphHeight <- 10
        leftMargin <- 10
    }
     
    # Calculate the offset of labels, based on the total y-span of the plot and the orientation
    # of the bar (toward the bottom or toward the top).
    labelOffsets <- (max(enrichPercent$EnrichPercent)-min(enrichPercent$EnrichPercent)) *       # y-span
                    ifelse(sign(enrichPercent$EnrichPercent)==1, 0.025, 0.05) *                 # Orientation of the bar
                    sign(enrichPercent$EnrichPercent)
    labelJust <- ifelse(sign(labelOffsets)==1, 0, 1)
    enrichPercent <- cbind(enrichPercent, Offset=enrichPercent$EnrichPercent + labelOffsets, Just=labelJust)

    # Create a named color vector so that colors will match those of the stacked-bar plots.
    colorVector <- getColorVector(enrichPercent$Category, appendWhite, singleColor)
    
    # Determine the span of the graph. Take the highest absolute enrichment value,
    # round it up to the closest 0.5 increment and use that or 1.5, whichever is larger.
    ratioEdge <- max(c(abs(min(enrichPercent$EnrichPercent)), max(enrichPercent$EnrichPercent))) + 0.3
    ratioRounded <- ceiling(ratioEdge/0.5)*0.5
    ratioLimit <- max(1.5, ratioRounded)
    
#    if(colorColumn=="") {
        enrichPercent <- cbind(enrichPercent, ColorInfo=enrichPercent$Category)
#    } else {
#        enrichPercent <- cbind(enrichPercent, ColorInfo=enrichPercent[,colorColumn])
#    }
    
    # Generate the main plot.
    gPlot <- ggplot(enrichPercent, aes(y=Category, x=EnrichPercent, fill=ColorInfo)) +  # Set data
             geom_hbar(colour="black") +                                               # Set type (bars)
             geom_vline(xintercept=0, linetype="solid", size=1) +                      # Draw line down the 0 line to delineate both sides.
             geom_vline(xintercept=baseline, linetype="dashed", size=0.25) +           # Draw the dotted "baseline". If 0, will draw over full middle line and be invisible.
             labs(y="", x="log2(Enrichment ratio)") +                                  # Set axis labels
             xlim(c(-ratioLimit, ratioLimit)) +                                        # Set axis limits
             theme( panel.grid.major.y = element_blank(),                              # Remove x grid lines
                    axis.text = element_text(colour="black", size=14),                 # Set axis text to black
                    plot.margin = unit(c(0,1,1,1), "lines"),
                    legend.position="none") +                                          # Remove legend
             scale_fill_manual(values=colorVector)                                     # Set colors               
    if(showCount) {
        gPlot <- gPlot + geom_text(aes(x=Offset,label=Count,hjust=Just))               # Set bar labels
    }

    heightSplit <- c(0.2, 0.8)
    if(combined) {
        gPlot <- gPlot + facet_grid(Categorization~., scale="free", space="free")
        heightSplit <- c(0.12, 0.88)
    }
             
    # Disable clipping in the main grob so that labels can overflow from the plot area.
    enrichGrob <- ggplot_gtable(ggplot_build(gPlot))
    enrichGrob$layout$clip[enrichGrob$layout$name == "panel"] <- "off"              

    # Generate the scale arrows bitmap grob for annotation.
    if(exists(divergentScalePath)) {
        divergentScalePath <- "DivergenceScaleNoLabel.png"
    }
    divergenceScale <- readPNG(divergentScalePath)
    gScale <- rasterGrob(divergenceScale, interpolate=TRUE)

    # Generate the top annotation, including the arrow and the labels.
    labelDF <- data.frame(Label=topLabels, Pos=c(-1,0,1), Y=c(2.5,2.5,2.5))  # Define label positions.
    annot <- ggplot(labelDF, aes(x=Pos, label=Label, y=Y)) +         # Set data.
             geom_text(size=4.3) +                                   # Set the type (text) and its font size.
             labs(y="", x="") + xlim(c(-1.2, 1.2)) + ylim(c(0,4)) +  # Remove axis labels, set axis limits.
             theme( panel.grid = element_blank(),                    # Remove grid lines
                    axis.text = element_blank(),                     # Set axis text to black
                    plot.margin = unit(c(0,1,0,leftMargin), "lines"),# Remove all margins except the left one
                    legend.position="none",                          # Remove the legend
                    panel.background=element_blank(),                # Remove the background
                    axis.ticks=element_blank()) +                    # Remove the ticks.
        annotation_custom(gScale, ymin=-1, ymax=1, xmin=-1.2, xmax=1.2) # Add the arrow scale.

    # Disable clipping so the arrow scale will draw close enough to the actual plot.
    annotGrob <- ggplot_gtable(ggplot_build(annot))
    annotGrob$layout$clip[annotGrob$layout$name == "panel"] <- "off"      

    # Draw both plots in a column.
    allPlots <- arrangeGrob(annotGrob, enrichGrob, nrow=2, heights=heightSplit)

    # Save plot: can't use ggsave for plots drawn through grid. Start a png graphical device.
    png(plotName, width = 7, height = graphHeight, units = "in", res=300)
    grid.draw(allPlots)       
    dev.off()
}

getDMREnrichmentLabels <- function() {
    return(c("Higher\nconservation\nof methylation",
             "Average\ndivergence\nof methylation",
             "Higher\ndivergence\nof methylation"))
}

# Produces a bar plot comparing the distribution of probes on the whole array with that
# of those within the list of differentially methylated probes.
doRelativeBarPlot <- function(enData, plotName, relativeOnly="", singleColor=FALSE) {
    enrichPercent <- mapToDataFrame(enData, c(EnrichPercent="Relative DMR", Count="Relative DMR Count"), TRUE, FALSE)

    topLabels <- getDMREnrichmentLabels()
    if(relativeOnly != "") {
        topLabels <- c(paste("Lower odds\nof methylation\nin", relativeOnly),
                       paste("Average odds\nof methylation\nin", relativeOnly),
                       paste("Higher odds\nof methylation\nin", relativeOnly))
    }
                                
    doRelativePlot(enrichPercent, topLabels, paste(plotName, "enrichment - Relative.png"), 0, FALSE, singleColor)                                
    
    return(enrichPercent)
}

getHyperMethylationLabels <- function(refCondition, otherCondition) {
    return(c(paste("More\nhypermethylation\nin", refCondition),
                   "Hypermethylation\nevenly spread",
                   paste("More\nhypermethylation\nin", otherCondition)))
}

# Produces a bar plot comparing the distribution of probes on the whole array with that
# of those within the list of differentially methylated probes.
doRelativeHyperRatioPlot <- function(enData, plotName, refCondition, otherCondition, singleColor=FALSE) {
    enData <- appendAllRow(enData)
    enrichPercent <- mapToDataFrame(enData, c(EnrichPercent="Relative Hyper", Count="Relative Hyper Count"), TRUE, FALSE)

    topLabels <- getHyperMethylationLabels(refCondition, otherCondition)
    fullName <-  paste(plotName, "enrichment - Relative Hypermethylation.png")
                               
    doRelativePlot(enrichPercent, topLabels, fullName, enrichPercent$EnrichPercent[enrichPercent$Category=="All"], appendWhite=TRUE, singleColor)
                               
    return(enrichPercent)
}


# Produces a bar plot comparing the distribution of probes on the whole array with that
# of those within the list of differentially methylated probes.
doCombinedRelativeBarPlot <- function(enrichDFList, columnNames, plotName, topLabels, addBaseline=FALSE, colorColumn="", showCount=TRUE) {
    # Add a "Categorization" column for facetting.
    for(i in 1:length(enrichDFList)) {
        enrichDFList[[i]] <- cbind(enrichDFList[[i]], Categorization=names(enrichDFList)[i])
    }

    # Concatenate all separate enrichment data.
    finalDF <- rbind(enrichDFList[[1]], enrichDFList[[2]])
    for(i in 3:length(enrichDFList)) {
        finalDF <- rbind(finalDF, enrichDFList[[i]])
    }
    
#    if(colorColumn=="") {
        enrichPercent <- mapToDataFrame(finalDF, c(EnrichPercent=columnNames[1], Count=columnNames[2], Categorization="Categorization"))
#    } else {
#        enrichPercent <- mapToDataFrame(finalDF, c(EnrichPercent=columnNames[1], Count=columnNames[2], Categorization="Categorization", ColorInfo="ColorInfo"))
#    }
    
    baseline <- 0
    if(addBaseline) {
        baselineDF <- appendAllRow(finalDF)
        baseline <- baselineDF[rownames(baselineDF)=="All", columnNames[1]]
    }
    
#    singleColor <- TRUE
#    if(colorColumn!="") {
#        singleColor <- FALSE
#    }
    
    doRelativePlot(enrichPercent, topLabels,  paste(plotName, "enrichment - Combined Relative Hypermethylation.png"),
                   baseline, appendWhite=FALSE, singleColor=TRUE, combined=TRUE, colorColumn="", showCount=showCount)
}

# Generate both a stacked bar plot and a relative bar plot for a set of data.
doPlots <- function(enrichData, legendName, refCondition, otherCondition, relativeOnly) {
    if(relativeOnly == "") {
        doStackedBarPlot(enrichData, legendName)
        doStackedHyperPlot(enrichData, legendName, refCondition, otherCondition)
        doRelativeBarPlot(enrichData, legendName)
        doDodgedRelativeHyperPlot(enrichData, legendName, refCondition, otherCondition)
        doRelativeHyperRatioPlot(enrichData, legendName, refCondition, otherCondition)
    }
}

# Plot all enrichment categories using stacked and side-by-side bars.
plotEnrichmentData <- function(enrich, refCondition, otherCondition, relativeOnly="") {
    # Plot all data types.
    doPlots(enrich$Proximity, "Distance from CpG Island", refCondition, otherCondition, relativeOnly)
    doPlots(enrich$Length, "CpG Island Length", refCondition, otherCondition, relativeOnly)
    doPlots(enrich$Density, "CpG Island Density", refCondition, otherCondition, relativeOnly)
    doPlots(enrich$GeneRegion, "Genic region", refCondition, otherCondition, relativeOnly)
                
    # Do the two graphs that can be directly applied to repeats:
    doRelativeBarPlot(enrich$RepeatClasses,  "Repeat", relativeOnly)    
    if(relativeOnly == "") {
        doStackedHyperPlot(enrich$RepeatClasses, "Repeat",  refCondition, otherCondition)
        doRelativeHyperRatioPlot(enrich$RepeatClasses, "Repeat", refCondition, otherCondition)

        # Now, instead of a stacked bar graph, do a dodged bar graph.
        enData <- enrich$RepeatClasses[order(enrich$RepeatClasses[,"% of possible successes"], decreasing=TRUE),]
        enDataSubset <- enData[,c("% of possible successes", "% of successes")]
        colnames(enDataSubset) <- c("Proportion within\nall EDMA probes", "Proportion within\ndifferentially methylated probes")
        
        mData <- melt(as.matrix(enDataSubset))
        colnames(mData) <- c("Type", "Category", "value") 
        mData$Type <- factor(mData$Type, levels = rownames(enDataSubset))
        ggplot(mData, aes(x=Type, y=value, fill=Category)) +
            geom_bar(stat="identity", colour="black", position="dodge") +
            labs(x="", y="") +
            theme( panel.grid.major.x = element_blank(),
                   axis.text = element_text(colour="black")) +
            scale_fill_manual(values=c("#FFCC00", "#3399FF"))
        ggsave(filename="Repeat enrichment - Absolute bars.png", width=par("din")*1.5)  
    }

    enrich[["TFBS"]] <- NULL
    for(i in 1:length(enrich)) {
        enrich[[i]] <- cbind(enrich[[i]], ColorInfo=ifelse(enrich[[i]][,"Relative Hyper"] < 0, refCondition, otherCondition))
    }
    doCombinedRelativeBarPlot(enrich, c("Relative Hyper", "Relative Hyper Count"), "DMRs", getHyperMethylationLabels(refCondition, otherCondition), TRUE, "ColorInfo")
    doCombinedRelativeBarPlot(enrich, c("Relative DMR", "Relative DMR Count"), "Hypermethylation", getDMREnrichmentLabels(), showCount=FALSE)
}


