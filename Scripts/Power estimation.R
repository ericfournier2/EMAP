library(limma)
library(ggplot2)

# Set seed to make results reproducible
set.seed(87619648)

# Load intensity values, IDs and fold-changes.
masterR <- read.table("DessieR.txt", sep="\t", header=TRUE)
masterG <- read.table("DessieG.txt", sep="\t", header=TRUE)
masterIDs = scan("DessieIDs.txt", what=character())
masterFC <- read.table("MasterFC.txt", sep="\t", header=TRUE)

# Estimate fold-change SD.
fcSD <- sd(unlist(masterFC), na.rm=TRUE)

# Build dummy target data.frame.
targetData <- data.frame(Cy3=rep(c("Treatment", "Control"), c(7, 7)), Cy5=rep(c("Control", "Treatment"), c(7, 7)))

# Calculate mean intensity of all probes for use later when simulating fold-changes.
meanIntensities <- apply(cbind(masterR, masterG), 1, mean)

resultsDF <- data.frame(N=integer(0),                           # Number of replicate in simulation run.
                        NumDMRs=numeric(0),                     # Number of DMR regions identified.
                        PercFalsePositive=numeric(0),           # Percentage of false positives within DMRs
                        PercTruePositiveIdentified=numeric(0),  # Percentage of simulated "true positive" which were identified.
                        Top100=numeric(0),                      # Percentage of true positives in Top 100 p-values.
                        Top500=numeric(0),                      # Percentage of true positives in Top 500 p-values.
                        Top1000=numeric(0),                     # Percentage of true positives in Top 1000 p-values.
                        Top2000=numeric(0),                     # Percentage of true positives in Top 2000 p-values.
                        Top4000=numeric(0),                      # Percentage of true positives in Top 4000 p-values.
                        Top8000=numeric(0))                      # Percentage of true positives in Top 8000 p-values.

# Perform simulation a 100 times.
for(i in 1:100) {                        
    # Use null data set as a starting point.
    rValues <- masterR
    gValues <- masterG
    
    # Generate a set of "true" biological fold-change by sampling past fold-changes.
    trueFC <- sample(unlist(masterFC), nrow(rValues))
    
    # Add variability to the "true" fold-change to model individual x treatment interaction.
    actualFC <- t(apply(matrix(trueFC, ncol=1), 1, function(x) { rnorm(14, x, fcSD) }))
    
    # Transform the simulated FC into absolute intensity values offsets.
    addIntensity <- ((2^actualFC)-1)*meanIntensities
    
    # Turn simulated fold-changes into a condition-dependant matrix so it can be added to the matrix of all raw values..
    rFC <- cbind(matrix(0, ncol=7, nrow=nrow(rValues)), addIntensity[,1:7])
    gFC <- cbind(addIntensity[,8:14], matrix(0, ncol=7, nrow=nrow(rValues)))
    
    # Add raw values and simulated fold-changes.
    rValues <- rValues + rFC
    gValues <- gValues + gFC

    # Cap intensity values to the maximum value of the scanner.
    rValues[rValues > 65535] <- 65535
    gValues[gValues > 65535] <- 65535
    
    # Create an RGList object for limma processing.
    rgObject <- new("RGList", list(R=rValues, G=gValues))
    
    # Normalize within arrays outside of the loop, since the results should not be affected by the number of replicates.
    normWithin <- normalizeWithinArrays(rgObject, method="loess")
    
    for(nRep in c(4, 6, 8, 10, 12, 14)) {
        # Define a set of samples to use dependant on the number of replicates.
        colIndices <- c(1:(nRep/2), 7:(6+(nRep/2)))
        
        # Normalize between arrays and apply linear fit.
        normBetween <- normalizeBetweenArrays(normWithin[,colIndices], method="quantile")
        fitResults <- lmFit(normBetween, modelMatrix(targetData[colIndices,], ref="Control"))
        ebayesFit <- eBayes(fitResults)
        
        # Determine which probes are DMRs, and how many there are.
        dmProbes <- abs(ebayesFit$coefficients) > log2(1.5) & (ebayesFit$p.value < 0.05)
        numDMRs <- sum(dmProbes, na.rm=TRUE)
        
        # Determine which probes yielded true and false positives.
        expectedPositives <- abs(trueFC)>log2(1.5)
        truePositives <- dmProbes & expectedPositives
        falsePositives <- dmProbes & !expectedPositives
        
        percFalsePositives <- sum(falsePositives, na.rm=TRUE)/numDMRs
        percTruePositivesFound <- sum(truePositives, na.rm=TRUE)/sum(expectedPositives, na.rm=TRUE)
        
        pOrder <- order(ebayesFit$p.value)
        
        getTopXTruePositivePerc <- function(x) {
            useInd <- pOrder[1:x]
            dmr <- abs(ebayesFit$coefficients[useInd]) > log2(1.5) & (ebayesFit$p.value[useInd] < 0.05)
            truePos <- abs(trueFC[useInd])>log2(1.5)
            
            return(sum(truePos, na.rm=TRUE)/sum(dmr, na.rm=TRUE))
        }
        
        top100 <- getTopXTruePositivePerc(100)
        top500 <- getTopXTruePositivePerc(500)
        top1000 <- getTopXTruePositivePerc(1000)
        top2000 <-  getTopXTruePositivePerc(2000)
        top4000 <-  getTopXTruePositivePerc(4000)
        top8000 <-  getTopXTruePositivePerc(8000)
        
        
        resultRow <- data.frame(N=nRep,
                                NumDMRs=numDMRs,
                                PercFalsePositive=percFalsePositives,
                                PercTruePositiveIdentified=percTruePositivesFound,
                                PercTruePositiveTop100=top100,
                                PercTruePositiveTop500=top500,
                                PercTruePositiveTop1000=top1000,
                                PercTruePositiveTop2000=top2000,
                                PercTruePositiveTop4000=top4000,
                                PercTruePositiveTop8000=top8000)
                                
        resultsDF <- rbind(resultsDF, resultRow)
        write.table(resultsDF, file="SimulationResults.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
    }
}
