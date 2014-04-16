# Function to give a 0 weight to the data which was marked as "Ignored" in the intensity files.
# Zeroes will be replaced with NAs after data is loaded.
# If your intensity file use a marker different from "Ignored cells" to mark cells which should
# be ignored, change the name of the column in this function.
weightFunction <- function(intensityInput) {
    result <- intensityInput[,"Ignored cells"]
    zeroIndices <- intensityInput[,"Ignored cells"]==0 
    result[zeroIndices] <- 1
    result[!zeroIndices] <- 0
     
    return(result)
}		

# Column names for the data within the files. If your columns have some other names, change these lines.
columnNames <- c(R="Raw intensity (mean) {635}",	# Red channel intensity
                 G="Raw intensity (mean) {532}",	# Green channel intensity
                 Rb="Background (mean) {635}",		# Red background intensity
                 Gb="Background (mean) {532}") 		# Green background intensity

# Load a set of data, regardless of whether it is transcriptomic or epigenetic data.
loadData <- function(dataPath, targetName) {
    # Read in sample information.
    targets <- readTargets(file.path(dataPath, targetName))

    # Read in data.
    intensityData <- read.maimages(files=targets$Filename, path=dataPath, source="generic", sep="\t",
                                   annotation=c("ID","Ignored cells"),	columns=columnNames, wt.fun=weightFunction)

    # Replace ignored cells (zeroes) with NAs.					   
    for(i in 1:ncol(intensityData$R)) {
        indices <- intensityData$weights[,i]==0
        intensityData$R[indices,i] <- NA
        intensityData$G[indices,i] <- NA
    }

    return(list(Target=targets, IntensityData=intensityData))
}

# Detects cases where the entries within the input files are unaligned,
# which will cause limma to assign intensity values to the wrong probes.
performSanityCheck <- function(dataObject, dataFolder) {
    # Sanity check
    allData <- list()
    for(i in dataObject$Target$Filename) {
        allData[[i]] <- read.table(file.path(dataFolder, i), header=TRUE, sep="\t")
    }

    hasError <- FALSE
    for(i in 2:length(allData)) {
        if(any(allData[[1]]$ID!=allData[[i]]$ID, na.rm=TRUE)) {
            write(paste("Mismatch between", names(allData)[1], "and", names(allData)[i]), stderr())
            hasError <- TRUE
        }
    }
    
    return(!hasError);
}
