
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
doLimmaAnalysisControl <- function(targetData, intensityData, foldChangeCutoff, pValueCutoff, refCondition) {
    # Create a "marker" file to make it clear that the generated results
    # have had control normalization applied.
    file.create("WARNING CONTROL NORMALIZATION IS BEING USED", showWarnings=FALSE)
    
    # If no reference condition was provided, pick one randomly (Red channel of first array).
    refCond <- refCondition
    if(is.na(refCondition) || refCondition=="") {
        refCond <- targetData$Cy5[1]
    }
    
    # Determine which probes are above the background.
    aboveBG <- listAboveBG(intensityData, targetData, refCond)
    
    # Perform normalization
    controlWeights <- rep(0, nrow(intensityData$R))
    controlWeights[substr(intensityData$gene$ID, 1, 8)=="EDMA_SPK"] <- 1
    controlWeights[substr(intensityData$gene$ID, 1, 8)=="EDMA_DIG"] <- 0.05
    
	Std_MA_Within <- normalizeWithinArrays(intensityData, method="control", bc.method="none",
                                           controlspots=(substr(intensityData$genes$ID,1,8)=="EDMA_SPK")|(substr(intensityData$genes$ID,1,8)=="EDMA_DIG"),
                                           weights=controlWeights)
	Std_MA_Between <- normalizeBetweenArrays(Std_MA_Within, method="quantile")

    # Perform linear fit and bayesian correction.
	fitDesign <- modelMatrix(targetData, ref=refCond)	
	fit <- lmFit(Std_MA_Between, design=fitDesign)
	ebayes_fit <- eBayes(fit)
	
    # Determine which probes are differentially expressed/methylated.
	coef <- abs(ebayes_fit$coefficients) > foldChangeCutoff
	pval <- ebayes_fit$p.value < pValueCutoff
	indices <- coef & pval
	indices[is.na(indices)] <- FALSE
    diffExpr <- subsetFit(ebayes_fit, indices)
    
    # Build the return object.
    return(list(Norm=Std_MA_Between, Fit=ebayes_fit, DiffExpr=diffExpr, AboveBG=aboveBG))    
}