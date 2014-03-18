
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
doLimmaAnalysis <- function(targetData, intensityData, foldChangeCutoff, pValueCutoff, refCondition) {
	cond1 <- targetData$Cy5[1]
	cond2 <- targetData$Cy3[1]
	
    refCond <- refCondition
    if(is.na(refCondition) || refCondition=="") {
        refCond <- targetData$Cy5[1]
    }
    
	fitDesign <- modelMatrix(targetData, ref=refCond)
	
    controlWeights <- rep(0, nrow(intensityData$R))
    controlWeights[substr(intensityData$gene$ID, 1, 8)=="EDMA_SPK"] <- 1
    controlWeights[substr(intensityData$gene$ID, 1, 8)=="EDMA_DIG"] <- 0.05
    
	Std_MA_Within <- normalizeWithinArrays(intensityData, method="control", bc.method="none",
                                           controlspots=(substr(intensityData$genes$ID,1,8)=="EDMA_SPK")|(substr(intensityData$genes$ID,1,8)=="EDMA_DIG"),
                                           weights=controlWeights)
	Std_MA_Between <- normalizeBetweenArrays(Std_MA_Within, method="quantile")
	
	fit <- lmFit(Std_MA_Between, design=fitDesign)
	ebayes_fit <- eBayes(fit)
	
	coef <- abs(ebayes_fit$coefficients) > foldChangeCutoff
	pval <- ebayes_fit$p.value < pValueCutoff
	indices <- coef & pval
	indices[is.na(indices)] <- FALSE
	
	diffExpr <- data.frame(ID=ebayes_fit$genes[indices,"ID"], Coef=ebayes_fit$coefficients[indices], PVal=ebayes_fit$p.value[indices])
    return(list(Norm=Std_MA_Between, Fit=ebayes_fit, DiffExpr=diffExpr, NormWithin=Std_MA_Within))
}