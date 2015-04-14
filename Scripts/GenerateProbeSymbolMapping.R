annotation <- read.table("EDMA.Annotation.txt", header=TRUE, sep="\t")

finalIDs <-rep("", nrow(annotation))
for(category in c("Exon", "Intron", "Promoter", "Proximal_Promoter")) {
    trimIDs <- gsub("^\\s+", "", annotation[,category])
    rawIDs <- strsplit(x=as.character(trimIDs), split=" ", fixed=TRUE)
    catIDs <- lapply(rawIDs, function(x) { if(length(x)>0) { return(paste(category, x, sep="-", collapse=",")) } else {return(character(0)) } })
    combinedIDs <- paste(finalIDs, catIDs, sep=",")
    combinedIDs <- gsub(",character\\(0\\)$", "", combinedIDs)
    finalIDs <- gsub("^,", "", combinedIDs)
}

write.table(data.frame(Probe=annotation$Probe, Gene_Symbol=finalIDs), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, file="ProbeSymbol.mapping")
