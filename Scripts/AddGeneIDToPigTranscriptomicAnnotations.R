# Load the GeneSymbol -> ENtrez ID map.
# Obtained from BioDBNet: http://biodbnet.abcc.ncifcrf.gov/tools/orgTaxon.php?org=sus+scrofa&request=orgTaxon&requestType=tools
genesymbolToEntrez <- read.table("Annotations/pig/GeneSymbolToEntrezID.txt", header=TRUE, sep="\t")
genesymbolToEntrez$Gene.ID <- gsub(";.*", "", genesymbolToEntrez$Gene.ID)

# Load the transcriptomic slide annotation obtained from ELMA
annotationTrans <- read.table("Annotations/pig/Trans.annotation.txt"), header=TRUE, sep="\t", quote="", comment.char = "")

# Match gene symbols to gene IDs.
annotationTrans$GeneID=genesymbolToEntrez$Gene.ID[match(annotationTrans$Gene_Symbol, genesymbolToEntrez$Gene.Symbol)]
annotationTrans$GeneID[is.na(annotationTrans$GeneID)] <- ""

# Write the new annotation file.
write.table(annotationTrans, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, file="../pig/Trans.annotation.txt")