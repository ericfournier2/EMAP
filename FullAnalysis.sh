#/bin/bash

EPIGENETIC_NAME=$1
EPIGENETIC_TARGET=$2
EPIGENETIC_FOLDER=$3
TRANSCRIPTOMIC_NAME=$4
TRANSCRIPTOMIC_TARGET=$5
TRANSCRIPTOMIC_FOLDER=$6
REFERENCE_CONDITION=$7
COMBINED_NAME=$8
FOLDCHANGE_THRESHOLD=$9
PVALUE_THRESHOLD=${10}
VERSION=${11}

# If version is unspecified or invalid, default to v2.
if [ "$VERSION" != "v1" ]
then
    VERSION="v2"
fi

# Generate both transcriptomic and epigenetic analysis
Rscript main.R "$EPIGENETIC_NAME" "$EPIGENETIC_TARGET" "$EPIGENETIC_FOLDER" \
               "$TRANSCRIPTOMIC_NAME" "$TRANSCRIPTOMIC_TARGET" "$TRANSCRIPTOMIC_FOLDER" \
               "$REFERENCE_CONDITION" "$COMBINED_NAME" $FOLDCHANGE_THRESHOLD $PVALUE_THRESHOLD $VERSION

perl IPAInput.pl Annotations/$VERSION/EDMA.Annotation.txt "Results/$EPIGENETIC_NAME/Limma Analysis/LimmaAnalysis.txt" > "Results/$EPIGENETIC_NAME/IPA.txt"
               
# Create circos input data structure
CIRCOS_FOLDER=`sed "s/\//./g" <(echo "$EPIGENETIC_NAME")`

if [ "$COMBINED_NAME" != "" ]
then
   CIRCOS_FOLDER=`sed "s/\//./g" <(echo "$COMBINED_NAME")`
   mkdir -p "CircosAutomation/Data/$CIRCOS_FOLDER"
   
   ln -f -s "`pwd`/Results/$TRANSCRIPTOMIC_NAME/Trans-FC.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
   ln -f -s "`pwd`/Results/$TRANSCRIPTOMIC_NAME/Trans-P-Value.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
   ln -f -s "`pwd`/Results/$TRANSCRIPTOMIC_NAME/DiffExpr-P-Value-Trans.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
   ln -f -s "`pwd`/Results/$COMBINED_NAME/Concordant.txt" "CircosAutomation/Data/$CIRCOS_FOLDER"
fi

mkdir -p "CircosAutomation/Data/$CIRCOS_FOLDER"
ln -f -s "`pwd`/Results/$EPIGENETIC_NAME/Bedgraph/DiffExpr-P-Value.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
ln -f -s "`pwd`/Results/$EPIGENETIC_NAME/Bedgraph/Fragment-Cond-Mean-A.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
ln -f -s "`pwd`/Results/$EPIGENETIC_NAME/Bedgraph/Fragment-Cond-Mean-B.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
ln -f -s "`pwd`/Results/$EPIGENETIC_NAME/Bedgraph/Fragment-Fold-Change.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
ln -f -s "`pwd`/Results/$EPIGENETIC_NAME/Bedgraph/Fragment-P-Value.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
ln -f -s "`pwd`/Results/$EPIGENETIC_NAME/Bedgraph/Fragment-Fold-Change-Significant.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
ln -f -s "`pwd`/Results/$EPIGENETIC_NAME/Bedgraph/Fragment-P-Value-Significant.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
ln -f -s "`pwd`/Results/$EPIGENETIC_NAME/Bedgraph/Fragment-Fold-Change-30K.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
ln -f -s "`pwd`/Results/$EPIGENETIC_NAME/Bedgraph/Fragment-P-Value-30K.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"
ln -f -s "`pwd`/Results/$EPIGENETIC_NAME/Bedgraph/Imprint-Fold-Change.bedgraph" "CircosAutomation/Data/$CIRCOS_FOLDER"

cd CircosAutomation
if [ "$COMBINED_NAME" == "" ]
then
    ./GeneratePlotEpi.sh "$CIRCOS_FOLDER"
    cp "Output/$CIRCOS_FOLDER"/*.legend.png "../Results/$EPIGENETIC_NAME/"
else
    ./GeneratePlotCombined.sh "$CIRCOS_FOLDER"
    cp "Output/$CIRCOS_FOLDER/"*.legend.png "../Results/$COMBINED_NAME/"
fi
cd ..


