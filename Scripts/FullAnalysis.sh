#/bin/bash

EPIGENETIC_NAME=$1
EPIGENETIC_TARGET=$2
EPIGENETIC_FOLDER=$3
TRANSCRIPTOMIC_NAME=$4
TRANSCRIPTOMIC_TARGET=$5
TRANSCRIPTOMIC_FOLDER=$6
REFERENCE_CONDITION=$7
COMBINED_NAME=$8
EPI_FOLDCHANGE_THRESHOLD=$9
EPI_PVALUE_THRESHOLD=${10}
EPI_PVALUE_ADJUSTED=${11}
TRANS_FOLDCHANGE_THRESHOLD=${12}
TRANS_PVALUE_THRESHOLD=${13}
TRANS_PVALUE_ADJUSTED=${14}
VERSION=${15}
EMAIL=${16}

# If version is unspecified or invalid, default to v2.
if [ "$VERSION" != "v1" -a "$VERSION" != "pigv1" -a "$VERSION" != "pigv2" ]
then
    VERSION="v2"
fi

# Select organism.
if [ "$VERSION" == "pigv1" -o "$VERSION" == "pigv2" ]
then
    ORGANISM="pig"
else
    ORGANISM="cow"
fi

# Determine location of output file:
OUTDIR="Results/$EPIGENETIC_NAME"
if [ "$EPIGENETIC_NAME" == "" ]
then
    OUTDIR="Results/$TRANSCRIPTOMIC_NAME/"
fi
mkdir -p "$OUTDIR"
OUTPUTLOG="$OUTDIR"/Log.txt


# Generate both transcriptomic and epigenetic analysis
Rscript Scripts/main.R "$EPIGENETIC_NAME" "$EPIGENETIC_TARGET" "$EPIGENETIC_FOLDER" \
                       "$TRANSCRIPTOMIC_NAME" "$TRANSCRIPTOMIC_TARGET" "$TRANSCRIPTOMIC_FOLDER" \
                       "$REFERENCE_CONDITION" "$COMBINED_NAME" \
                       $EPI_FOLDCHANGE_THRESHOLD $EPI_PVALUE_THRESHOLD $EPI_PVALUE_ADJUSTED\
                       $TRANS_FOLDCHANGE_THRESHOLD $TRANS_PVALUE_THRESHOLD $TRANS_PVALUE_ADJUSTED\
                       $VERSION > "$OUTPUTLOG"
rc=$?
if [[ $rc != 0 ]]
then
    echo "main.R exited with error code "$rc". Aborting." >> "$OUTPUTLOG"
else
    # Generate file for importing data into IPA.
    perl -I /home/efournier/UCSCMirror/lib/ Scripts/IPAInput.pl Annotations/$VERSION/EDMA.Annotation.txt "Results/$EPIGENETIC_NAME/Limma Analysis/LimmaAnalysis.txt" > "Results/$EPIGENETIC_NAME/IPA.txt"
                   
    # Create circos input data structure
    NAME_EPI=`sed "s/\//./g" <(echo "$EPIGENETIC_NAME")`
    CIRCOS_FOLDER=$NAME_EPI

    if [ "$COMBINED_NAME" != "" ]
    then
       NAME_COMBI=`sed "s/\//./g" <(echo "$COMBINED_NAME")`
       CIRCOS_FOLDER=$NAME_COMBI
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

    ./GeneratePlotEpi.sh "$CIRCOS_FOLDER" "$NAME_EPI" "$ORGANISM" >> ../"$OUTPUTLOG" 2>&1
    cp "Output/$NAME_EPI"/*.legend.png "../Results/$EPIGENETIC_NAME/"

    if [ "$COMBINED_NAME" != "" ]
    then
        ./GeneratePlotCombined.sh "$CIRCOS_FOLDER" "$NAME_COMBI" "$ORGANISM" >> ../"$OUTPUTLOG" 2>&1
        cp "Output/$NAME_COMBI/"*.legend.png "../Results/$COMBINED_NAME/"
    fi
    cd ..
fi

# Zip the results of the epigenetic/transcriptomic/combined analysis, depending on which was performed.
if [ "$COMBINED_NAME" != "" ]
then
	zip -r "$CIRCOS_FOLDER".zip "Results/$EPIGENETIC_NAME/" "Results/$COMBINED_NAME/" "Results/$TRANSCRIPTOMIC_NAME/"
elif [ "$EPIGENETIC_NAME" != "" ]
then
	zip -r "$CIRCOS_FOLDER".zip "Results/$EPIGENETIC_NAME/"
else
	zip -r "$CIRCOS_FOLDER".zip "Results/$TRANSCRIPTOMIC_NAME/"
fi

ln -s `pwd`/"$CIRCOS_FOLDER".zip /var/www/bioinfo/html/epigenetics/Analysis/
if [ "${16}" != "" ]
then 
    echo "EMAP has finished processing your data. You can obtain your results at http://emb-bioinfo.fsaa.ulaval.ca/bioinfo/html/epigenetics/Analysis/""$CIRCOS_FOLDER".zip | mail -s "EMAP processing is complete." $EMAIL
fi
