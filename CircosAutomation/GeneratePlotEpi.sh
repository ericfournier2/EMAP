rm -r ProcessedData
mkdir -p ProcessedData

# Process all fragment bedgraphs.
for f in Data/"$1"/Fragment-*
do
    outbase=`basename "$f"`
    sed "1d" "$f" | grep -v -F "NA" | grep -v "^gi" | perl AverageWindow.pl 5000000 20 > "ProcessedData/$outbase"
done

for f in Data/"$1"/Fragment-*-Significant*
do
    outbase=`basename "$f"`
    sed "1d" "$f" | grep -v -F "NA" | grep -v "^gi" | perl AverageWindow.pl 5000000 1 > "ProcessedData/$outbase"
done

for f in Data/"$1"/Fragment-*-30K*
do
    outbase=`basename "$f"`
    sed "1d" "$f" | grep -v -F "NA" | grep -v "^gi" | perl AverageWindow.pl 5000000 5 > "ProcessedData/$outbase"
done

sed "1d" "Data/$1/Imprint-Fold-Change.bedgraph" | grep -v -F "NA" | grep -v "^gi" | perl AverageWindowImprint.pl 5000000 > "ProcessedData/Imprint-Fold-Change.bedgraph"


# Process DiffExpr values for scatterplot.
sed "1d" "Data/$1/DiffExpr-P-Value.bedgraph" | grep -v "^gi" | sort -k 4nr | head -n 200 > ProcessedData/DiffExpr-P-Value.bedgraph

# Create output directory.
mkdir -p Output/"$2"

sed -e "s/OUTPUT_FILE/$2/g" Input/conf/epi.circos.conf.template |
    sed -e "s/ORGANISM/$3/g"  > Output/temp
sed -e "s/PROBE_SUBSET//g" Output/temp > "Output/$2/$2.epi.circos.conf"
sed -e "s/PROBE_SUBSET/-Significant/g" Output/temp > "Output/$2/$2.epi.circos.sig.conf"
sed -e "s/PROBE_SUBSET/-30K/g" Output/temp > "Output/$2/$2.epi.circos.30K.conf"
rm Output/temp

# Generate The circos plot.
circos -conf "Output/$2/$2.epi.circos.conf"
circos -conf "Output/$2/$2.epi.circos.sig.conf"
circos -conf "Output/$2/$2.epi.circos.30K.conf"

# Add the legend through composition.
composite -gravity center Input/Epi_Legend.png "Output/$2/$2.png" "Output/$2/$2.legend.png"
composite -gravity center Input/Epi_Legend_Significant.png "Output/$2/$2-Significant.png" "Output/$2/$2-Significant.legend.png"
composite -gravity center Input/Epi_Legend_Significant.png "Output/$2/$2-30K.png" "Output/$2/$2-30K.legend.png"