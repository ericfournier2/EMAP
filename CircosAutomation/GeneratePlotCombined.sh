#rm -r ProcessedData
mkdir -p ProcessedData/$2

# Process all fragment bedgraphs.
for f in Data/"$1"/Fragment-* 
do
    outbase=`basename "$f"`
    sed "1d" "$f" | grep -v -F "NA" | grep -v "^gi" | perl AverageWindow.pl 5000000 20 > "ProcessedData/$2/$outbase"
done

for f in Data/"$1"/Fragment-*-Significant* Data/"$1"/Trans-*
do
    outbase=`basename "$f"`
    sed "1d" "$f" | grep -v -F "NA" | grep -v "^gi" | perl AverageWindow.pl 5000000 1 > "ProcessedData/$2/$outbase"
done

for f in Data/"$1"/Fragment-*-30K*
do
    outbase=`basename "$f"`
    sed "1d" "$f" | grep -v -F "NA" | grep -v "^gi" | perl AverageWindow.pl 5000000 5 > "ProcessedData/$2/$outbase"
done

# Process DiffExpr values for scatterplot.
sed "1d" "Data/$1/DiffExpr-P-Value.bedgraph" | grep -v "^gi" | sort -k 4nr | head -n 200 > ProcessedData/$2/DiffExpr-P-Value.bedgraph
sed "1d" "Data/$1/DiffExpr-P-Value-Trans.bedgraph" | grep -v "^gi" | sort -k 4nr | head -n 200 > ProcessedData/$2/DiffExpr-P-Value-Trans.bedgraph

# Copy "Concordant" text file.
cp Data/"$1"/Concordant.txt ProcessedData/$2/Concordant.txt

# Create output directory.
mkdir -p Output/"$2"

# Replace output file name within configuration file, and make a copy in the output directory.
sed -e "s/OUTPUT_FILE/$2/g" Input/conf/circos.conf.template |
    sed -e "s/ORGANISM/$3/g"  > Output/temp
sed -e "s/PROBE_SUBSET//g" Output/temp > "Output/$2/$2.circos.conf"
sed -e "s/PROBE_SUBSET/-Significant/g" Output/temp > "Output/$2/$2.circos.sig.conf"
sed -e "s/PROBE_SUBSET/-30K/g" Output/temp > "Output/$2/$2.circos.30K.conf"
rm Output/temp

# Generate The ciscos plot.
circos -silent -conf "Output/$2/$2.circos.conf"
circos -silent -conf "Output/$2/$2.circos.sig.conf"
circos -silent -conf "Output/$2/$2.circos.30K.conf"

# Add the legend through composition.
composite -gravity center Input/Epi_Trans_Legend.png "Output/$2/$2.png" "Output/$2/$2.legend.png"
composite -gravity center Input/Epi_Trans_Legend_Significant.png "Output/$2/$2-Significant.png" "Output/$2/$2-Significant.legend.png"
composite -gravity center Input/Epi_Trans_Legend_Significant.png "Output/$2/$2-30K.png" "Output/$2/$2-30K.legend.png"