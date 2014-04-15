rm -r ProcessedData
mkdir -p ProcessedData

# Process all fragment bedgraphs.
for f in Data/"$1"/Fragment-* 
do
    outbase=`basename "$f"`
    sed "1d" "$f" | grep -v -F "NA" | perl AverageWindow.pl 5000000 20 > "ProcessedData/$outbase"
done

for f in Data/"$1"/Fragment-*-Significant* Data/"$1"/Trans-*
do
    outbase=`basename "$f"`
    sed "1d" "$f" | grep -v -F "NA" | perl AverageWindow.pl 5000000 1 > "ProcessedData/$outbase"
done

for f in Data/"$1"/Fragment-*-30K*
do
    outbase=`basename "$f"`
    sed "1d" "$f" | grep -v -F "NA" | perl AverageWindow.pl 5000000 5 > "ProcessedData/$outbase"
done

# Process DiffExpr values for scatterplot.
sed "1d" "Data/$1/DiffExpr-P-Value.bedgraph" | sort -k 4nr | head -n 200 > ProcessedData/DiffExpr-P-Value.bedgraph
sed "1d" "Data/$1/DiffExpr-P-Value-Trans.bedgraph" | sort -k 4nr | head -n 200 > ProcessedData/DiffExpr-P-Value-Trans.bedgraph

# Copy "Concordant" text file.
cp Data/"$1"/Concordant.txt ProcessedData/Concordant.txt

# Create output directory.
mkdir -p Output/"$2"

# Replace output file name within configuration file, and make a copy in the output directory.
sed -e "s/OUTPUT_FILE/$2/g" Input/conf/circos.conf.template > temp
sed -e "s/PROBE_SUBSET//g" temp > "Output/$2/$2.circos.conf"
sed -e "s/PROBE_SUBSET/-Significant/g" temp > "Output/$2/$2.circos.sig.conf"
sed -e "s/PROBE_SUBSET/-30K/g" temp > "Output/$2/$2.circos.30K.conf"
rm temp

# Generate The ciscos plot.
circos -conf "Output/$2/$2.circos.conf"
circos -conf "Output/$2/$2.circos.sig.conf"
circos -conf "Output/$2/$2.circos.30K.conf"

# Add the legend through composition.
composite -gravity center Input/Epi_Trans_Legend.png "Output/$2/$2.png" "Output/$2/$2.legend.png"
composite -gravity center Input/Epi_Trans_Legend_Significant.png "Output/$2/$2-Significant.png" "Output/$2/$2-Significant.legend.png"
composite -gravity center Input/Epi_Trans_Legend_Significant.png "Output/$2/$2-30K.png" "Output/$2/$2-30K.legend.png"