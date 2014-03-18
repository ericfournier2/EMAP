rm -r ProcessedData
mkdir -p ProcessedData

# Process all fragment bedgraphs.
for f in Data/"$1"/Fragment-*
do
    outbase=`basename "$f"`
    sed "1d" "$f" | grep -v -F "NA" | grep -F "chr29" | perl AverageWindow.pl 100000 > "ProcessedData/$outbase"
done

sed "1d" "Data/$1/Imprint-Fold-Change.bedgraph" | grep -v -F "NA" | grep -F "chr29" | perl AverageWindowImprint.pl 100000 > "ProcessedData/Imprint-Fold-Change.bedgraph"


# Process DiffExpr values for scatterplot.
sed "1d" "Data/$1/DiffExpr-P-Value.bedgraph" | grep -F "chr29" | sort -k 4nr | head -n 100 > ProcessedData/DiffExpr-P-Value.bedgraph

# Create output directory.
mkdir -p Output/"$1"

# Replace output file name within configuration file, and make a copy in the output directory.
sed -e "s/OUTPUT_FILE/$1/g" Control-CSF2-Chr29.conf > "Output/$1/$1.epi.circos.conf"

# Generate The ciscos plot.
circos -conf "Output/$1/$1.epi.circos.conf"

# Add the legend through composition.
composite -gravity center Input/Epi_Legend.png "Output/$1/$1.png" "Output/$1/$1.legend.png"