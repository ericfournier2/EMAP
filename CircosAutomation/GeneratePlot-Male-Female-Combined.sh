rm -r ProcessedData
mkdir -p ProcessedData

# Process all fragment bedgraphs.
for f in Data/"$1"/*-Fragment-*
do
    outbase=`basename "$f"`
    sed "1d" "$f" | grep -v -F "NA" | perl AverageWindow.pl 5000000 > "ProcessedData/$outbase"
done

# Process DiffExpr values for scatterplot.
sed "1d" "Data/$1/Male-DiffExpr-P-Value.bedgraph" | sort -k 4nr | head -n 100 > ProcessedData/Male-DiffExpr-P-Value.bedgraph
sed "1d" "Data/$1/Female-DiffExpr-P-Value.bedgraph" | sort -k 4nr | head -n 100 > ProcessedData/Female-DiffExpr-P-Value.bedgraph

# Create output directory.
mkdir -p Output/"$1"

# Replace output file name within configuration file, and make a copy in the output directory.
sed -e "s/OUTPUT_FILE/$1/g" Male-Female-Combined.conf > "Output/$1/$1.circos.conf"

# Generate The ciscos plot.
circos -conf "Output/$1/$1.circos.conf"

# Add the legend through composition.
composite -gravity center Input/Epi_Trans_Legend.png "Output/$1/$1.png" "Output/$1/$1.legend.png"