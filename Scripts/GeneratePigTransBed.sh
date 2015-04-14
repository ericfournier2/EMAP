awk '{ print $14 "\t" $16 "\t" $17 "\t" $10 }' ~/Projects/Pig/Probes-Final.psl > probes.bed
cut -f 4 probes.bed | sort | uniq -d > duplicated.bed
grep -v -F -f duplicated.bed probes.bed > Annotations/pig/Trans.bed

rm probes.bed duplicated.bed
