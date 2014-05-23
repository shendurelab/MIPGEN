# first argument is file of list of genes
# second argument is refGene.txt file
# standard out is BED file of coding exons for listed gene symbols
awk '{if (ARGV[1]==FILENAME) { genes[$0] = 1;} else { if ($13 in genes) {split($10,strt, ","); split($11,stp,",");
for (entry in stp){if (strt[entry]==0) { } else{current_start = strt[entry]; 
current_stop = stp[entry];
if ($7 > stp[entry]) {}
else if ($8 < strt[entry]) {}
else {if ($7 > strt[entry]) {current_start = $7} if ($8 < stp[entry]) {current_stop = $8} print $3 "\t" current_start "\t" current_stop "\t" $13 "/" $2}}} delete strt; delete stp;}}}' $1 $2
