# first argument is file of list of genes
awk '{if (FILENAME==ARGV[1]) { genes[$0] = 1;} else { if ($13 in genes) {split($10,strt, ","); split($11,stp,",");for (entry in stp) {if (length(strt[entry])==0) { } else{print $3 "\t" strt[entry] - 1 "\t" stp[entry] "\t" $13 "/" $2}}}}}' $1 refGene.txt
