#bin/bash/!

for id in $(cat $1);
do
    total_fragments=$(wc -l ${id}-kraken2-out.txt | cut -f 1 -d " ")
    fragments_retained=$(grep -w -m 1 "unclassified" ${id}-kraken2-report.txt | cut -f 2)
    perc_removed=$(printf "%.2f\n" $(echo "scale=4; 100 - ${fragments_retained} / ${total_fragments} * 100" | bc -l) )
    printf "${id}\t${total_fragments}\t${fragments_retained}\t${perc_removed}\n" >> building.tmp
done
cat <( printf "id\ttotal_reads_before\ttotal_reads_after\tpercent_reads_removed\n" ) building.tmp > kraken2-read-removal-summary.tsv
