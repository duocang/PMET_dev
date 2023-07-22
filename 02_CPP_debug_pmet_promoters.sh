output=results/PMET
gene_input_file=data/gene.txt


mkdir -p $output

# remove genes not present in pre-computed pmet index
grep -Ff results/PMETindex/universe.txt $gene_input_file > $gene_input_file"temp"

scripts/pmet \
    -d . \
    -g $gene_input_file"temp" \
    -i 4 \
    -p results/PMETindex/promoter_lengths.txt \
    -b results/PMETindex/binomial_thresholds.txt \
    -c results/PMETindex/IC.txt \
    -f results/PMETindex/fimohits \
    -o $output

rm $gene_input_file"temp"