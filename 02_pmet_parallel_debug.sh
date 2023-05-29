output=results/PMET_parallel
indexoutput=results/PMETindex

gene_input_file=data/gene.txt


# # remove genes not present in promoter_lengths.txt
# awk -F"\t" '{print $1"\t"}' indexoutput/promoter_lengths.txt > indexoutput/temp_genes_list.txt
# cat indexoutput/temp_genes_list.txt | while read line; do
#     grep $line $gene_input_file
# done > genes/temp_${task}.txt
# rm indexoutput/temp_genes_list.txt

mkdir -p $output

# remove genes not present in pre-computed pmet index
grep -Ff $indexoutput/universe.txt $gene_input_file > $gene_input_file"temp"

scripts/pmetParallel \
    -d . \
    -g $gene_input_file"temp" \
    -i 4 \
    -p $indexoutput/promoter_lengths.txt \
    -b $indexoutput/binomial_thresholds.txt \
    -c $indexoutput/IC.txt \
    -f $indexoutput/fimohits \
    -o $output \
    -t 2

cat $output/*.txt > $output/motif_output.txt
rm $output/temp*.txt


rm $gene_input_file"temp"

mkdir -p performance
mv pmetParallel.prof performance
cp scripts/pmetParallel performance