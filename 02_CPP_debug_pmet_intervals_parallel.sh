output=results/PMET_intervals_parallel
indexoutput=results/PMETindex_intervals

gene_input_file=data/PMETindex_intervals/intervals.txt


# # remove genes not present in promoter_lengths.txt
# awk -F"\t" '{print $1"\t"}' indexoutput/promoter_lengths.txt > indexoutput/temp_genes_list.txt
# cat indexoutput/temp_genes_list.txt | while read line; do
#     grep $line $gene_input_file
# done > genes/temp_${task}.txt
# rm indexoutput/temp_genes_list.txt

mkdir -p $output


scripts/pmetParallel_linux \
    -d . \
    -g $gene_input_file \
    -i 4 \
    -p $indexoutput/promoter_lengths.txt \
    -b $indexoutput/binomial_thresholds.txt \
    -c $indexoutput/IC.txt \
    -f $indexoutput/fimohits \
    -o $output \
    -t 1

cat $output/*.txt > $output/motif_output.txt
rm $output/temp*.txt


# mkdir -p performance
# mv pmetParallel.prof performance
# cp scripts/pmetParallel performance