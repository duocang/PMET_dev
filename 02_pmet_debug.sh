promlength=1000
maxk=5
topn=5000
task=genes_cycle_me_sig001_pcta07_equal_size

output=results/${task}_${promlength}_${maxk}_${topn}
indexoutput=results/PMETindex_parallel_output



gene_input_file=data/${task}.txt


# # remove genes not present in promoter_lengths.txt
# awk -F"\t" '{print $1"\t"}' indexoutput/promoter_lengths.txt > indexoutput/temp_genes_list.txt
# cat indexoutput/temp_genes_list.txt | while read line; do
#     grep $line $gene_input_file
# done > genes/temp_${task}.txt
# rm indexoutput/temp_genes_list.txt



mkdir -p $output

# remove genes not present in pre-computed pmet index
grep -Ff $indexoutput/universe.txt $gene_input_file > $gene_input_file"temp"

scripts/pmet \
    -d . \
    -g $gene_input_file"temp" \
    -i 4 \
    -p $indexoutput/promoter_lengths.txt \
    -b $indexoutput/binomial_thresholds.txt \
    -c $indexoutput/IC.txt \
    -f $indexoutput/fimohits \
    -o $output

cat $output/*.txt > $output/motif_output.txt
rm $output/temp*.txt


cp  ../../01_PMETDEV-code/PMETdev/ui/heatmap.html $output
cp  ../../01_PMETDEV-code/PMETdev/ui/running.gif  $output


rm $gene_input_file"temp"