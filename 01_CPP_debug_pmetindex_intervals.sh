
# check if fimo ready for PMETindex
directory="results/PMETindex_intervals/fimo"
mkdir -p $directory

txt_files=$(find "$directory" -name "*.txt")

rm -rf results/PMETindex_intervals/binomial_thresholds.txt
rm -rf  results/fimohits

if [ -n "$txt_files" ]; then
    echo "txt files exist in $directory."
else
    echo "No txt files found in $directory. Running fimo.sh..."
    
    scripts/needed_by_debug_PMETindex_intervals.sh  \
    -r scripts \
    -o results/PMETindex_intervals \
    -n 5000 \
    -k 5 \
    -f 0.05 \
    -t 8 \
    data/PMETindex_intervals/intervals.fa \
    data/PMETindex_intervals/motif.meme
fi


mkdir -p results/PMETindex_intervals/fimohits
# run pmet index
scripts/pmetindex \
    -f results/PMETindex_intervals/fimo \
    -k 5 -n 5000 \
    -p results/PMETindex_intervals/promoter_lengths.txt \
    -o results/PMETindex_intervals/

# # mkdir -p performance
# # mv pmetindex.prof performance
# # cp scripts/pmetindex performance

rm pmetindex.prof

exit 0;