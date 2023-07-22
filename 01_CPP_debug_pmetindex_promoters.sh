# download data
mkdir -p data/PMETindex_promoters
cd data
if [ -f "anno.gff3" ]; then
    echo "anno.gff3 exists."
else
    echo "anno.gff3 does not exist. Fetching data..."
    ./fetch_data.sh
    mv anno.gff3 PMETindex_promoters/
    mv genome.fasta PMETindex_promoters/
    rm anno.gff3
    rm genome.fasta
fi
cd ..


# check if fimo ready for PMETindex
directory="results/PMETindex/fimo"

mkdir -p $directory

txt_files=$(find "$directory" -name "*.txt")

if [ -n "$txt_files" ]; then
    echo "txt files exist in $directory."
else
    echo "No txt files found in $directory. Running fimo.sh..."

    scripts/cpp_debug_needed/needed_by_PMETindex.sh \
    -r scripts \
    -o results/PMETindex_promoters \
    -i gene_id= \
    -k 5 \
    -n 5000 \
    -p 1000 \
    -v NoOverlap \
    -u Yes \
    -t 4 \
    data/PMETindex_promoters/genome.fasta \
    data/PMETindex_promoters/anno.gff3 \
    data/PMETindex_promoters/motif.meme
fi

# run pmet index
mkdir -p results/PMETindex_promoters/fimohits

scripts/pmetindex \
    -f results/PMETindex_promoters/fimo \
    -k 5 -n 5000 \
    -p results/PMETindex_promoters/promoter_lengths.txt \
    -o results/PMETindex_promoters/

# # mkdir -p performance
# # mv pmetindex.prof performance
# # cp scripts/pmetindex performance
# rm pmetindex.prof

exit 0;