
mkdir -p data
mkdir -p results/PMETindex


# download data
cd data
if [ -f "anno.gff3" ]; then
    echo "anno.gff3 exists."
else
    echo "anno.gff3 does not exist. Fetching data..."
    ./fetch_data.sh
fi
cd ..


# check if fimo ready for PMETindex
directory="results/PMETindex/fimo"
txt_files=$(find "$directory" -name "*.txt")

if [ -n "$txt_files" ]; then
    echo "txt files exist in $directory."
else
    echo "No txt files found in $directory. Running fimo.sh..."
    
    scripts/needed_by_PMETindex.sh \
    -r scripts \
    -o results/PMETindex \
    -i gene_id= \
    -k 5 \
    -n 5000 \
    -p 1000 \
    -v NoOverlap \
    -u Yes \
    -t 16 \
    data/genome.fasta \
    data/anno.gff3 \
    data/motif.meme
fi



# run pmet index
scripts/pmetindex \
    -f results/PMETindex/fimo \
    -k 5 -n 5000 \
    -p results/PMETindex/promoter_lengths.txt \
    -o results/PMETindex/

mkdir -p performance
mv pmetindex.prof performance
cp scripts/pmetindex performance
