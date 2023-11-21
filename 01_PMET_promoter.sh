#!/bin/bash

function error_exit() {
    echo "ERROR: $1" >&2
    usage
    exit 1
}

print_red(){
    RED='\033[0;31m'
    NC='\033[0m' # No Color
    printf "${RED}$1${NC}\n"
}

print_green(){
    GREEN='\033[0;32m'
    NC='\033[0m' # No Color
    printf "${GREEN}$1${NC}\n"
}

print_orange(){
    ORANGE='\033[0;33m'
    NC='\033[0m' # No Color
    printf "${ORANGE}$1${NC}\n"
}

print_fluorescent_yellow(){
    FLUORESCENT_YELLOW='\033[1;33m'
    NC='\033[0m' # No Color
    printf "${FLUORESCENT_YELLOW}$1${NC}\n"
}

print_white(){
    WHITE='\033[1;37m'
    NC='\033[0m' # No Color
    printf "${WHITE}$1${NC}"
}

# Give execute permission to all users for the file.
chmod a+x scripts/cpp_debug_needed/homotypic_promoters.sh
chmod a+x scripts/gff3sort/gff3sort.pl

########################## Downloading data #######################################
cd data
if [ -f "homotypic_promoters/anno.gff3" ]; then
    echo ""
else
    print_green "Downloading genome and annotation...\n"
    mkdir -p data/homotypic_promoters

    chmod a+x ./fetch_data.sh
    bash ./fetch_data.sh
    mv anno.gff3 homotypic_promoters/
    mv genome.fasta homotypic_promoters/
    rm anno.gff3
    rm genome.fasta
fi
cd ..


########################## Parameters #######################################
print_green "Searching for homotypic motif hits..."
indexoutput=results/01_PMET_promoter/homotypic

genome=data/homotypic_promoters/genome.fasta
anno=data/homotypic_promoters/anno.gff3
meme=data/homotypic_promoters/motif.meme


output=results/01_PMET_promoter/heterotypic
gene_input_file=data/gene.txt

################################ Running FIMO ###################################
directory=$indexoutput/fimo
mkdir -p $indexoutput/fimo

scripts/PMETindex_promoters.sh                  \
    -r scripts                                  \
    -o $indexoutput                             \
    -i gene_id=                                 \
    -k 5                                        \
    -n 5000                                     \
    -p 1000                                     \
    -v NoOverlap                                \
    -u Yes                                      \
    -t 8                                        \
    $genome                                     \
    $anno                                       \
    $meme


rm -rf $indexoutput/fimo

##################################### Heterotypic ##################################
print_green "\n\nSearching for heterotypic motif hits..."
mkdir -p $output

# remove genes not present in pre-computed pmet index
grep -Ff $indexoutput/universe.txt $gene_input_file > $gene_input_file"temp"

scripts/pmetParallel                        \
    -d .                                    \
    -g $gene_input_file"temp"               \
    -i 4                                    \
    -p $indexoutput/promoter_lengths.txt    \
    -b $indexoutput/binomial_thresholds.txt \
    -c $indexoutput/IC.txt                  \
    -f $indexoutput/fimohits                \
    -o $output                              \
    -t 2

cat $output/*.txt > $output/motif_output.txt
rm $output/temp*.txt
rm $gene_input_file"temp"
