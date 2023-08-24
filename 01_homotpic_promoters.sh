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


########################## Running FIMO #######################################
# download data
mkdir -p data/promoters_homotypic
cd data
if [ -f "promoters_homotypic/anno.gff3" ]; then
    echo "anno.gff3 exists."
else
    echo "anno.gff3 does not exist. Fetching data..."
    ./fetch_data.sh
    mv anno.gff3 promoters_homotypic/
    mv genome.fasta promoters_homotypic/
    rm anno.gff3
    rm genome.fasta
fi
cd ..

################################ Running FIMO ###################################
# check if fimo ready for PMETindex
print_orange "check if fimo ready for PMETindex"

directory="results/promoters_homotypic/fimo"
mkdir -p $directory

txt_files=$(find "$directory" -name "*.txt")
if [ -n "$txt_files" ]; then
    print_green "(FIMO result) txt files exist in $directory."
else
    print_red "(FIMO result) no txt files found in $directory."
    print_fluorescent_yellow "Running fimo.sh..."

    scripts/cpp_debug_needed/homotypic_promoters.sh \
        -r scripts \
        -o results/promoters_homotypic \
        -i gene_id= \
        -k 5 \
        -n 5000 \
        -p 1000 \
        -v NoOverlap \
        -u Yes \
        -t 4 \
        data/promoters_homotypic/genome.fasta \
        data/promoters_homotypic/anno.gff3 \
        data/promoters_homotypic/motif.meme
fi

########################## Running pmet indexing ##################################
print_green "Running pmet indexing..."
# run pmet index
mkdir -p results/promoters_homotypic/fimohits

scripts/pmetindex \
    -f results/promoters_homotypic/fimo \
    -k 5 -n 5000 \
    -p results/promoters_homotypic/promoter_lengths.txt \
    -o results/promoters_homotypic/

print_green "done"
