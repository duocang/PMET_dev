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
# download data

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

################################ Running FIMO ###################################
# check if fimo ready for PMETindex
print_fluorescent_yellow "Checking if fimo result (txt files) ready for PMET index\n"

directory="results/homotypic_promoters/fimo"
mkdir -p $directory

txt_files=$(find "$directory" -name "*.txt")
if [ -n "$txt_files" ]; then
    print_green "Yes, (FIMO result) txt files exist in $directory.\n"
else
    print_red   "   No FIMO result (txt files) found in $directory.\n"
    print_green "Running fimo...\n"

    scripts/cpp_debug_needed/homotypic_promoters.sh \
        -r scripts                                  \
        -o results/homotypic_promoters              \
        -i gene_id=                                 \
        -k 5                                        \
        -n 5000                                     \
        -p 1000                                     \
        -v NoOverlap                                \
        -u Yes                                      \
        -t 4                                        \
        data/homotypic_promoters/genome.fasta       \
        data/homotypic_promoters/anno.gff3          \
        data/homotypic_promoters/motif.meme
fi

########################## Running pmet indexing ##################################
print_green "Running pmet indexing...\n"
# run pmet index
mkdir -p results/homotypic_promoters/fimohits
start=$(date +%s)

scripts/pmetindex                                       \
    -f results/homotypic_promoters/fimo                 \
    -k 5 -n 5000                                        \
    -p results/homotypic_promoters/promoter_lengths.txt \
    -o results/homotypic_promoters/

end=$(date +%s)
time_taken=$((end - start))
print_red "Time taken: $time_taken seconds"

print_fluorescent_yellow "\n\nYou many want to run '02_heterotypic_promoters.sh' now..."
