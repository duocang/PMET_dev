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

print_green "Searching for heterotypic motif hits..."

output=results/heterotypic_promoters_single_CPU
gene_input_file=data/gene.txt


mkdir -p $output

# remove genes not present in pre-computed pmet index
grep -Ff results/promoters_homotypic/universe.txt $gene_input_file > $gene_input_file"temp"

scripts/pmet \
    -d . \
    -g $gene_input_file"temp" \
    -i 4 \
    -p results/promoters_homotypic/promoter_lengths.txt \
    -b results/promoters_homotypic/binomial_thresholds.txt \
    -c results/promoters_homotypic/IC.txt \
    -f results/promoters_homotypic/fimohits \
    -o $output

rm $gene_input_file"temp"
