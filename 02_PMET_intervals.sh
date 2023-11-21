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
chmod a+x scripts/cpp_debug_needed/homotypic_intervals.sh

##################################### Parameters #######################################
indexOutput=results/02_PMET_intervals/homotypic

output=results/02_PMET_intervals/heterotypic
gene_input_file=data/homotypic_intervals/intervals.txt

##################################### Running FIMO #####################################
print_green "Running fimo...\n"

scripts/cpp_debug_needed/homotypic_intervals.sh  \
    -r scripts                                   \
    -o $indexOutput                              \
    -n 5000                                      \
    -k 5                                         \
    -f 0.05                                      \
    -t 8                                         \
    data/homotypic_intervals/intervals.fa        \
    data/homotypic_intervals/motif_more.meme



################################### Heterotypic #####################################

# # remove genes not present in promoter_lengths.txt
# awk -F"\t" '{print $1"\t"}' indexOutput/promoter_lengths.txt > indexOutput/temp_genes_list.txt
# cat indexOutput/temp_genes_list.txt | while read line; do
#     grep $line $gene_input_file
# done > genes/temp_${task}.txt
# rm indexOutput/temp_genes_list.txt

print_green "Searching for heterotypic motif hits..."
mkdir -p $output
scripts/pmetParallel                        \
    -d .                                    \
    -g $gene_input_file                     \
    -i 4                                    \
    -p $indexOutput/promoter_lengths.txt    \
    -b $indexOutput/binomial_thresholds.txt \
    -c $indexOutput/IC.txt                  \
    -f $indexOutput/fimohits                \
    -o $output                              \
    -t 1

cat $output/*.txt > $output/motif_output.txt
rm $output/temp*.txt
