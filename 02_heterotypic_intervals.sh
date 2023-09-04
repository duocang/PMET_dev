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



output=results/heterotypic_intervals
indexoutput=results/homotypic_intervals
gene_input_file=data/homotypic_intervals/intervals.txt

# # remove genes not present in promoter_lengths.txt
# awk -F"\t" '{print $1"\t"}' indexoutput/promoter_lengths.txt > indexoutput/temp_genes_list.txt
# cat indexoutput/temp_genes_list.txt | while read line; do
#     grep $line $gene_input_file
# done > genes/temp_${task}.txt
# rm indexoutput/temp_genes_list.txt

if [[ -f "$indexoutput/promoter_lengths.txt" ]]; then
    print_green "Searching for heterotypic motif hits..."
    mkdir -p $output
    scripts/pmetParallel \
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

else
    print_red "Nothing found in $indexoutput.\n"
    print_fluorescent_yellow "Please run '01_homotypic_intervals.sh' first."
fi