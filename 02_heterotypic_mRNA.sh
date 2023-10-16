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



id=gene_cortex_epidermis_pericycle
gene_input_file=data/genes/$id.txt

output=results/02_heterotypic_mRNA/Franco-Zorrilla_et_al_2014/02_heterotypic_promoters_$id
indexoutput=results/01_homotypic_mRNA/Franco-Zorrilla_et_al_2014/01_homotypic_index



>$output/new_genes.txt


while IFS=$'\r\n' read -r line; do
    line=$(echo -n "$line")
    IFS=' ' read -ra parts <<< "$line"
    cluster="${parts[0]}"
    gene="${parts[1]}"
    printf "%s.%d\n" "$line" "1"  >> $output/new_genes.txt

    awk -v pattern="^${gene}\\." -v clus="$cluster" '$0 ~ pattern {print clus, $0}' "$indexoutput/universe_isoform.txt" >> "$output/new_genes.txt"
done < "$gene_input_file"







# ---------------------------- gene name to transcript name ----------------------

if [[ -f "$indexoutput/universe.txt" ]]; then
    print_green "Searching for heterotypic motif hits..."
    mkdir -p $output

    # remove genes not present in pre-computed pmet index
    grep -Ff $indexoutput/universe.txt $output/new_genes.txt > $output/new_genes_temp.txt
    rm $output/new_genes.txt

    scripts/pmetParallel\
        -d .                                    \
        -g $output/new_genes_temp.txt           \
        -i 4                                    \
        -p $indexoutput/promoter_lengths.txt    \
        -b $indexoutput/binomial_thresholds.txt \
        -c $indexoutput/IC.txt                  \
        -f $indexoutput/fimohits                \
        -o $output                              \
        -t 8

    rm $output/new_genes_temp.txt

    cat $output/*.txt > $output/motif_output.txt
    rm $output/temp*.txt

else
    print_red "Nothing found in $indexoutput.\n"
    print_fluorescent_yellow "Please run 01_homotypic_promoters.sh first."
fi
