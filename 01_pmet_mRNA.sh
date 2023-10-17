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
find . -type f \( -name "*.sh" -o -name "*.pl" \) -exec chmod a+x {} \;


################################ 1. input parameters ####################################
# homotypic
genomic_element=mRNA
delete_temp=no
homotypic_output=results/mRNA/01_homotypic

# heterotypic
gene_input_file=data/gene.txt
heterotypic_output=results/mRNA/02_heterotypic
################################ input parameters ####################################



########################## 2. Downloading data #######################################
cd data
if [ -f "homotypic_promoters/anno.gff3" ]; then
    print_green "Genome and annotation are ready!"
else
    print_fluorescent_yellow "Downloading genome and annotation...\n"
    mkdir -p data/homotypic_promoters
    bash ./fetch_data.sh
    mv anno.gff3 homotypic_promoters/
    mv genome.fasta homotypic_promoters/
    rm anno.gff3
    rm genome.fasta
fi
cd ..


############################## 3. Running homotypic #################################
print_green "Running homotypic searching...\n"

scripts/PMETindex_mRNA.sh                             \
    -r scripts                                        \
    -o $homotypic_output                              \
    -e $genomic_element                               \
    -i ID=transcript:                                 \
    -k 5                                              \
    -n 5000                                           \
    -p 1000                                           \
    -v NoOverlap                                      \
    -u Yes                                            \
    -t 8                                              \
    -d $delete_temp                                   \
    data/homotypic_promoters/genome.fasta             \
    data/homotypic_promoters/anno.gff3                \
    data/Franco-Zorrilla_et_al_2014.meme


############################ 4. Running heterotypic ###############################
# gene name to transcript name
if [[ -f "$homotypic_output/universe.txt" ]]; then
    mkdir -p $heterotypic_output

    print_green "Preparing transcripts ID in all clusters..."

    >$heterotypic_output/new_genes.txt
    while IFS=$'\r\n' read -r line; do
        line=$(echo -n "$line")
        IFS=' ' read -ra parts <<< "$line"
        cluster="${parts[0]}"
        gene="${parts[1]}"
        printf "%s.%d\n" "$line" "1"  >> $heterotypic_output/new_genes.txt

        awk -v pattern="^${gene}." \
            -v clus="$cluster" '$0 ~ pattern {print clus, $0}' "$homotypic_output/universe_isoform.txt" \
            >> "$heterotypic_output/new_genes.txt"
    done < "$gene_input_file"

    # remove genes not present in pre-computed pmet index
    grep -Ff $homotypic_output/universe.txt $heterotypic_output/new_genes.txt > $heterotypic_output/new_genes_temp.txt
    rm $heterotypic_output/new_genes.txt

    print_green "Searching for heterotypic motif hits..."

    scripts/pmetParallel\
        -d .                                         \
        -g $heterotypic_output/new_genes_temp.txt    \
        -i 4                                         \
        -p $homotypic_output/promoter_lengths.txt    \
        -b $homotypic_output/binomial_thresholds.txt \
        -c $homotypic_output/IC.txt                  \
        -f $homotypic_output/fimohits                \
        -o $heterotypic_output                       \
        -t 8 > $heterotypic_output/pmet.log

    rm $heterotypic_output/new_genes_temp.txt
    # merge pmet result
    cat $heterotypic_output/*.txt > $heterotypic_output/motif_output.txt
    rm $heterotypic_output/temp*.txt

    end_time=$SECONDS
    elapsed_time=$((end_time - start_time))
    days=$((elapsed_time/86400))
    hours=$(( (elapsed_time%86400)/3600 ))
    minutes=$(( (elapsed_time%3600)/60 ))
    seconds=$((elapsed_time%60))
    print_green "Time taken: $days day $hours hour $minutes minute $seconds second"
else
    print_red "Nothing found in $homotypic_output.\n"
fi
