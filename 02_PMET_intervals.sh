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
# tool
toolDir=scripts
HOMOTYPIC=$toolDir/cpp_debug_needed/homotypic_intervals.sh
HETEROTYPIC=$toolDir/pmetParallel

chmod a+x $HOMOTYPIC
chmod a+x $HETEROTYPIC

threads=24
res_dir=results/02_PMET_intervals

# homotypic
gff3id="gene_id="
overlap="NoOverlap"
utr="No"
topn=5000
maxk=5
length=1000
fimothresh=0.05
distance=1000
gff3id="gene_id="
delete_temp=yes


# data
genome=data/homotypic_intervals/intervals.fa
meme=data/homotypic_intervals/motif_more.meme

# output
homotypic_output=$res_dir/01_homotypic


# heterotypic
task=gene
gene_input_file=data/homotypic_intervals/intervals.txt
heterotypic_output=$res_dir/02_heterotypic
icthresh=4

mkdir -p $homotypic_output
mkdir -p $heterotypic_output


##################################### Running FIMO #####################################
print_green "Running fimo...\n"

$HOMOTYPIC               \
    -r $toolDir          \
    -o $homotypic_output \
    -k $maxk             \
    -n $topn             \
    -f $fimothresh       \
    -t $threads          \
    $genome              \
    $meme



################################### Heterotypic #####################################

# # remove genes not present in promoter_lengths.txt
# awk -F"\t" '{print $1"\t"}' homotypic_output/promoter_lengths.txt > homotypic_output/temp_genes_list.txt
# cat homotypic_output/temp_genes_list.txt | while read line; do
#     grep $line $gene_input_file
# done > genes/temp_${task}.txt
# rm homotypic_output/temp_genes_list.txt

print_green "Searching for heterotypic motif hits..."
mkdir -p $output
$HETEROTYPIC                                     \
    -d .                                         \
    -g $gene_input_file                          \
    -i $icthresh                                 \
    -p $homotypic_output/promoter_lengths.txt    \
    -b $homotypic_output/binomial_thresholds.txt \
    -c $homotypic_output/IC.txt                  \
    -f $homotypic_output/fimohits                \
    -o $heterotypic_output                       \
    -t $threads

cat $heterotypic_output/*.txt > $heterotypic_output/motif_output.txt
rm $heterotypic_output/temp*.txt


##################################### Heatmap ##################################

Rscript 05_heatmap.R                     \
    Overlap                              \
    $heterotypic_output/heatmap.png      \
    $heterotypic_output/motif_output.txt
