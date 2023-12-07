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



################################ 1. input parameters ###################################
# tool
toolDir=scripts
HOMOTYPIC=$toolDir/PMETindex_promoters_with_distance_tss_fimo_intergrated.sh
HETEROTYPIC=$toolDir/pmetParallel

chmod a+x $HOMOTYPIC
chmod a+x $HETEROTYPIC

threads=24
res_dir=results/07_homotypic_promoters_with_distance_to_tss

# homotypic
gff3id="gene_id="
overlap="NoOverlap"
utr="No"
topn=5000
maxk=5
length=300
fimothresh=0.05
gap=100
promlengthlimit=300
gff3id="gene_id="
delete_temp=yes

# data
genome=data/homotypic_promoters/genome.fasta
anno=data/homotypic_promoters/anno.gff3
meme=data/Franco-Zorrilla_et_al_2014.meme

# output
homotypic_output=$res_dir/01_homotypic


# heterotypic
task=gene
gene_input_file=data/$task.txt
heterotypic_output=$res_dir/02_heterotypic
icthresh=4

mkdir -p $homotypic_output
mkdir -p $heterotypic_output

################################ 2. Downloading data #######################################
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

################################ Running homotypic ###################################
print_green "Running fimo...\n"

$HOMOTYPIC                \
    -r $toolDir           \
    -o $homotypic_output  \
    -i $gff3id            \
    -k $maxk              \
    -n $topn              \
    -p $length            \
    -g $gap               \
    -l $promlengthlimit   \
    -v $overlap           \
    -u $utr               \
    -f $fimothresh        \
    -t $threads           \
    -d $delete_temp       \
    $genome               \
    $anno                 \
    $meme


# ########################## Running heterotypic ##################################
print_green "\n\nSearching for heterotypic motif hits..."

# remove genes not present in pre-computed pmet index
grep -Ff $homotypic_output/universe.txt $gene_input_file > $gene_input_file"temp"

$HETEROTYPIC                                     \
    -d .                                         \
    -g $gene_input_file"temp"                    \
    -i $icthresh                                 \
    -p $homotypic_output/promoter_lengths.txt    \
    -b $homotypic_output/binomial_thresholds.txt \
    -c $homotypic_output/IC.txt                  \
    -f $homotypic_output/fimohits                \
    -o $heterotypic_output                       \
    -t $threads

cat $heterotypic_output/*.txt > $heterotypic_output/motif_output.txt
rm $heterotypic_output/temp*.txt
rm $gene_input_file"temp"


##################################### Heatmap ##################################

Rscript 05_heatmap.R                     \
    Overlap                              \
    $heterotypic_output/heatmap.png      \
    $heterotypic_output/motif_output.txt
