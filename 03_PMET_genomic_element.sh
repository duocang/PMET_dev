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

start_time=$SECONDS
################################ 1. input parameters ###################################
# tool
toolDir=scripts
HOMOTYPIC=$toolDir/PMETindex_genomic_element.sh
HETEROTYPIC=$toolDir/pmetParallel

chmod a+x $HOMOTYPIC
chmod a+x $HETEROTYPIC

threads=24
res_dir=results/03_PMET_genomic_element


# homotypic
overlap="NoOverlap"
utr="Yes"
topn=5000
maxk=5
length=1000
fimothresh=0.05
distance=1000
gff3id="gene_id="
delete_temp=yes

mrnaFull=No

# data
genome=data/homotypic_promoters/genome.fasta
anno=data/homotypic_promoters/anno.gff3
meme=data/Franco-Zorrilla_et_al_2014.meme

# genomic element
echo -e "Select the genomic element:\n    1. Three Prime UTR\n    2. Five Prime UTR\n    3. mRNA"
read -p "Enter your choice (1/2/3): " choice
case $choice in
    1) genomic_element="three_prime_UTR"; gff3id='Parent=transcript:' ;;
    2) genomic_element="five_prime_UTR"; gff3id='Parent=transcript:' ;;
    3) genomic_element="mRNA"; gff3id='ID=transcript:' ;;
    *) echo "Invalid choice. Please enter 1, 2, or 3."; exit 1 ;;
esac
print_fluorescent_yellow "You selected: $genomic_element"
print_fluorescent_yellow "GFF3 ID format: $gff3id"

# output
homotypic_output=$res_dir/01_homotypic_$genomic_element

# heterotypic
task=gene
gene_input_file=data/$task.txt
heterotypic_output=$res_dir/02_heterotypic_$genomic_element
icthresh=4

# plot output
plot_output=$res_dir/plot

mkdir -p $homotypic_output
mkdir -p $heterotypic_output
mkdir -p $plot_output


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

$HOMOTYPIC               \
    -r $toolDir          \
    -o $homotypic_output \
    -e $genomic_element  \
    -m $mrnaFull         \
    -i $gff3id           \
    -k $maxk             \
    -n $topn             \
    -p $length           \
    -v $overlap          \
    -u $utr              \
    -f $fimothresh       \
    -t $threads          \
    -d $delete_temp      \
    $genome              \
    $anno                \
    $meme


############################ 4. Running heterotypic ###############################
print_green "\n\nSearching for heterotypic motif hits...\n"

# remove genes not present in pre-computed pmet index
grep -Ff $homotypic_output/universe.txt $gene_input_file > $heterotypic_output/new_genes_temp.txt


$HETEROTYPIC                                     \
    -d .                                         \
    -g $heterotypic_output/new_genes_temp.txt    \
    -i $icthresh                                 \
    -p $homotypic_output/promoter_lengths.txt    \
    -b $homotypic_output/binomial_thresholds.txt \
    -c $homotypic_output/IC.txt                  \
    -f $homotypic_output/fimohits                \
    -o $heterotypic_output                       \
    -t $threads > $heterotypic_output/pmet.log

rm $heterotypic_output/new_genes_temp.txt
# merge pmet result
cat $heterotypic_output/*.txt > $heterotypic_output/motif_output.txt
rm $heterotypic_output/temp*.txt

# remove genes not present in pre-computed pmet index
grep -Ff $homotypic_output/universe.txt $gene_input_file > $heterotypic_output/new_genes_temp.txt

end_time=$SECONDS
elapsed_time=$((end_time - start_time))
days=$((elapsed_time/86400))
hours=$(( (elapsed_time%86400)/3600 ))
minutes=$(( (elapsed_time%3600)/60 ))
seconds=$((elapsed_time%60))
print_orange "      Time taken: $days day $hours hour $minutes minute $seconds second\n"


print_green "DONE: heterotypic search"

##################################### Heatmap ##################################
print_green "\n\nCreating heatmap...\n"

Rscript 05_heatmap.R                        \
    All                                  \
    $plot_output/heatmap.png             \
    $heterotypic_output/motif_output.txt \
    15                                    \
    3                                    \
    6                                    \
    FALSE

Rscript 05_heatmap.R                           \
    Overlap                                 \
    $plot_output/heatmap_overlap_unique.png \
    $heterotypic_output/motif_output.txt    \
    15                                       \
    3                                       \
    6                                       \
    TRUE

Rscript 05_heatmap.R                        \
    Overlap                              \
    $plot_output/heatmap_overlap.png     \
    $heterotypic_output/motif_output.txt \
    15                                    \
    3                                    \
    6                                    \
    FALSE

# method
# filename
# pmet.out
# topn
# histgram_ncol
# histgram_width
# unique_cmbination