#!/bin/bash

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

########################## Downloading data #######################################
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





################################# Parameters #######################################
# tool
toolDir=scripts
HOMOTYPIC=$toolDir/PMETindex_promoters_benchmark_parameters.sh
HETEROTYPIC=$toolDir/pmetParallel

chmod a+x $HOMOTYPIC
chmod a+x $HETEROTYPIC

threads=24
res_dir=results/06_PMET_promoter_benchmark_parameters

# Homotypic
gff3id="gene_id="
noOverlap="NoOverlap"
utr="Yes"
topnRangeStr="2500-5000"
maxkRangeStr="5-6"
promlengthRangeStr="1000-1500"

# data
genome=data/homotypic_promoters/genome.fasta
anno=data/homotypic_promoters/anno.gff3
meme=data/homotypic_promoters/motif.meme

# output
homotypic_output=$res_dir/01_homotypic


# Heterotypic
task=gene_cell_identities
gene_input_file=data/$task.txt
fimothresh=005
icthresh=4

heterotypic_output=$res_dir/02_heterotypic
logDir=$res_dir/logs

mkdir -p $homotypic_output
mkdir -p $heterotypic_output
mkdir -p $logDir


################################## Homotypic #######################################
chmod a+x $HOMOTYPIC
$HOMOTYPIC                 \
    -r $toolDir            \
    -o $homotypic_output   \
    -i $gff3id             \
    -v $noOverlap          \
    -u $utr                \
    -t $threads            \
    -n $topnRangeStr       \
    -k $maxkRangeStr       \
    -p $promlengthRangeStr \
    $genome                \
    $anno                  \
    $meme




################################## Heterotypic #######################################
# 定义函数以将逗号分隔的字符串转换为数组
string_to_array() {
    local IFS='-' # 设置字段分隔符为逗号
    read -ra arr <<< "$1" # 读取字符串并分割到数组
    echo "${arr[@]}" # 返回数组
}
topnRange=($(string_to_array "$topnRangeStr"))
maxkRange=($(string_to_array "$maxkRangeStr"))
promlengthRange=($(string_to_array "$promlengthRangeStr"))

for promlength in ${promlengthRange[@]}; do
    for maxk in ${maxkRange[@]}; do
        for topn in ${topnRange[@]}; do

            print_orange             "    parameters:"
            print_fluorescent_yellow "        genes     : $gene_input_file"
            print_fluorescent_yellow "        promlength: $promlength"
            print_fluorescent_yellow "        maxk      : $maxk"
            print_fluorescent_yellow "        topn      : $topn"
            print_fluorescent_yellow "        fimothresh: $fimothresh"
            print_green              "    Searching for heterotypic motifs pairs with ${threads} CPUS..."

            indexOutput=$homotypic_output/fimoThresh${fimothresh}_promLength${promlength}_maxK${maxk}_topN${topn}
            output=$heterotypic_output/${task}_fimoThresh${fimothresh}_promLength${promlength}_maxK${maxk}_topN${topn}

            mkdir -p $output

            # remove genes not present in pre-computed pmet index
            grep -Ff $indexOutput/universe.txt $gene_input_file > $output/available_genes.txt

            $HETEROTYPIC                                 \
                -d .                                     \
                -g $output/available_genes.txt           \
                -i $icthresh                             \
                -p $indexOutput/promoter_lengths.txt     \
                -b $indexOutput/binomial_thresholds.txt  \
                -c $indexOutput/IC.txt                   \
                -f $indexOutput/fimohits                 \
                -o $output                               \
                -t $threads > $logDir/${task}_fimoThresh${fimothresh}_promLength${promlength}_maxK${maxk}_topN${topn}.log 2>&1

            rm  $output/available_genes.txt
            cat $output/*.txt > $output/motif_output.txt
            rm  $output/temp*.txt

            Rscript 05_heatmap.R         \
                Overlap                  \
                $output/heatmap.png      \
                $output/motif_output.txt

            grep -Ff $indexOutput/universe.txt $gene_input_file > $output/available_genes.txt
        done
    done
done
