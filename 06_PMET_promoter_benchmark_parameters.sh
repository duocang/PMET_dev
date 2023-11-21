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
# Homotypic
genome=data/homotypic_promoters/genome.fasta
anno=data/homotypic_promoters/anno.gff3
meme=data/homotypic_promoters/motif.meme


topnRangeStr="2500-5000"
maxkRangeStr="5-6"
promlengthRangeStr="1000-1500"

indexOutput=results/06_PMET_promoter_benchmark_parameters/homotypic
mkdir -p $indexOutput

# Heterotypic
task=gene_cell_identities
gene_file=data/$task.txt
fimothresh=005

output=results/06_PMET_promoter_benchmark_parameters/heterotypic
logDir=results/06_PMET_promoter_benchmark_parameters/logs
mkdir -p $output
mkdir -p $logDir

ncpu=16
################################## Homotypic #######################################
chmod a+x scripts/PMETindex_promoters_benchmark_parameters.sh
scripts/PMETindex_promoters_benchmark_parameters.sh \
    -r scripts                                      \
    -o $indexOutput                                 \
    -i gene_id=                                     \
    -v NoOverlap                                    \
    -u Yes                                          \
    -t $ncpu                                        \
    -n $topnRangeStr                                \
    -k $maxkRangeStr                                \
    -p $promlengthRangeStr                          \
    $genome                                         \
    $anno                                           \
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
            print_fluorescent_yellow "        genes     : $gene_file"
            print_fluorescent_yellow "        promlength: $promlength"
            print_fluorescent_yellow "        maxk      : $maxk"
            print_fluorescent_yellow "        topn      : $topn"
            print_fluorescent_yellow "        fimothresh: $fimothresh"
            print_green              "    Searching for heterotypic motifs pairs with ${ncpu} CPUS..."

            indexOutputTemp=$indexOutput/fimoThresh${fimothresh}_promLength${promlength}_maxK${maxk}_topN${topn}
            outputTemp=$output/${task}_fimoThresh${fimothresh}_promLength${promlength}_maxK${maxk}_topN${topn}

            mkdir -p $outputTemp

            # remove genes not present in pre-computed pmet index
            grep -Ff $indexOutputTemp/universe.txt $gene_file > $outputTemp/available_genes.txt

            scripts/pmetParallel                             \
                -d .                                         \
                -g $outputTemp/available_genes.txt           \
                -i 4                                         \
                -p $indexOutputTemp/promoter_lengths.txt     \
                -b $indexOutputTemp/binomial_thresholds.txt  \
                -c $indexOutputTemp/IC.txt                   \
                -f $indexOutputTemp/fimohits                 \
                -o $outputTemp                               \
                -t $ncpu > $logDir/${task}_fimoThresh${fimothresh}_promLength${promlength}_maxK${maxk}_topN${topn}.log 2>&1

            rm  $outputTemp/available_genes.txt
            cat $outputTemp/*.txt > $outputTemp/motif_output.txt
            rm  $outputTemp/temp*.txt

            grep -Ff $indexOutputTemp/universe.txt $gene_file > $outputTemp/available_genes.txt
        done
    done
done
