#!/bin/bash

# cd src/meme-5.5.3

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

################################ 1. input parameters ####################################
chmod a+x scripts/PMETindex_promoters_benchmark_parameters.sh

# homotypic
ncpu=16
indexOutput=results/06_benchmark_parameters/homotypic_output

genome=data/homotypic_promoters/genome.fasta
anno=data/homotypic_promoters/anno.gff3
meme=data/Franco-Zorrilla_et_al_2014.meme

# heterotypic
fimothresh=005
task=gene_cell_identities
gene_file=data/$task.txt

output=results/06_benchmark_parameters/heterotypic_output

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

start=$(date +%s)

############################## 3. Running homotypic #################################
mkdir -p $indexOutput

scripts/PMETindex_promoters_benchmark_parameters.sh \
    -r scripts                                      \
    -o $indexOutput                                 \
    -i gene_id=                                     \
    -v NoOverlap                                    \
    -u Yes                                          \
    -t $ncpu                                        \
    -n "5000"                                       \
    -k "1-2-3-4-5-6-7-8-9"                          \
    -p "50-200-500-1000-1500-2500-5000"             \
    $genome                                         \
    $anno                                           \
    $meme


############################ 4. Running heterotypic ###############################

mkdir -p $output/logs

for promlength in 50 200 500 1000 1500 2500 5000; do
    for maxk in 1 2 3 4 5 6 7 8 9; do
        for topn in 5000; do

            print_orange             "    parameters:"
            print_fluorescent_yellow "        genes     : $gene_file"
            print_fluorescent_yellow "        promlength: $promlength"
            print_fluorescent_yellow "        maxk      : $maxk"
            print_fluorescent_yellow "        topn      : $topn"
            print_fluorescent_yellow "        fimothresh: $fimothresh"
            print_red                "    Searching for heterotypic motifs pairs with ${ncpu} CPUS..."

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
                -t $ncpu > $output/logs/${task}_fimoThresh${fimothresh}_promLength${promlength}_maxK${maxk}_topN${topn}.log 2>&1

            rm  $outputTemp/available_genes.txt
            cat $outputTemp/*.txt > $outputTemp/motif_output.txt
            rm  $outputTemp/temp*.txt

            grep -Ff $indexOutputTemp/universe.txt $gene_file > $outputTemp/available_genes.txt
        done
    done
done

print_green("DONE")
