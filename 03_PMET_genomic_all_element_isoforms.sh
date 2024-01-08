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
print_middle(){
    FLUORESCENT_YELLOW='\033[1;33m'
    NC='\033[0m' # No Color
    # 获取终端的宽度
    COLUMNS=$(tput cols)
    # 遍历每一行
    while IFS= read -r line; do
        # 计算需要的空格数来居中文本
        padding=$(( (COLUMNS - ${#line}) / 2 ))
        printf "%${padding}s" ''
        printf "${FLUORESCENT_YELLOW}${line}${NC}\n"
    done <<< "$1"
}


# Give execute permission to all users for the file.
find . -type f \( -name "*.sh" -o -name "*.pl" \) -exec chmod a+x {} \;


echo -e "\n\n"
print_middle "The purpose of this script is to                                      \n"
print_middle "  1. List all isoforms of genes in each cluster you provide             "
print_middle "  2. Search motif pairs on each isoform of any genomic elemet           "
print_middle "                                                                      \n"



########################## 1. Downloading data #######################################
cd data
if [ -f "TAIR10.gff3" ]; then
    print_green "Genome and annotation are ready!"
else
    print_fluorescent_yellow "Downloading genome and annotation...\n"
    bash ./fetch_data.sh
fi
cd ..

start_time=$SECONDS
################################ 2. input parameters ###################################
# tool
toolDir=scripts
HOMOTYPIC=$toolDir/PMETindex_genomic_element_all_isoforms.sh
HETEROTYPIC=$toolDir/pmetParallel

chmod a+x $HOMOTYPIC
chmod a+x $HETEROTYPIC

threads=8
res_dir=results/03_genomic_element_isoforms

# homotypic
overlap="NoOverlap"
utr="Yes"
topn=5000
maxk=5
length=1000
fimothresh=0.05
distance=1000
gff3id="gene_id="
delete_temp=no

# data
genome=data/TAIR10.fasta
anno=data/TAIR10.gff3
meme=data/Franco-Zorrilla_et_al_2014.meme


# genomic element
echo -e "Select the genomic element:\n    1. 3' UTR\n    2. 5' UTR\n    3. mRNA\n    4. CDS\n    5. Exon"
read -p "Enter your choice (1/2/3/4): " choice
case $choice in
    1) genomic_element="three_prime_UTR"; gff3id='Parent=transcript:' ;;
    2) genomic_element="five_prime_UTR"; gff3id='Parent=transcript:' ;;
    3) genomic_element="mRNA"; gff3id='ID=transcript:';;
    4) genomic_element="CDS"; gff3id='Parent=transcript:' ;;
    5) genomic_element="exon"; gff3id='Parent=transcript:' ;;
    *) echo "Invalid choice. Please enter 1, 2, or 3."; exit 1 ;;
esac
print_fluorescent_yellow "Chosen Genomic Element: $genomic_element"
print_fluorescent_yellow "GFF3 ID format: $gff3id"

# output
homotypic_output=$res_dir/01_homotypic

# heterotypic
task=salt_top300
gene_input_file=data/genes/$task.txt
heterotypic_output=$res_dir/02_heterotypic
icthresh=4

# plot output
plot_output=$res_dir/03_plot

mkdir -p $homotypic_output
# mkdir -p $heterotypic_output
# mkdir -p $plot_output

# ############################## 3. Running homotypic #################################
# print_green "Running homotypic searching...\n"

# $HOMOTYPIC               \
#     -r $toolDir          \
#     -o $homotypic_output \
#     -e $genomic_element  \
#     -i $gff3id           \
#     -k $maxk             \
#     -n $topn             \
#     -p $length           \
#     -v $overlap          \
#     -u $utr              \
#     -f $fimothresh       \
#     -t $threads          \
#     -d $delete_temp      \
#     $genome              \
#     $anno                \
#     $meme


############################ 4. Running heterotypic ###############################
print_green "\nSearching for heterotypic motif hits..."

# 创建一个临时的 R 脚本 Creating a Temporary R Script
# R 脚本的目的是列出所有的基因isoform The purpose of the R script is to extend all gene isoform
temp_r_script=$(mktemp)
# 将 R 代码写入临时文件 Write R code to a temporary file
cat <<EOF >"$temp_r_script"
suppressPackageStartupMessages(library(dplyr))

args            <- commandArgs(trailingOnly = TRUE)
input_file      <- args[1]
output_file     <- args[2]
max_num_isoform <- as.integer(args[3])
universe        <- args[4]

full_genes    <- read.table(universe)\$V1
genes_cluster <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)

# 扩展基因列表 Expand the gene isoforms
expanded_genes <- lapply(genes_cluster\$V2, function(gene) {
    paste0(gene, ".", seq_len(max_num_isoform))
}) %>% unlist()

expanded_clustter <- lapply(genes_cluster\$V1, function(x){rep(x, max_num_isoform)}) %>% unlist()

df <- data.frame(expanded_clustter, expanded_genes) %>%
  filter(expanded_genes %in% full_genes)

write.table(df, file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
EOF

# 从universe.txt中获取所有基因中isoform的最大数值
# Get the maximum value of isoform in all genes from universe.txt
max_num_isoform=0
while read -r line; do
    num=${line##*.}
    if (( num > max_num_isoform )); then
        max_num_isoform=$num
    fi
done < $homotypic_output/universe.txt

# Perform pmet analysis on different types of gene clusters
for task in "salt_top300" "random_genes_300" "genes_cell_type_treatment" "gene_cortex_epidermis_pericycle" "heat_top300"; do
    print_orange "Gens: $task"

    heterotypic_output=${heterotypic_output}_${task}
    plot_output=${plot_output}_${task}
    gene_input_file=data/genes/$task.txt
    mkdir -p $heterotypic_output
    mkdir -p $plot_output
    # create gene with all possible isofrms in clusters
    Rscript "$temp_r_script"                     \
        "$gene_input_file"                       \
        "$heterotypic_output/new_genes_temp.txt" \
         $max_num_isoform                        \
         $homotypic_output/universe.txt

    ##################################### PMET ##################################
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

    end_time=$SECONDS
    elapsed_time=$((end_time - start_time))
    days=$((elapsed_time/86400))
    hours=$(( (elapsed_time%86400)/3600 ))
    minutes=$(( (elapsed_time%3600)/60 ))
    seconds=$((elapsed_time%60))
    print_orange "      Time taken: $days day $hours hour $minutes minute $seconds second\n"
    # print_green "DONE: heterotypic search"

    ##################################### Heatmap ##################################
    print_green "Creating heatmap..."

    Rscript 05_heatmap.R                     \
        All                                  \
        $plot_output/heatmap.png             \
        $heterotypic_output/motif_output.txt \
        15                                   \
        3                                    \
        6                                    \
        FALSE > $plot_output/motifs.txt

    Rscript 05_heatmap.R                        \
        Overlap                                 \
        $plot_output/heatmap_overlap_unique.png \
        $heterotypic_output/motif_output.txt    \
        15                                      \
        3                                       \
        6                                       \
        TRUE > $plot_output/log.txt

    Rscript 05_heatmap.R                     \
        Overlap                              \
        $plot_output/heatmap_overlap.png     \
        $heterotypic_output/motif_output.txt \
        15                                   \
        3                                    \
        6                                    \
        FALSE > $plot_output/log.txt
    rm -rf $plot_output/log.txt
done
# method
# filename
# pmet.out
# topn
# histgram_ncol
# histgram_width
# unique_cmbination

rm "$temp_r_script"
