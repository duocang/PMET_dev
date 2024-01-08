#!/bin/bash
set -e

function usage () {
    cat >&2 <<EOF
USAGE: PMETindexgenome [options] <genome> <gff3> <memefile>

Creates PMET index for Paired Motif Enrichment Test using genome files.
Required arguments:
-r <PMETindex_path>	: Full path of python scripts called from this file. Required.
-i <gff3_identifier> : gene identifier in gff3 file e.g. gene_id=

Optional arguments:
-o <output_directory> : Output directory for results
-n <topn>	: How many top promoter hits to take per motif. Default=5000
-k <max_k>	: Maximum motif hits allowed within each promoter.  Default: 5
-p <promoter_length>	: Length of promoter in bp used to detect motif hits default: 1000
-v <include_overlaps> :  Remove promoter overlaps with gene sequences. AllowOverlap or NoOverlap, Default : AllowOverlap
-u <include_UTR> : Include 5' UTR sequence? Yes or No, default : No
-f <fimo_threshold> : Specify a minimum quality for hits matched by fimo. Default: 0.05

EOF
}

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


# set up defaults
topn=5000
maxk=5
promlength=1000
fimothresh=0.05
element=mRNA
mrnaFull=No
delete=yes
overlap="AllowOverlap"
utr="No"
gff3id='transcript:'
pmetroot="scripts"
threads=4
icthreshold=24

# set up empty variables

indexingOutputDir=
genomefile=
gff3file=
memefile=

# deal with arguments
# if none, exit
if [ $# -eq 0 ]
    then
        echo "No arguments supplied"  >&2
        usage
        exit 1
fi

while getopts ":r:i:o:n:k:p:f:v:u:e:m:t:d:" options; do
    case $options in
        r) print_white "Full path of PMET_index                : "; print_orange "$OPTARG" >&2
        pmetroot=$OPTARG;;
        i) print_white "GFF3 feature identifier                : "; print_orange "$OPTARG" >&2
        gff3id=$OPTARG;;
        o) print_white "Output directory for results           : "; print_orange "$OPTARG" >&2
        indexingOutputDir=$OPTARG;;
        n) print_white "Top n promoter hits to take per motif  : "; print_orange "$OPTARG" >&2
        topn=$OPTARG;;
        k) print_white "Top k motif hits within each promoter  : "; print_orange "$OPTARG" >&2
        maxk=$OPTARG;;
        p) print_white "Promoter length                        : "; print_orange "$OPTARG" >&2
        promlength=$OPTARG;;
        f) print_white "Fimo threshold                         : "; print_orange "$OPTARG" >&2
        fimothresh=$OPTARG;;
        v) print_white "Remove promoter overlaps with sequences: "; print_orange "$OPTARG" >&2
        overlap=$OPTARG;;
        u) print_white "Include 5' UTR sequence?               : "; print_orange "$OPTARG" >&2
        utr=$OPTARG;;
        e) print_white "Genomic element                        : "; print_orange "$OPTARG" >&2
        element=$OPTARG;;
        m) print_white "Include complete mRNA (5' and 3' UTR)? : "; print_orange "$OPTARG" >&2
        mrnaFull=$OPTARG;;
        t) print_white "Number of threads                      : "; print_orange "$OPTARG" >&2
        threads=$OPTARG;;
        d) print_white "Delete unnecssary files                : "; print_orange "$OPTARG" >&2
        delete=$OPTARG;;
        \?) print_red  "Invalid option: -$OPTARG" >&2
        exit 1;;
        :)  print_red "Option -$OPTARG requires an argument." >&2
        exit 1;;
    esac
done


shift $((OPTIND - 1))
genomefile=$1
gff3file=$2
memefile=$3
universefile=$indexingOutputDir/universe.txt
bedfile=$indexingOutputDir/${element}.bed

print_white "Genome file                            : "; print_orange $genomefile
print_white "Annotation file                        : "; print_orange $gff3file
print_white "Motif meme file                        : "; print_orange $memefile

mkdir -p $indexingOutputDir

start=$SECONDS

print_green "Preparing data for FIMO and PMET index..."

# -------------------------------------------------------------------------------------------
# 1. sort annotaion by gene coordinates
print_fluorescent_yellow "     1. Sorting annotation by gene coordinates"
$pmetroot/gff3sort/gff3sort.pl $gff3file > $indexingOutputDir/sorted.gff3

# -------------------------------------------------------------------------------------------
# 2. extract gene line from annoitation
print_fluorescent_yellow "     2. Extracting gene line from annoitation"
# grep -P '\mRNA\t' $indexingOutputDir/sorted.gff3 > $indexingOutputDir/genelines.gff3
if [[ "$(uname)" == "Linux" ]]; then
    grep -P "\t${element}\t" $indexingOutputDir/sorted.gff3 > $indexingOutputDir/genelines.gff3
elif [[ "$(uname)" == "Darwin" ]]; then
    grep    "\t${element}\t" $indexingOutputDir/sorted.gff3 > $indexingOutputDir/genelines.gff3
else
    print_red "Unsupported operating system."
fi

# -------------------------------------------------------------------------------------------
# 3. extract chromosome , start, end, gene ('gene_id' for input) ...
print_fluorescent_yellow "     3. Extracting chromosome, start, end, gene..."
python3 $pmetroot/parse_mRNAlines.py $gff3id $indexingOutputDir/genelines.gff3 $bedfile

# -------------------------------------------------------------------------------------------
# 4. filter invalid genes: start should be smaller than end
invalidRows=$(awk '$2 >= $3' $bedfile)
if [[ -n "$invalidRows" ]]; then
    echo "$invalidRows" > $indexingOutputDir/invalid_lines.bed
fi
print_fluorescent_yellow "     4. Extracting genes coordinates: start should be smaller than end (genelines.bed)"
# 在BED文件格式中，无论是正链（+）还是负链（-），起始位置总是小于终止位置。
# In the BED file format, the start position is always less than the end position for both positive (+) and negative (-) chains.
# 起始和终止位置是指定基因上的物理位置，而不是表达或翻译的方向。
# start and end positions specify the physical location of the gene, rather than the direction of expression or translation.
# starting site < stopped site in bed file
awk '$2 <  $3' $bedfile > temp.bed && mv temp.bed $bedfile


# # -------------------------------------------------------------------------------------------
# # 5. remove rows with same chromosome and coordinates
print_fluorescent_yellow "     5. Removing rows with same chromosome and coordinates"
awk -v output_dir="$indexingOutputDir" '{
    key = $1 "\t" $2 "\t" $3
    if (seen[key]++) {
        print $0 > (output_dir "/duplicated.bed")
    } else {
        print $0
    }
}' "$bedfile" > "$indexingOutputDir/unique.bed"

rm -rf $bedfile
mv "$indexingOutputDir/unique.bed" $bedfile



# -------------------------------------------------------------------------------------------
# 6. list of all genes found
print_fluorescent_yellow "     6. Extracting all genes found (universe.txt)"
cut -f 4 $bedfile > $universefile
awk '{
    # 使用 split 函数去除点和后面的数字，只保留基因名
    split($1, arr, ".");
    gene = arr[1];

    # 如果基因名还未出现过，则打印并标记为已出现
    if (!seen[gene]++) {
        print gene;
    }
}' "$universefile" > $indexingOutputDir/cleaned_universe.txt
rm -rf "$universefile"
mv $indexingOutputDir/cleaned_universe.txt "$universefile"

# -------------------------------------------------------------------------------------------
# 7. promoter lenfths from promoters.bed
print_fluorescent_yellow "     7. Promoter lengths from genelines.bed"
awk '{print $4 "\t" ($3 - $2)}' $bedfile      \
    > $indexingOutputDir/promoter_lengths.txt

# -------------------------------------------------------------------------------------------
# 8. strip the potential FASTA line breaks. creates genome_stripped.fa
print_fluorescent_yellow "     8. Removing potential FASTA line breaks (genome_stripped.fa)"
awk '/^>/ { if (NR!=1) print ""; printf "%s\n",$0; next;} \
    { printf "%s",$0;} \
    END { print ""; }'  $genomefile > $indexingOutputDir/genome_stripped.fa


# -------------------------------------------------------------------------------------------
# 9. create promoters fasta
print_fluorescent_yellow "     9. Creating mRNA FASTA file (promoter_rought.fa)";
bedtools getfasta \
        -fi  $indexingOutputDir/genome_stripped.fa \
        -bed $bedfile                              \
        -fo  $indexingOutputDir/promoter_rought.fa \
        -name -s

# -------------------------------------------------------------------------------------------
# 11. replace the id of each seq with gene names
print_fluorescent_yellow "    11. Replacing the id of each sequences' with gene names (promoter.fa)"
sed 's/::.*//g' $indexingOutputDir/promoter_rought.fa > $indexingOutputDir/promoter.fa

# check if any duplicated id
grep "^>" $indexingOutputDir/promoter.fa > $indexingOutputDir/ids.txt
sort $indexingOutputDir/ids.txt | uniq -d > $indexingOutputDir/duplicate_ids.txt
if [ ! -s $indexingOutputDir/duplicate_ids.txt ]; then
    rm -rf $indexingOutputDir/ids.txt
    rm -rf $indexingOutputDir/duplicate_ids.txt
fi

# -------------------------------------------------------------------------------------------
# 11. promoter.bg from promoter.fa
print_fluorescent_yellow "    11. fasta-get-markov estimates a Markov model from promoter.fa. (promoter.bg)"
fasta-get-markov $indexingOutputDir/promoter.fa > $indexingOutputDir/promoter.bg

# -------------------------------------------------------------------------------------------
# 12. IC.txt
print_fluorescent_yellow "    12. Generating information content (IC.txt)"
[ ! -d $indexingOutputDir/memefiles ] && mkdir $indexingOutputDir/memefiles
python3 $pmetroot/parse_memefile.py $memefile $indexingOutputDir/memefiles/
python3 $pmetroot/calculateICfrommeme_IC_to_csv.py \
    $indexingOutputDir/memefiles/                  \
    $indexingOutputDir/IC.txt
rm -rf $indexingOutputDir/memefiles/*

# -------------------------------------------------------------------------------------------
# 13. individual motif files from user's meme file
print_fluorescent_yellow "    13. Spliting motifs into individual meme files (folder memefiles)"
python3 $pmetroot/parse_memefile_batches.py $memefile $indexingOutputDir/memefiles/ $threads

# -------------------------------- Run fimo and pmetindex --------------------------
[ ! -d $indexingOutputDir/fimohits ] && mkdir $indexingOutputDir/fimohits

print_green "Running FIMO and PMET index..."
runFimoIndexing () {
    memefile=$1
    indexingOutputDir=$2
    fimothresh=$3
    pmetroot=$4
    maxk=$5
    topn=$6
    filename=`basename $memefile .txt`

    $pmetroot/fimo                              \
        --no-qvalue                             \
        --text                                  \
        --thresh $fimothresh                    \
        --verbosity 1                           \
        --bgfile $indexingOutputDir/promoter.bg \
        --topn $topn                            \
        --topk $maxk                            \
        --oc $indexingOutputDir/fimohits        \
        $memefile                               \
        $indexingOutputDir/promoter.fa          \
        $indexingOutputDir/promoter_lengths.txt
}
export -f runFimoIndexing


nummotifs=$(grep -c '^MOTIF' "$memefile")
print_orange "    $nummotifs motifs found"

find $indexingOutputDir/memefiles -name \*.txt \
    | parallel --progress --jobs=$threads \
        "runFimoIndexing {} $indexingOutputDir $fimothresh $pmetroot $maxk $topn"

mv $indexingOutputDir/fimohits/binomial_thresholds.txt $indexingOutputDir/


print_fluorescent_yellow "    14. Adding up lengths"
# add up lengths and keep only one genes for multiple mRAN fragments
awk '{
    # 分割基因名和亚型，只保留基因名 Split gene names and isoforms, keep only gene names
    split($1, arr, ".");
    gene = arr[1];

    # 累加长度 Cumulative length
    total_length[gene] += $2;
}
END {
    # 打印最终结果 Print the final result
    for (gene in total_length) {
        print gene "\t" total_length[gene];
    }
}' $indexingOutputDir/promoter_lengths.txt > $indexingOutputDir/summed_lengths.txt

rm -rf $indexingOutputDir/promoter_lengths.txt
mv $indexingOutputDir/summed_lengths.txt $indexingOutputDir/promoter_lengths.txt


print_fluorescent_yellow "    15. Merging and filtering fimo hits out of multiple mRNA fragments"
# 创建一个临时的 R 脚本 Creating a Temporary R Script
# R 脚本的目的是列出所有的基因isoform The purpose of the R script is to extend all gene isoform
temp_r_script=$(mktemp)
# 将 R 代码写入临时文件 Write R code to a temporary file
cat <<EOF >"$temp_r_script"
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: script_name.R binomial_thresholds_file fimohits_dir output_dir\n")
  quit(status = 1)
}

binomial_thresholds_file <- args[1]
fimohits_dir             <- args[2]
output_dir               <- args[3]

dir.create(output_dir, showWarnings = FALSE)

binomial_thresholds <- fread(binomial_thresholds_file) %>% as.data.frame()

for (i in 1:nrow(binomial_thresholds)) {
  motif <- binomial_thresholds[i, "V1"]
  binomial_threshold <- binomial_thresholds[i, "V2"]
  a <- file.path(fimohits_dir, paste0(motif, ".txt")) %>%
    fread() %>% as.data.frame()
  # str_remove("AT4G25700.2", "\\.\\d+")
  a\$V2 <- str_remove(a\$V2, "\\\\.\\\\d+")
  a <- a %>%
    arrange(V2, V7) %>%
    group_by(V2) %>%
    arrange(V7) %>%
    slice_head(n = 5) %>%
    filter(V7 < binomial_threshold)

  write.table(a, file.path(output_dir, paste0(motif, ".txt")),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t")
}
EOF

Rscript $temp_r_script \
    $indexingOutputDir/binomial_thresholds.txt    \
    $indexingOutputDir/fimohits                   \
    $indexingOutputDir/fimohits_
rm -rf $indexingOutputDir/fimohits
mv $indexingOutputDir/fimohits_ $indexingOutputDir/fimohits
rm "$temp_r_script"

# -------------------------------------------------------------------------------------------
# Deleting unnecessary files
if [[ $delete == "yes" || $delete == "YES" || $delete == "Y" || $delete == "y" ]]; then
    print_green "Deleting unnecessary files...\n\n"
    # rm -rf $indexingOutputDir/IC.txt
    # rm -rf $indexingOutputDir/binomial_thresholds.txt
    rm -rf $indexingOutputDir/filtered_bedfile.bed
    # rm -rf $indexingOutputDir/fimohits
    rm -rf $indexingOutputDir/fimohits_
    rm -rf $indexingOutputDir/five_prime_UTR.bed
    rm -rf $indexingOutputDir/genelines.gff3
    rm -rf $indexingOutputDir/genome_stripped.fa
    rm -rf $indexingOutputDir/genome_stripped.fa.fai
    rm -rf $indexingOutputDir/id_duplicated.txt
    rm -rf $indexingOutputDir/mRNA.bed
    # rm -rf $indexingOutputDir/memefiles
    rm -rf $indexingOutputDir/modified_bedfile.bed
    rm -rf $indexingOutputDir/promoter.bg
    rm -rf $indexingOutputDir/promoter.fa
    # rm -rf $indexingOutputDir/promoter_lengths.txt
    rm -rf $indexingOutputDir/promoter_rought.fa
    rm -rf $indexingOutputDir/sorted.gff3
    rm -rf $indexingOutputDir/three_prime_UTR.bed
    # rm -rf $indexingOutputDir/universe.txt
    rm -rf $indexingOutputDir/with_overlapping.bed
fi

# -------------------------------------------------------------------------------------------
# Checking results
print_green "Checking results...\n\n"
# 计算 $indexingOutputDir/fimohits 目录下 .txt 文件的数量
# Count the number of .txt files in the $indexingOutputDir/fimohits directory
file_count=$(find "$indexingOutputDir/fimohits" -maxdepth 1 -type f -name "*.txt" | wc -l)

# 检查文件数量是否等于 meotif的数量 （$nummotifs）
# Check if the number of files equals the number of meotifs ($nummotifs)
if [ "$file_count" -eq "$nummotifs" ]; then
    end=$SECONDS
    elapsed_time=$((end - start))
    days=$((elapsed_time/86400))
    hours=$(( (elapsed_time%86400)/3600 ))
    minutes=$(( (elapsed_time%3600)/60 ))
    seconds=$((elapsed_time%60))
    print_orange "      Time take: $days day $hours hour $minutes minute $seconds seconds"

    print_green "DONE: homotypic search"
else
    print_red "\nError: there are $file_count fimohits files, it should be $nummotifs."
fi

# # next stage needs the following inputs

# #   promoter_lengths.txt        made by parse_promoter_lengths.py from .bed file
# #   bimnomial_thresholds.txt    made by PMETindex
# #   IC.txt                      made by calculateICfrommeme.py from meme file
# #   gene input file             supplied by user
