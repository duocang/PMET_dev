#!/bin/bash
set -e

# called when user selects Promoters on web UI
# Takes as inputs from web page : fasta file , gff3 file, meme file, gene clusters file
# gfff3 identifier, N, k, promoter length, overlap included?, utr included?, [fimo threshold]

#This version requires outputdir, the 3 input files will be put there


# ICdict.pickle file, binomial_thresholds.txt file
# a directory called fimohits contaning files called motifname.txt fo reach mtoif in meme file input

# This is all the input needed by pmet tool, along with genes list file (clusters) that user will provide via web UI



# 22.1.18 Charlotte Rich
# last edit: 7.2.18 - removed the make 1 big fimohits files
#Edited by PEB June 2020. Will be called by run_pmet_min.php to create inputs for call to pmet binary
#Take progress file from 1% to 95%

#PEB DEc 2020. Called by run_pmet.php with params
# -r ./scripts -o jobdir/indexoutput/ -i gene_id= -k k -n N -p promlength -u Yes|No -v NoOverlap|AllowOverlap  fastafile gtffile memefile

# fimo threshold is default

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

while getopts ":r:i:o:n:k:p:f:g:v:u:e:t:d:" options; do
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
bedfile=$indexingOutputDir/genelines.bed

print_white "Genome file                            : "; print_orange $genomefile
print_white "Annotation file                        : "; print_orange $gff3file
print_white "Motif meme file                        : "; print_orange $memefile

mkdir -p $indexingOutputDir

start=$SECONDS

print_green "Preparing data for FIMO and PMET index..."


# -------------------------------------------------------------------------------------------
# 1. sort annotaion by gene coordinates
print_fluorescent_yellow "     1. Sorting annotation by gene coordinates"
chmod a+x $pmetroot/gff3sort/gff3sort.pl
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
    echo "$invalidRows" > $indexingOutputDir/invalid_mRNAlines.bed
fi
# awk '$2 >= $3' $bedfile > $indexingOutputDir/invalid_mRNAlines.bed

print_fluorescent_yellow "     4. Extracting genes coordinates: start should be smaller than end (genelines.bed)"
awk '$2 <  $3' $bedfile > temp.bed && mv temp.bed $bedfile
# 在BED文件格式中，无论是正链（+）还是负链（-），起始位置总是小于终止位置。
# In the BED file format, the start position is always less than the end position for both positive (+) and negative (-) chains.
# 起始和终止位置是指定基因上的物理位置，而不是表达或翻译的方向。
# start and end positions specify the physical location of the gene, rather than the direction of expression or translation.
# starting site < stopped site in bed file


# -------------------------------------------------------------------------------------------
# 5. list of all genes found
print_fluorescent_yellow "     5. Extracting genes names: complete list of all genes found (universe.txt)"
cut -f 4 $bedfile > $universefile

# -------------------------------------------------------------------------------------------
# 6. strip the potential FASTA line breaks. creates genome_stripped.fa
print_fluorescent_yellow "     6. Removing potential FASTA line breaks (genome_stripped.fa)"
awk '/^>/ { if (NR!=1) print ""; printf "%s\n",$0; next;} \
    { printf "%s",$0;} \
    END { print ""; }'  $genomefile > $indexingOutputDir/genome_stripped.fa

# -------------------------------------------------------------------------------------------
# 7. promoter lenfths from promoters.bed
print_fluorescent_yellow "     7. Promoter lengths from genelines.bed (mRNA_lengths_all.txt)"
awk '{print $4 "\t" ($3 - $2)}' $indexingOutputDir/genelines.bed \
    > $indexingOutputDir/promoter_lengths_all.txt



# -------------------------------------------------------------------------------------------
# 8. promoter lenfths from promoters.bed
print_fluorescent_yellow "     8. mRNA isoform with longest lengths from genelines.bed (promoter_lengths_max.txt)"
# one gene has several isoforms. We select the mRNA with longest sequence.
awk '
BEGIN { FS=OFS="\t" }
{
    gene = $1
    sub(/\..*$/, "", gene)
    if (gene in max_length) {
        if ($2 > max_length[gene]) {
            max_length[gene] = $2
            isoform[gene] = $1
        }
    } else {
        max_length[gene] = $2
        isoform[gene] = $1
    }
}
END {
    for (gene in isoform) {
        print isoform[gene], max_length[gene]
    }
}' $indexingOutputDir/promoter_lengths_all.txt | sort > $indexingOutputDir/promoter_lengths_max.txt

rm -rf $indexingOutputDir/promoter_lengths_all.txt
mv $indexingOutputDir/promoter_lengths_max.txt $indexingOutputDir/promoter_lengths_all.txt


# -------------------------------------------------------------------------------------------
# 9. filters out the rows with NEGATIVE lengths
print_fluorescent_yellow "     9. Filtering out the rows of mRNA_lengths_all.txt with NEGATIVE lengths (promoter_lengths.txt)"
while read -r gene length; do
    # Check if the length is a positive number
    if (( length >= 0 )); then
        # Append rows with positive length to the output file
        echo "$gene $length" >> $indexingOutputDir/promoter_lengths.txt
    else
        # Append rows with negative length to the deleted file
        echo "$gene $length" >> $indexingOutputDir/promoter_lengths_deleted.txt
    fi
done < $indexingOutputDir/promoter_lengths_all.txt


# -------------------------------------------------------------------------------------------
# 10. find NEGATIVE genes
if [ -f "$indexingOutputDir/promoter_lengths_deleted.txt" ]; then
    print_fluorescent_yellow "    10. Finding genes with NEGATIVE promoter lengths (promoter_negative.txt)"
    cut -d " " \
        -f1  $indexingOutputDir/promoter_lengths_deleted.txt \
        > $indexingOutputDir/promoter_negative.txt
else
    print_fluorescent_yellow "    10. (skipped) Finding genes with NEGATIVE promoter lengths (promoter_negative.txt)"
fi

# -------------------------------------------------------------------------------------------
# 11. update gene list (no NEGATIVE genes)
print_fluorescent_yellow "    11. Updating gene list without NEGATIVE genes (universe.txt)";
cut -d " " -f1  $indexingOutputDir/promoter_lengths.txt > $universefile


# -------------------------------------------------------------------------------------------
# 12. filter promoter annotation with negative length
print_fluorescent_yellow "    12. Extracting bed coordinates with longest length (genelines.bed)"
awk 'NR==FNR{genes[$1]=1; next} genes[$4]'     \
    $universefile                              \
    $indexingOutputDir/genelines.bed           \
    > $indexingOutputDir/matched_promoterlines.bed


# -------------------------------------------------------------------------------------------
# 13. convert isoform name to it gene
print_fluorescent_yellow "    13. Converting isoform name to its gene (cleaned_matched_promoterlines.bed)"
awk -v OFS="\t" '{split($4, arr, "."); $4 = arr[1]; print $0}' \
    $indexingOutputDir/matched_promoterlines.bed \
    > $indexingOutputDir/cleaned_matched_promoterlines.bed

awk '{split($1, arr, "."); $1 = arr[1]; print $0}'    \
    $indexingOutputDir/promoter_lengths.txt           \
    > $indexingOutputDir/cleaned_promoter_lengths.txt
rm $indexingOutputDir/promoter_lengths.txt
mv $indexingOutputDir/cleaned_promoter_lengths.txt $indexingOutputDir/promoter_lengths.txt

awk -F'.' '{print $1}' $universefile > "${universefile}_cleaned"
rm -rf $universefile
mv "${universefile}_cleaned" $universefile

# -------------------------------------------------------------------------------------------
# 14. create promoters fasta
print_fluorescent_yellow "    14. Creating mRNA FASTA file (promoter_rought.fa)";
bedtools getfasta \
        -fi  $indexingOutputDir/genome_stripped.fa                \
        -bed $indexingOutputDir/cleaned_matched_promoterlines.bed \
        -fo  $indexingOutputDir/promoter_rought.fa                \
        -name -s

# -------------------------------------------------------------------------------------------
# 15. replace the id of each seq with gene names
print_fluorescent_yellow "    15. Replacing the id of each sequences' with gene names (promoter.fa)"
sed 's/::.*//g' $indexingOutputDir/promoter_rought.fa > $indexingOutputDir/promoter.fa

# -------------------------------------------------------------------------------------------
# 16. promoter.bg from promoter.fa
print_fluorescent_yellow "    16. fasta-get-markov estimates a Markov model from promoter.fa. (promoter.bg)"
fasta-get-markov $indexingOutputDir/promoter.fa > $indexingOutputDir/promoter.bg

# -------------------------------------------------------------------------------------------
# 17. IC.txt
print_fluorescent_yellow "    17. Generating information content (IC.txt)"
[ ! -d $indexingOutputDir/memefiles ] && mkdir $indexingOutputDir/memefiles
python3 $pmetroot/parse_memefile.py $memefile $indexingOutputDir/memefiles/
python3 $pmetroot/calculateICfrommeme_IC_to_csv.py \
    $indexingOutputDir/memefiles/ \
    $indexingOutputDir/IC.txt
rm -rf $indexingOutputDir/memefiles/*

# -------------------------------------------------------------------------------------------
# 18. individual motif files from user's meme file
print_fluorescent_yellow "    18. Spliting motifs into individual meme files (folder memefiles)"
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


# # -------------------------------- Transcript to gene --------------------------
# print_fluorescent_yellow "\n\nChangingg transcript names to gene names...\n";

# for file in "$indexingOutputDir"/fimohits/*.txt; do
#     awk 'BEGIN {FS=OFS="\t"} {sub(/\.[0-9]+$/, "", $2)} 1' "$file" > "${file}.tmp"
#     mv "${file}.tmp" "$file"
# done


# -------------------------------------------------------------------------------------------
# Deleting unnecessary files
if [[ $delete == "yes" || $delete == "YES" || $delete == "Y" || $delete == "y" ]]; then
    print_green "Deleting unnecessary files...\n\n"
    rm -f $indexingOutputDir/bedgenome.genome
    rm -f $indexingOutputDir/genome_stripped.fa
    rm -f $indexingOutputDir/genome_stripped.fa.fai
    mr -f $indexingOutputDir/memefiles
    rm -f $indexingOutputDir/promoter_rought.fa
    rm -f $indexingOutputDir/promoter.bg
    rm -f $indexingOutputDir/genelines.gff3
    rm -f $bedfile
    rm -f $indexingOutputDir/promoter_length_deleted.txt
    rm -r $indexingOutputDir/memefiles
    rm -f $indexingOutputDir/promoter.fa
    rm -f $indexingOutputDir/sorted.gff3
    rm -f $indexingOutputDir/pmetindex.log
    rm -f $indexingOutputDir/promoter_lengths_all.txt
    rm -f $indexingOutputDir/promoters_before_filter.bed
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
