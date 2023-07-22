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
        -t <threads>: Number of threads : 8
EOF
}

function error_exit() {
    echo "ERROR: $1" >&2
    usage
    exit 1
}

# set up defaults
topn=5000
maxk=5
promlength=1000
fimothresh=0.05
overlap="AllowOverlap"
utr="No"
gff3id='gene_id'
pmetroot="scripts"
progFile="progress/progress"
threads=8

# set up empty variables
outputdir=
genomefile=
gff3file=
memefile=
gene_input_file=

# deal with arguments
# if none, exit
if [ $# -eq 0 ]; then
    echo "No arguments supplied"  >&2
    usage
    exit 1
fi

while getopts ":r:i:o:n:k:p:f:g:v:u:t:" options; do
    case $options in
        r) echo "Full path of PMET_index:  $OPTARG" >&2
        pmetroot=$OPTARG;;
        i) echo "GFF3 feature identifier: $OPTARG" >&2
        gff3id=$OPTARG;;
        o) echo "Output directory for results: $OPTARG" >&2
        outputdir=$OPTARG;;
        n) echo "Top n promoter hits to take per motif: $OPTARG" >&2
        topn=$OPTARG;;
        k) echo "Top k motif hits within each promoter: $OPTARG" >&2
        maxk=$OPTARG;;
        p) echo "Promoter length: $OPTARG" >&2
        promlength=$OPTARG;;
        f) echo "Fimo threshold: $OPTARG" >&2
        fimothresh=$OPTARG;;
        v) echo "Remove promoter overlaps with gene sequences: $OPTARG" >&2
        overlap=$OPTARG;;
        u) echo "Include 5' UTR sequence?: $OPTARG" >&2
        utr=$OPTARG;;
        t) echo "Number of threads: $OPTARG" >&2
        threads=$OPTARG;;
        \?) echo "Invalid option: -$OPTARG" >&2
        exit 1;;
        :)  echo "Option -$OPTARG requires an argument." >&2
        exit 1;;
    esac
done

shift $((OPTIND - 1))
genomefile=$1
gff3file=$2
memefile=$3
gene_input_file=$4
universefile=$outputdir/universe.txt
bedfile=$outputdir/genelines.bed

[ ! -d $outputdir ] && mkdir $outputdir

echo "Preparing sequences...";

# --------------------------- sort annotaion by coordinates ----------------------------
$pmetroot/gff3sort/gff3sort.pl $gff3file > ${gff3file}temp

# ----------------------------------- promoters.bed ------------------------------------
if [[ "$(uname)" == "Linux" ]]; then
    grep -P '\tgene\t' ${gff3file}temp > $outputdir/genelines.gff3
elif [[ "$(uname)" == "Darwin" ]]; then
    grep '\tgene\t' ${gff3file}temp > $outputdir/genelines.gff3
else
    echo "Unsupported operating system."
fi
# grep -P '\tgene\t' ${gff3file}temp > $outputdir/genelines.gff3

# parse up the .bed for promoter extraction, 'gene_id'
python3 $pmetroot/parse_genelines.py $gff3id $outputdir/genelines.gff3 $bedfile
# rm $outputdir/genelines.gff3

# 在BED文件格式中，无论是正链（+）还是负链（-），起始位置总是小于终止位置。
# 这是因为起始和终止位置是指定基因或基因组特性在基因组上的物理位置，而不是表达或翻译的方向。
# starting site < stopped site in bed file
awk '$2 >= $3' $bedfile > $outputdir/invalid_genelines.bed
awk '$2 <  $3' $bedfile > temp.bed && mv temp.bed $bedfile

#list of all genes found
cut -f 4 $bedfile > $universefile



### Make bedgenome.genome and genome_stripped.fa
echo "Creating genome file...";

# ------------------------------------ extract genome --------------------------------------
# strip the potential FASTA line breaks. creates genome_stripped.fa
# python3 $pmetroot/strip_newlines.py $genomefile $outputdir/genome_stripped.fa
awk '/^>/ { if (NR!=1) print ""; printf "%s\n",$0; next;} \
    { printf "%s",$0;} \
    END { print ""; }'  $genomefile > $outputdir/genome_stripped.fa

# produces ouputdir/genome_stripped.fa
# create the .genome file which contains coordinates for each chromosome start
samtools faidx $outputdir/genome_stripped.fa
cut -f 1-2 $outputdir/genome_stripped.fa.fai > $outputdir/bedgenome.genome


### Use genelines.bed and bedgenome.genome to make promoters.bed
echo "Preparing promoter region information...";

# 在bedtools中，flank是一个命令行工具，用于在BED格式的基因组坐标文件中对每个区域进行扩展或缩短。
# 当遇到负链（negative strand）时，在区域的右侧进行扩展或缩短，而不是左侧。
bedtools flank \
    -l $promlength \
    -r 0-s -i $bedfile \
    -g $outputdir/bedgenome.genome \
    > $outputdir/promoters.bed
rm $outputdir/bedgenome.genome

echo "Removing overlapping promoter chunks (if true)..."
# remove overlapping promoter chunks
if [ $overlap == 'NoOverlap' ]; then
	echo "Removing overlaps";
	sleep 0.1
	bedtools subtract -a $outputdir/promoters.bed -b $bedfile > $outputdir/promoters2.bed
	mv $outputdir/promoters2.bed $outputdir/promoters.bed
fi
rm $bedfile


### Update promoters.bed using gff3file and universe file

# check that we have no split promoters. if so, keep the bit closer to the TSS
# Updates promoters.bed
python3 $pmetroot/assess_integrity.py $outputdir/promoters.bed
# possibly add 5' UTR
if [ $utr == 'Yes' ]; then
    echo "Adding UTRs...";
	python3 $pmetroot/parse_utrs.py $outputdir/promoters.bed ${gff3file}temp $universefile
fi

rm ${gff3file}temp

# -------------------- promoter_lengths file from promoters.bed ----------------------------
# python3 $pmetroot/parse_promoter_lengths.py $outputdir/promoters.bed $outputdir/promoter_lengths.txt
awk '{print $4 "\t" ($3 - $2)}' $outputdir/promoters.bed > $outputdir/promoter_lengths_all.txt

# filters out the rows with NEGATIVE lengths
# Process the data line by line
while read -r gene length; do
    # Check if the length is a positive number
    if (( length >= 0 )); then
        # Append rows with positive length to the output file
        echo "$gene $length" >> $outputdir/promoter_lengths.txt
    else
        # Append rows with negative length to the deleted file
        echo "$gene $length" >> $outputdir/promoter_length_deleted.txt
    fi
done < $outputdir/promoter_lengths_all.txt

# remove rows from "promoters.bed" that contain NEGATIVE genes (promoter_length_deleted.txt)
# get genes with negative length
cut -d " " -f1  $outputdir/promoter_length_deleted.txt > $outputdir/genes_negative.txt

grep -v -w -f \
    $outputdir/genes_negative.txt \
    $outputdir/promoters.bed \
    > $outputdir/filtered_promoters.bed

mv $outputdir/promoters.bed $outputdir/promoters_before_filter.bed
mv $outputdir/filtered_promoters.bed $outputdir/promoters.bed
rm $outputdir/genes_negative.txt

# update gene list (no NEGATIVE genes)
cut -d " " -f1  $outputdir/promoter_lengths.txt > $universefile

# ----------------------------------- promoters.fa -----------------------------------
## Make promoters.fa from promoters.bed and genome_stripped.fa
echo "Creating promoters file";

# get promoters
bedtools getfasta \
        -fi $outputdir/genome_stripped.fa \
        -bed $outputdir/promoters.bed \
        -s -fo $outputdir/promoters_rough.fa
rm $outputdir/genome_stripped.fa
rm $outputdir/genome_stripped.fa.fai


# creates promoters.fa
# replace the id of each seq with gene names
# python3 $pmetroot/parse_promoters.py $outputdir/promoters_rough.fa $outputdir/promoters.bed $outputdir/promoters.fa
# rm $outputdir/promoters.bed
# rm $outputdir/promoters_rough.fa
awk 'BEGIN{OFS="\t"} NR==FNR{a[NR]=$4; next} /^>/{$0=">"a[++i]} 1' \
    $outputdir/promoters.bed \
    $outputdir/promoters_rough.fa \
    > $outputdir/promoters.fa
rm $outputdir/promoters.bed
rm $outputdir/promoters_rough.fa

#------------------------- promoters.bg from promoters.fa ----------------------------
fasta-get-markov $outputdir/promoters.fa > $outputdir/promoters.bg
# FIMO barfs ALL the output. that's not good. time for individual FIMOs
# on individual MEME-friendly motif files too

echo "Processing motifs...";
# mkdir $outputdir/memefiles

# Make individual motif files from user's meme file
[ ! -d $outputdir/memefiles ] && mkdir $outputdir/memefiles
python3 $pmetroot/parse_memefile.py $memefile $outputdir/memefiles/

# ----------------------------------- IC.txt ---------------------------------------
python3 $pmetroot/calculateICfrommeme_IC_to_csv.py $outputdir/memefiles/ $outputdir/IC.txt

# -------------------------------- Run fimo and pmetindex --------------------------
echo "Processing PMET indexing ...";
# Create a fimo hits file form each motif using promoters.bg and promoters.fa -----------------
[ ! -d $outputdir/fimo ] && mkdir $outputdir/fimo
[ ! -d $outputdir/fimohits ] && mkdir $outputdir/fimohits

# shopt -s nullglob # prevent loop produncing '*.txt'

# Run fimo and pmetindex on each mitif (parallel version)
runFimoIndexing () {
    memefile=$1
    outputdir=$2
    fimothresh=$3
    pmetroot=$4
    maxk=$5
    topn=$6
    # echo $memefile
    filename=`basename $memefile .txt`
    echo $filename

    mkdir -p $outputdir/fimo/$filename

    fimo \
        --no-qvalue \
        --text \
        --thresh $fimothresh \
        --verbosity 1 \
        --bgfile $outputdir/promoters.bg\
        $memefile \
        $outputdir/promoters.fa \
        > $outputdir/fimo/$filename/$filename.txt
    $pmetroot/pmetindex \
        -f $outputdir/fimo/$filename \
        -k $maxk \
        -n $topn \
        -o $outputdir \
        -p $outputdir/promoter_lengths.txt > $outputdir/pmetindex.log
    rm -rf $outputdir/fimo/$filename
}
export -f runFimoIndexing

find $outputdir/memefiles -name \*.txt \
    | parallel --jobs=$threads \
        "runFimoIndexing {} $outputdir $fimothresh $pmetroot $maxk $topn"

numfiles=$(ls -l $outputdir/memefiles/*.txt | wc -l)
echo $numfiles" motifs found"

echo "Delete unnecessary files"
rm -r $outputdir/memefiles
rm $outputdir/promoters.bg
rm $outputdir/promoters.fa
# touch ${outputdir}_FLAG



# next stage needs the following inputs

#   promoter_lengths.txt        made by parse_promoter_lengths.py from .bed file
#   bimnomial_thresholds.txt    made by PMETindex
#   IC.txt                      made by calculateICfrommeme.py from meme file
#   gene input file             supplied by user
