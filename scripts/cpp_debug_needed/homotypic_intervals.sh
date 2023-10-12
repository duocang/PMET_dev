#!/bin/bash
set -e

# disable warnings
set -o errexit
set -o pipefail

# 22.1.18 Charlotte Rich
# last edit: 7.2.18 - removed the make 1 big fimohits files

# cl_index_wrapper.sh
# mac -> server Version differences
# ggrep = grep

# Called when user selects 'Genomic Intervals'
# Input files are genomic intevals fasta file, meme file location, gene clusters file
# Other inputs N and k

function usage () {
    cat >&2 <<EOF
        USAGE: PMETindexgenome [options] <genome> <gff3> <memefile>

        Creates PMET index for Paired Motif Enrichment Test using genome files.
        Required arguments:
        -r <PMETindex_path>	: Full path of python scripts called from this file. Required.
        -i <gff3_identifier> : gene identifier in gff3 file e.g. gene_id=

        Optional arguments:
        -o <output_directory> : Output directory for results.
        -n <topn>	: How many top promoter hits to take per motif. Default=5000
        -k <max_k>	: Maximum motif hits allowed within each promoter.  Default: 5
        -f <fimo_threshold> : Specify a minimum quality for hits matched by fimo. Default: 0.05
        -t <threads>: Number of threads. Default: 4
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


# set up arguments
topn=5000
maxk=5
fimothresh=0.05
pmetroot="scripts"
threads=8

outputdir=
genomefile=
memefile=

# check if arguments have been specified
if [ $# -eq 0 ]
then
    echo "No arguments supplied"  >&2
    usage
    exit 1
fi

# bring in arguments
while getopts ":r:o:k:n:f:t:" options; do
    case $options in
        r) print_white "Full path of PMET_index                : "; print_orange "$OPTARG" >&2
        pmetroot=$OPTARG;;
        o) print_white "Output directory for results           : "; print_orange "$OPTARG" >&2
        outputdir=$OPTARG;;
        n) print_white "Top n promoter hits to take per motif  : "; print_orange "$OPTARG" >&2
        topn=$OPTARG;;
        k) print_white "Top k motif hits within each promoter  : "; print_orange "$OPTARG" >&2
        maxk=$OPTARG;;
        f) print_white "Fimo threshold                         : "; print_orange "$OPTARG" >&2
        fimothresh=$OPTARG;;
        t) print_white "Number of threads                      : "; print_orange "$OPTARG" >&2
        threads=$OPTARG;;
        \?) print_red "Invalid option: -$OPTARG" >&2
        exit 1;;
        :)  print_red "Option -$OPTARG requires an argument." >&2
        exit 1;;
    esac
done

# rename input file variable
shift $((OPTIND - 1))
genomefile=$1
memefile=$2

print_white "Genomic interval file                  : "; print_orange $genomefile
print_white "Motif meme file                        : "; print_orange $memefile

[ ! -d $outputdir ] && mkdir $outputdir
# cd $outputdir

print_green "Preparing sequences...";

# final pmet binary requires the universe file. Need to create this if validation scrip didnt.
# In promoters version, this is initially all genes in gff3 file. This version is used to add UTRs if
# requested, but any genes not in promoter_lengths file are filtered out before we get to PMET binary stage
# In this version we can just take a copy of all IDs in promoter lengths as we dont to UTR stuff

universefile=$outputdir/universe.txt

if [[ ! -f "$universefile" || ! -f "$outputdir/promoter_lengths.txt" ]]; then
    # should have been done by consistency checker
    # *** ADD THE DEPUPLICATION OF THE FASTA FILE HERE ****
    python3 $pmetroot/deduplicate.py \
            $genomefile \
            $outputdir/no_duplicates.fa

    # generate the promoter lengths file from the fasta file
    python3 $pmetroot/parse_promoter_lengths_from_fasta.py \
            $outputdir/no_duplicates.fa \
            $outputdir/promoter_lengths.txt
    # rm -f $outputdir/no_duplicates.fa

    cut -f 1  $outputdir/promoter_lengths.txt > $universefile
fi

# now we can actually FIMO our way to victory
fasta-get-markov $genomefile > $outputdir/genome.bg
# FIMO barfs ALL the output. that's not good. time for individual FIMOs
# on individual MEME-friendly motif files too

print_green "Processing motifs...";

### Make motif  files from user's meme file
[ ! -d $outputdir/memefiles ] && mkdir $outputdir/memefiles

python3 $pmetroot/parse_memefile.py \
        $memefile \
        $outputdir/memefiles/

### creates IC.txt tsv file from, motif files
python3 $pmetroot/calculateICfrommeme_IC_to_csv.py \
        $outputdir/memefiles/ \
        $outputdir/IC.txt

### Create a fimo hits file form each motif file
[ ! -d $outputdir/fimo ] && mkdir $outputdir/fimo
[ ! -d $outputdir/fimohits ] && mkdir $outputdir/fimohits

shopt -s nullglob # prevent loop produncing '*.txt'

numfiles=$(ls -l $outputdir/memefiles/*.txt | wc -l)
print_orange "$numfiles motifs found"

n=0
# paralellise this loop
for memefile in $outputdir/memefiles/*.txt; do
    let n=$n+1
    fimofile=`basename $memefile`
    echo $fimofile

    fimo    --text                        \
            --thresh $fimothresh          \
            --verbosity 1                 \
            --bgfile $outputdir/genome.bg \
            $memefile                     \
            $genomefile                   \
            > $outputdir/fimo/$fimofile &
    [ `expr $n % $threads` -eq 0 ] && wait
done

# echo "Delete unnecessary files"
# rm -r $outputdir/memefiles
# rm $outputdir/genome.bg


# next stage needs the following inputs

#   promoter_lengths.txt        made by parse_promoter_lengths.py from .bed file
#   bimnomial_thresholds.txt    made by PMETindex
#   IC.txt                      made by calculateICfrommeme.py from meme file
#   gene input file             supplied by user
