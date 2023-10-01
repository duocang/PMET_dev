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
chmod a+x scripts/cpp_debug_needed/homotypic_intervals.sh

################################ Running FIMO #####################################
# check if fimo ready for PMETindex
print_fluorescent_yellow "Checking if fimo result (txt files) is ready for PMET index"

output=results/01_homotypic_intervals

directory=$output/fimo
mkdir -p $directory

txt_files=$(find "$directory" -name "*.txt")
if [ -n "$txt_files" ]; then
    print_green "   Yes, FIMO result (txt files) exist in $directory.\n\n"
else
    print_red   "   No, no FIMO result (txt files) exist in $directory.\n\n"
    print_green "Running fimo...\n"

    scripts/cpp_debug_needed/homotypic_intervals.sh  \
    -r scripts                                       \
    -o $output                                       \
    -n 5000                                          \
    -k 5                                             \
    -f 0.05                                          \
    -t 8                                             \
    data/homotypic_intervals/intervals.fa            \
    data/homotypic_intervals/motif_more.meme
fi

########################## Running pmet indexing ##################################
print_green "Running pmet index"

mkdir -p $output/fimohits
# run pmet index
scripts/pmetindex                                       \
    -f $output/fimo                 \
    -k 5                                                \
    -n 5000                                             \
    -p $output/promoter_lengths.txt \
    -o $output/

print_green "done"
print_fluorescent_yellow "\n\nYou many want to run '02_heterotypic_intervals.sh' now..."
