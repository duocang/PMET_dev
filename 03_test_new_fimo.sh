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

print_fluorescent_yellow "New fimo is running ..."
print_fluorescent_yellow "You can view results in 'results/new_fimo_out'"

mkdir -p results/new_fimo_out
start=$(date +%s)


# Give execute permission to all users for the file.
chmod a+x scripts/fimo

scripts/fimo                                       \
    --topk 5                                       \
    --topn 5000                                    \
    --text                                         \
    --no-qvalue                                    \
    --thresh 0.05                                  \
    --verbosity 3                                  \
    --oc results/new_fimo_out                      \
    --bgfile src/meme-5.5.3/fimo_test/promoters.bg \
    src/meme-5.5.3/fimo_test/memefiles/a.txt       \
    src/meme-5.5.3/fimo_test/promoters.fa          \
    src/meme-5.5.3/fimo_test/promoter_lengths.txt

end=$(date +%s)
time_taken=$((end - start))
print_red "Time taken: $time_taken seconds"

print_fluorescent_yellow "You can view results in 'results/new_fimo_out'"

print_green "done"
