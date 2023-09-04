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

# fimo wht pmet index
print_fluorescent_yellow "Compiling FIMO with PMET homotopic (index) binary..."
cd src/meme-5.5.3

make distclean

currentDir=$(pwd)
chmod a+x ./configure
echo $currentDir/build
./configure --prefix=$currentDir/build  --enable-build-libxml2 --enable-build-libxslt
make
make install
cp build/bin/fimo ../../scripts/
make distclean
print_fluorescent_yellow "make distclean finished...\n"


# pmetindex
print_fluorescent_yellow "Compiling PMET homotopic (index) binary...\n"
cd ../indexing
chmod a+x build.sh
bash build.sh
mv bin/pmetindex ../../scripts/

# pmetParallel
print_fluorescent_yellow "Compiling PMET heterotypic (pair) binary...\n"
cd ../pmetParallel
chmod a+x build.sh
bash build.sh
mv bin/pmetParallel ../../scripts/

# pmet
print_fluorescent_yellow "Compiling PMET heterotypic (pair) binary...\n"
cd ../pmet
chmod a+x build.sh
bash build.sh
mv bin/pmet ../../scripts/


# Give execute permission to all users for the file.
chmod a+x scripts/pmetindex
chomd a+x scripts/pmetParallel
chmod a+x scripts/pmet
chmod a+x scripts/fimo

print_green "\n\npmet, pmetParallel, pmetindex and NEW fimo are ready in 'scripts' folder.\n"

print_green "DONE"
