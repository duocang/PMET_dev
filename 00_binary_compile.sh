#!/bin/bash

print_fluorescent_yellow(){
    FLUORESCENT_YELLOW='\033[1;33m'
    NC='\033[0m' # No Color
    printf "${FLUORESCENT_YELLOW}$1${NC}\n"
}

# pmetindex
print_fluorescent_yellow "Compiling PMET homotopic (index) binary..."
cd src/indexing

chmod a+x build.sh
bash build.sh

mv bin/pmetindex ../../scripts/

# pmetParallel
print_fluorescent_yellow "Compiling PMET heterotypic (pair) binary..."
cd ../pmetParallel

chmod a+x build.sh
bash build.sh

mv bin/pmetParallel ../../scripts/

# pmet
print_fluorescent_yellow "Compiling PMET heterotypic (pair) binary..."
cd ../pmet

chmod a+x build.sh
bash build.sh

mv bin/pmet ../../scripts/


# fimo wht pmet index
print_fluorescent_yellow "Compiling FIMO with PMET homotopic (index) binary..."
cd ../meme-5.5.3

currentDir=$(pwd)
chmod a+x ./configure

echo $currentDir/build

./configure --prefix=$currentDir/build  --enable-build-libxml2 --enable-build-libxslt
make
make install

cp build/bin/fimo ../../scripts/
