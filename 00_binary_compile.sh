#!/bin/bash

print_fluorescent_yellow(){
    FLUORESCENT_YELLOW='\033[1;33m'
    NC='\033[0m' # No Color
    printf "${FLUORESCENT_YELLOW}$1${NC}\n"
}

print_fluorescent_yellow "Compiling PMET homotopic (index) binary..."
cd src/indexing

chmod a+x build.sh
bash build.sh

mv bin/pmetindex ../../scripts/


print_fluorescent_yellow "Compiling PMET heterotypic (pair) binary..."
cd ../pmetParallel

chmod a+x build.sh
bash build.sh

mv bin/pmetParallel ../../scripts/


print_fluorescent_yellow "Compiling PMET heterotypic (pair) binary..."
cd ../pmet

chmod a+x build.sh
bash build.sh

mv bin/pmet ../../scripts/