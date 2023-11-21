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

print_green_no_br(){
    GREEN='\033[0;32m'
    NC='\033[0m' # No Color
    printf "${GREEN}$1${NC}"
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
print_fluorescent_yellow_no_br(){
    FLUORESCENT_YELLOW='\033[1;33m'
    NC='\033[0m' # No Color
    printf "${FLUORESCENT_YELLOW}$1${NC}"
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

echo -e "\n\n"
print_middle "The purpose of this script is to                                      \n"
print_middle "  1. assign execute permissions to all users for bash and perl files    "
print_middle "  2. compile binaries needed by Shiny app                               "
print_middle "  3. install python package                                             "
print_middle "  4. check needed tools                                               \n"
print_middle "                                                                      \n"

if [ -d .git ]; then
    git config core.fileMode false
fi

############################ 1. assign execute permissions #############################

print_green_no_br "\n1. Would you like to assign execute permissions to all users for bash and perl files? [Y/n]: "
read -p "" answer
answer=${answer:-Y} # Default to 'Y' if no input provided

if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then
    # 遍历 scripts 目录及其所有子目录中的 .sh 和 .pl 文件
    find . -type f \( -name "*.sh" -o -name "*.pl" \) -exec chmod a+x {} \;
else
    print_orange "No assignment"
fi


################################## 2. compile binary #################################
print_green_no_br "\n2. Would you like to compile binaries? [y/N]:"
read -p " " answer
answer=${answer:-N} # Default to 'N' if no input provided

if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then

    print_fluorescent_yellow "Compiling... It takes minutes."

    rm -f scripts/pmetindex
    rm -f scripts/pmetParallel
    rm -f scripts/pmet
    rm -f scripts/fimo

    ############################# 2.1 fimo with pmet index
    print_orange "2.1 Compiling FIMO with PMET homotopic (index) integtated..."
    cd src/meme-5.5.3

    make distclean > /dev/null 2>&1

    currentDir=$(pwd)

    if [ -d "$currentDir/build" ]; then
        rm -rf "$currentDir/build"
    fi
    mkdir -p $currentDir/build

    # update congifure files according to different system
    aclocal  > /dev/null 2>&1
    automake > /dev/null 2>&1

    chmod a+x ./configure
    ./configure --prefix=$currentDir/build  --enable-build-libxml2 --enable-build-libxslt  > /dev/null 2>&1
    make          > /dev/null 2>&1
    make install  > /dev/null 2>&1
    cp build/bin/fimo ../../scripts/
    make distclean > /dev/null 2>&1
    rm -rf build
    # print_orange "make distclean finished...\n"


    ################################### 2.2 pmetindex
    print_orange "2.2 Compiling PMET homotypic (index) binary..."
    cd ../indexing

    bash build.sh > /dev/null 2>&1
    mv bin/pmetindex ../../scripts/
    rm -rf bin/*

    ################################## 2.3 pmetParallel
    print_orange "2.3 Compiling PMET heterotypic (pair) binary..."
    cd ../pmetParallel

    bash build.sh > /dev/null 2>&1
    mv bin/pmetParallel ../../scripts/
    rm -rf bin/*

    # pmet
    print_orange "2.4 Compiling PMET heterotypic (pair) binary..."
    cd ../pmet

    bash build.sh > /dev/null 2>&1
    mv bin/pmet ../../scripts/
    rm -rf bin/*

    # back to home directory
    cd ../..

    ################### 2.4 Check if the compilation was successful
    exists=""
    not_exists=""

    for file in scripts/pmetindex scripts/pmetParallel scripts/pmet scripts/fimo; do
        if [ -f "$file" ]; then
            exists="$exists\n    $file"
        else
            not_exists="$not_exists\n    $file"
        fi
    done

    if [ ! -z "$exists" ]; then
        echo -e "\n"
        print_green "Compilation Success:$exists"
    fi
    if [ ! -z "$not_exists" ]; then
        echo -e "\n"
        print_red "Compilation Failure:$not_exists"
    fi

    ############# 2.5 Give execute permission to all users for the file
    chmod a+x scripts/pmetindex
    chmod a+x scripts/pmetParallel
    chmod a+x scripts/pmet
    chmod a+x scripts/fimo
else
    print_orange "No tools compiled"
fi


############################# 3. install python packages ##############################
print_green_no_br "\n3. Would you like to install python packages? [y/N]: "
read -p " " answer
answer=${answer:-N} # Default to 'N' if no input provided

if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then
    print_orange "Installing python packages... It takes minutes."
    pip install numpy     > /dev/null 2>&1 || echo "Failed to install numpy"
    pip install pandas    > /dev/null 2>&1 || echo "Failed to install pandas"
    pip install scipy     > /dev/null 2>&1 || echo "Failed to install scipy"
    pip install bio       > /dev/null 2>&1 || echo "Failed to install bio"
    pip install biopython > /dev/null 2>&1 || echo "Failed to install biopython"
else
    print_orange "No python packages installed"
fi

################################ 4. check needed tools #################################
print_green "\n7. Checking the existence of GNU Parallel, bedtools, samtools and MEME Suite "
# List of tools to check
tools=("parallel" "bedtools" "samtools" "fimo")

# Assume all tools are installed until one is not found
all_tools_found=true

# Iterate over each tool and check if it is installed
for tool in "${tools[@]}"; do
    if ! command -v $tool &> /dev/null
    then
        print_red "$tool could not be found"
        all_tools_found=false
        # Optionally exit or continue to check other tools
        # exit 1
    fi
done

# If all tools were found, print a positive message
if $all_tools_found; then
    print_green "All tools were found!"
else
    print_red "Please install them and rerun the script"
fi

print_green "\nDONE"

