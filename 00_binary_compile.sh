#!/bin/bash

if [ true ]; then
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
    print_orange_no_br(){
        ORANGE='\033[0;33m'
        NC='\033[0m' # No Color
        printf "${ORANGE}$1${NC}"
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
fi
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


############################ 1. linux dependences #############################
if [ true ]; then
    print_green "\n1.Adding some necessary dependencies, please wait a moment..."
    print_orange_no_br "  Do you want to install Ubuntu dependencies? [y/N]: "
    read dependency
    dependency=${dependency:-N} # Default to 'N' if no input provided

    if [ "$dependency" == "Y" ] || [ "$dependency" == "y" ]; then
        print_orange "Installing dependencies..."
        # Check for MacOS
        if [[ "$OSTYPE" == "darwin"* ]]; then
            echo "This is MacOS. No support!"
        # Check for Linux
        elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
            # Source the os-release to get OS information
            if [ -f /etc/os-release ]; then
                # Load the os-release file
                . /etc/os-release

                if [[ "$ID" == "ubuntu" || "$ID_LIKE" == "debian" ]]; then
                    echo "This is a Debian/Ubuntu based distribution. Version: $PRETTY_NAME"

                    sudo apt update && sudo apt upgrade -y > /dev/null
                    sudo apt -y install openjdk-21-jdk     > /dev/null 2>&1 || print_red "    Failed to install openjdk-21-jdk"

                    # export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
                    # export PATH=$PATH:$JAVA_HOME/bin
                    # 询问用户是否处于调试模式

                    sudo add-apt-repository -y ppa:alex-p/tesseract-ocr-devel  > /dev/null 2>&1 || print_red "    Failed to add ppa:alex-p/tesseract-ocr-devel"
                    sudo apt -y install dpkg                > /dev/null 2>&1 || print_red "    Failed to install dpkg"
                    sudo apt -y install zip                 > /dev/null 2>&1 || print_red "    Failed to install zip"
                    sudo apt -y install build-essential     > /dev/null 2>&1 || print_red "    Failed to install build-essential"
                    sudo apt -y install zlib1g-dev          > /dev/null 2>&1 || print_red "    Failed to install zlib1g-dev"
                    sudo apt -y install libgdbm-dev         > /dev/null 2>&1 || print_red "    Failed to install libgdbm-dev"
                    sudo apt -y install autotools-dev       > /dev/null 2>&1 || print_red "    Failed to install autotools-dev"
                    sudo apt -y install automake            > /dev/null 2>&1 || print_red "    Failed to install automake"
                    sudo apt -y install nodejs              > /dev/null 2>&1 || print_red "    Failed to install nodejs"
                    sudo apt -y install npm                 > /dev/null 2>&1 || print_red "    Failed to install npm"
                    sudo apt -y install libxml2-dev         > /dev/null 2>&1 || print_red "    Failed to install libxml2-dev"
                    sudo apt -y install libxml2             > /dev/null 2>&1 || print_red "    Failed to install libxml2"
                    sudo apt -y install libxslt-dev         > /dev/null 2>&1 || print_red "    Failed to install libxslt-dev"
                    sudo apt -y install cmake               > /dev/null 2>&1 || print_red "    Failed to install cmake"
                    sudo apt -y install libjpeg-dev         > /dev/null 2>&1 || print_red "    Failed to install libjpeg-dev "
                    sudo apt -y install librsvg2-dev        > /dev/null 2>&1 || print_red "    Failed to install librsvg2-dev"
                    sudo apt -y install libpoppler-cpp-dev  > /dev/null 2>&1 || print_red "    Failed to install libpoppler-cpp-dev"
                    sudo apt -y install freetype2-demos     > /dev/null 2>&1 || print_red "    Failed to install freetype2-demos"
                    sudo apt -y install cargo               > /dev/null 2>&1 || print_red "    Failed to install cargo"
                    sudo apt -y install libharfbuzz-dev     > /dev/null 2>&1 || print_red "    Failed to install libharfbuzz-dev"
                    sudo apt -y install libfribidi-dev      > /dev/null 2>&1 || print_red "    Failed to install libfribidi-dev"
                    sudo apt -y install libtesseract-dev    > /dev/null 2>&1 || print_red "    Failed to install libtesseract-dev"
                    sudo apt -y install libleptonica-dev    > /dev/null 2>&1 || print_red "    Failed to install libleptonica-dev"
                    sudo apt -y install tesseract-ocr-eng   > /dev/null 2>&1 || print_red "    Failed to install libfribidi-dev"
                    sudo apt -y install libmagick++-dev     > /dev/null 2>&1 || print_red "    Failed to install libmagick++-dev"
                    sudo apt -y install libavfilter-dev     > /dev/null 2>&1 || print_red "    Failed to install libavfilter-dev"
                    sudo apt -y install libncurses5-dev     > /dev/null 2>&1 || print_red "    Failed to install libncurses5-dev"
                    sudo apt -y install libncursesw5-dev    > /dev/null 2>&1 || print_red "    Failed to install libncursesw5-dev"
                    sudo apt -y install libbz2-dev          > /dev/null 2>&1 || print_red "    Failed to install libbz2-dev"
                    sudo apt -y install libpcre2-dev        > /dev/null 2>&1 || print_red "    Failed to install libpcre2-dev"
                    sudo apt -y install libreadline-dev     > /dev/null 2>&1 || print_red "    Failed to install libreadline-dev"
                    sudo apt -y install libssl-dev          > /dev/null 2>&1 || print_red "    Failed to install libssl-dev"
                    sudo apt -y install glibc-source        > /dev/null 2>&1 || print_red "    Failed to install glibc-source"
                    sudo apt -y install libssl-dev          > /dev/null 2>&1 || print_red "    Failed to install libssl-dev"
                    sudo apt -y install libstdc++6          > /dev/null 2>&1 || print_red "    Failed to install libstdc++6"

                    echo "deb http://security.ubuntu.com/ubuntu focal-security main" | sudo tee /etc/apt/sources.list.d/focal-security.list
                    sudo apt update
                    sudo apt -y install libssl1.1             > /dev/null 2>&1 || print_red "    Failed to install libssl1.1"
                    sudo apt -y install libcurl4-openssl-dev  > /dev/null 2>&1 || print_red "    Failed to install libcurl4-openssl-dev"
                    sudo apt -y install libopenblas-dev       > /dev/null 2>&1 || print_red "    Failed to install libopenblas-dev"
                    sudo apt -y install gfortran              > /dev/null 2>&1 || print_red "    Failed to install gfortran"
                    sudo update-alternatives --config libblas.so.3-$(arch)-linux-gnu
                    sudo apt -y install ufw > /dev/null 2>&1 || print_red "    Failed to install ufw"

                elif [[ "$ID" == "fedora" || "$ID_LIKE" == "fedora" || "$ID" == "rhel" || "$ID_LIKE" == "rhel" ]]; then
                    echo "This is a Red Hat based distribution. Version: $PRETTY_NAME"
                else
                    echo "This is a Linux distribution but not specifically Debian/Ubuntu or Red Hat based. Version: $PRETTY_NAME"
                fi # if [[ "$ID" == "ubuntu" || "$ID_LIKE" == "debian" ]]; then
            else # if [ -f /etc/os-release ]; then
                echo "Linux distribution, but /etc/os-release not found."
            fi # if [ -f /etc/os-release ]; then
        # Check for Windows (Cygwin or MinGW)
        elif [[ "$OSTYPE" == "cygwin" ]]; then
            echo "This is Windows with Cygwin. No support!"
        elif [[ "$OSTYPE" == "msys" ]]; then
            echo "This is Windows with MinGW. No support!"
        # Other unknown OS
        else
            echo "Unknown Operating System. No support!"
        fi
    else # if [ "$dependency" == "Y" ] || [ "$dependency" == "y" ]; then
        print_fluorescent_yellow "    No Linux dependencies installed"
    fi
fi


############################ 2. assign execute permissions #############################

print_green_no_br "\n2. Would you like to assign execute permissions to all users for bash and perl files? [Y/n]: "
read -p "" answer
answer=${answer:-Y} # Default to 'Y' if no input provided

if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then
    # 遍历 scripts 目录及其所有子目录中的 .sh 和 .pl 文件
    find . -type f \( -name "*.sh" -o -name "*.pl" \) -exec chmod a+x {} \;
else
    print_orange "No assignment"
fi

################################## 3. install R #################################
# 检查 R 是否已安装
if [ true ]; then
    if ! command -v R >/dev/null 2>&1; then
        print_green " \n3. R is not installed. Installing R..."

        print_green_no_br "\n1. Would you like to install R (/opt/R/R-4.3.2/bin/R) as Root? [y/N]: "
        read -p "" answer
        answer=${answer:-N} # Default to 'Y' if no input provided

        if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then
            tool_dir=R-4.3.2
            tool_link=https://cran.rstudio.com/src/base/R-4/R-4.3.2.tar.gz

            wget $tool_link          > /dev/null 2>&1 || echo "    Failed to download R tar file"
            tar -xzvf R-4.3.2.tar.gz > /dev/null 2>&1 || echo "    Failed to unzip tar file"
            rm R-4.3.2.tar.gz

            cd $tool_dir
            ./configure                    \
                --prefix=/opt/R/$tool_dir  \
                --enable-R-shlib=yes        \
                --enable-memory-profiling  \
                --with-blas                \
                --with-lapack              \
                --with-readline=yes        \
                --with-x=no  > /dev/null 2>&1 || echo "    Failed to install configure"

            make              > /dev/null 2>&1 || echo "    Failed to make"
            sudo make install > /dev/null 2>&1 || echo "    Failed to make install"

            sudo ln -s /opt/R/$tool_dir/bin/R /usr/local/bin/R
            sudo ln -s /opt/R/$tool_dir/bin/Rscript /usr/local/bin/Rscript
            cd ..
            rm -rf $tool_dir
        else
            print_orange "    Please install R later."
        fi #if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then

    else
        print_green "\n2. R is already installed."
    fi
fi

################################## 4. install python #################################
# 检查 pyhon 是否已安装
if [ true ]; then
    # 检测并安装 python3
    if ! dpkg -l | grep -q "^ii\s*python3\s"; then
        print_green_no_br "\n4. Would you like to install python as Root? [y/N]: "
        read -p "" answer
        answer=${answer:-N} # Default to 'Y' if no input provided

        if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then
            print_green "  Installing python3..."
            sudo apt -y install python3 > /dev/null     2>&1 || print_red "    Failed to install python3"
            sudo apt -y install python3-pip > /dev/null 2>&1 || print_red "    Failed to install python3-pip"
        else
            print_orange "    Please install python later!"
        fi
    else
        print_green "\n3. python3 is already installed."
    fi
fi


################################## 5. compile binary #################################
print_green_no_br "\n5. Would you like to compile PMET and PMETindex binaries? [y/N]:"
if [ true ]; then
    read -p " " answer
    answer=${answer:-N} # Default to 'N' if no input provided
    if [[ "$answer" =~ ^[Yy]$ ]]; then
        if [ -f "scripts/pmetindex" ]; then
            chmod a+x scripts/pmetindex
            chmod a+x scripts/pmet
            chmod a+x scripts/fimo
            while true; do
                print_orange       "  PMET and PMETindex exist."
                print_orange_no_br "  Do you want to recompile PMET and PMETindex? [y/N]: "
                read -p " " recompile
                recompile=${recompile:-N} # Default to 'N' if no input provided
                # 检查用户输入是否为 Y, y, N, 或 n
                if [[ "$recompile" =~ ^[YyNn]$ ]]; then
                    break  # 如果输入有效，跳出循环
                else
                    print_red "    Invalid input. Please enter Y/y for yes or N/n for no."
                fi
            done
        else
            print_red "  PMETindex and PMET not found. Compile will start..."
            recompile=y
        fi

        if [[ "$recompile" =~ ^[Yy]$ ]]; then
            print_orange_no_br "  Do you still want to combile PMET and PMETindex? [y/N]: "
            read -p "" compile_still
        fi
        compile_still=${compile_still:-N}


        if [[ "$recompile" =~ ^[Yy]$ &&  "$compile_still" =~ ^[Yy]$ ]]; then

            print_fluorescent_yellow "Compiling... It takes minutes."

            rm -f scripts/pmetindex
            rm -f scripts/pmetParallel
            rm -f scripts/pmet
            rm -f scripts/fimo

            ############################# 2.1 fimo with pmet index
            print_orange "2.1 Compiling FIMO with PMET homotopic (index) integtated..."
            cd src/meme-5.5.3

            make distclean > /dev/null 2>&1 #|| print_red "    Failed to make distclean"
            aclocal        > /dev/null 2>&1 || print_red "    Failed to aclocal "
            automake       > /dev/null 2>&1 || print_red "    Failed to automake"

            currentDir=$(pwd)
            if [ -d "$currentDir/build" ]; then
                rm -rf "$currentDir/build"
            fi
            mkdir -p $currentDir/build

            chmod a+x ./configure
            ./configure --prefix=$currentDir/build  --enable-build-libxml2 --enable-build-libxslt  > /dev/null 2>&1 || print_red "    Failed to configure"
            make           > /dev/null 2>&1 || print_red "    Failed to make"
            make install   > /dev/null 2>&1 || print_red "    Failed to make install"
            cp build/bin/fimo ../../scripts/
            make distclean > /dev/null 2>&1 || print_red "    Failed to make distclean"
            rm -rf build
            # print_orange "make distclean finished...\n"


            ################################### 2.2 pmetindex
            print_orange "2.2 Compiling PMET homotypic (index) binary..."
            cd ../indexing

            bash build.sh > /dev/null 2>&1
            mv bin/pmetindex ../../scripts/
            rm -rf bin/*

            ################################## 2.3 pmetParallel
            print_orange "2.3 Compiling PMET heterotypic parallel (pair) binary..."
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
        fi # if [[ "$recompile" =~ ^[Yy]$ &&  "$compile_still" =~ ^[Yy]$ ]]; then
    else
        print_orange "    No tools compiled"
    fi # if [[ "$answer" =~ ^[Yy]$ ]]; then
fi

############################# 6. install python packages ##############################
if [ true ]; then
    print_green_no_br "\n6. Would you like to install python packages? [y/N]: "
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
        print_orange "    No python packages installed"
    fi
fi

################################ 7. check needed tools #################################
if [ true ]; then
    print_green "\n7. Checking the execution of GNU Parallel, bedtools, samtools and MEME Suite"
    tools=("parallel --version" "bedtools --version" "samtools --version" "fimo -h") # List of tools and their version check commands
    missing_tools=()  # Initialize an empty array to store tools that cannot be executed
    all_tools_found=true

    # Iterate over each tool and check if it can be executed
    for tool in "${tools[@]}"; do
        # Use 'command -v' to check if the tool command is available in the PATH
        if ! command -v $tool &> /dev/null; then
            echo "    $tool could not be found"
            all_tools_found=false
            missing_tools+=($tool)  # Add the tool to the missing tools array
        fi
    done

    # If all tools were found, print a positive message
    if $all_tools_found; then
        print_orange "    All tools were found!"
    else
        print_fluorescent_yellow_no_br "Would you like to install missing tools? [Y/n]: "
        read -p " " answer
        answer=${answer:-Y} # Default to 'N' if no input provided

        if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then
            # Install the missing tools
            for tool in "${missing_tools[@]}"; do
                if [ "$tool" == "parallel" ]; then
                    if [[ "$(uname)" == "Darwin" ]]; then
                        brew install parallel > /dev/null 2>&1 || print_red "    Failed to install bedtools"
                    elif [[ -f /etc/os-release ]]; then
                        . /etc/os-release
                        case $ID in
                            centos|fedora|rhel)
                                yum install parallel  > /dev/null 2>&1 || print_red "    Failed to install bedtools"
                                ;;
                            debian|ubuntu)
                                sudo apt update               > /dev/null 2>&1 || print_red "    Failed to apt update"
                                sudo apt -y install parallel  > /dev/null 2>&1 || print_red "    Failed to install bedtools"
                                ;;
                            *)
                                echo "Unsupported Linux distribution"
                                ;;
                        esac
                        parallel --citation              > /dev/null 2>&1 || print_red "    Failed to run parallel --citation"
                    else
                        echo "Unsupported OS"
                    fi
                fi

                if [ "$tool" == "bedtools" ]; then
                    if [[ "$(uname)" == "Darwin" ]]; then
                        brew install bedtools > /dev/null 2>&1 || print_red "    Failed to install bedtools"
                    elif [[ -f /etc/os-release ]]; then
                        . /etc/os-release
                        case $ID in
                            centos|fedora|rhel)
                                yum install BEDTools  > /dev/null 2>&1 || print_red "    Failed to install bedtools"
                                ;;
                            debian|ubuntu)
                                sudo apt update               > /dev/null 2>&1 || print_red "    Failed to apt update"
                                sudo apt -y install bedtools  > /dev/null 2>&1 || print_red "    Failed to install bedtools"
                                ;;
                            *)
                                echo "Unsupported Linux distribution"
                                ;;
                        esac
                    else
                        echo "Unsupported OS"
                    fi
                fi

                if [ "$tool" == "samtools" ]; then
                    print_green_no_br "\n1. Would you like to install samtools (~/tools/samtools-1.17/)? [y/N]: "
                    read -p "" answer
                    answer=${answer:-Y} # Default to 'Y' if no input provided

                    if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then
                        # 创建 ~/tools 目录，如果运行脚本的用户不是 shiny，则需要修改路径
                        mkdir -p ~/tools

                        tool_dir=samtools-1.17
                        tool_link=https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2

                        # 下载并解压 samtools
                        wget $tool_link -O $tool_dir.tar.bz2 > /dev/null 2>&1 || echo "    Failed to download samtools tar file"
                        tar -xjf $tool_dir.tar.bz2           > /dev/null 2>&1 || echo "    Failed to unzip tar file"
                        rm $tool_dir.tar.bz2

                        # 进入 samtools 目录并编译
                        cd $tool_dir && mkdir -p build
                        ./configure --prefix=$(pwd)/build    > /dev/null 2>&1 || echo "    Failed to configure"
                        make                                 > /dev/null 2>&1 || echo "    Failed to make"
                        sudo make install                    > /dev/null 2>&1 || echo "    Failed to make install"

                        # 返回到原始目录并将 samtools 目录移动到 ~/tools 下
                        cd ..
                        mv $tool_dir ~/tools/

                        # 更新 ~/.bashrc 和 ~/.zshrc，添加 samtools 到 PATH
                        echo "export PATH=~/tools/$tool_dir/build/bin:\$PATH" >> ~/.bashrc
                        echo "export PATH=~/tools/$tool_dir/build/bin:\$PATH" >> ~/.zshrc
                        echo "export PATH=~/tools/$tool_dir/build/bin:\$PATH" >> ~/.bash_profile
                        sudo chown -R shiny:shiny-apps ~/tools/$tool_dir/build/
                        print_green "    samtools successfully installed"
                    else
                        print_red "    Please install samtools later!"
                    fi # if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then
                fi # if [ "$tool" == "samtools" ]; then

                if [ "$tool" == "fimo" ]; then
                    print_green_no_br "\n1. Would you like to install MEME suite (~/tools/meme-5.5.2/)? [y/N]: "
                    read -p "" answer
                    answer=${answer:-Y} # Default to 'Y' if no input provided

                    if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then
                        mkdir -p ~/tools

                        tool_dir=meme-5.5.2
                        tool_link=https://meme-suite.org/meme/meme-software/5.5.2/meme-5.5.2.tar.gz

                        wget $tool_link           > /dev/null 2>&1 || print_red "    Failed to download fimo tar file"
                        tar zxf meme-5.5.2.tar.gz > /dev/null 2>&1 || print_red "    Failed to unzip tar file"
                        rm meme-5.5.2.tar.gz

                        # 进入 samtools 目录并编译
                        cd $tool_dir && mkdir -p build
                        ./configure                \
                            --prefix=$(pwd)/build  \
                            --enable-build-libxml2 \
                            --enable-build-libxslt > /dev/null 2>&1 || echo "    Failed to install configure"
                        make                       > /dev/null 2>&1 || print_red "    Failed to make"
                        make install               > /dev/null 2>&1 || print_red "    Failed to make install"

                        cd ..
                        mv $tool_dir ~/tools/

                        echo "export PATH=~/tools/$tool_dir/build/bin:\$PATH"                >> ~/.bashrc
                        echo "export PATH=~/tools/$tool_dir/build/bin:\$PATH"                >> ~/.zshrc
                        echo "export PATH=~/tools/$tool_dir/build/bin:\$PATH"                >> ~/.bash_profile
                        echo "export PATH=~/tools/$tool_dir/build/libexec/meme-5.5.2:\$PATH" >> ~/.bashrc
                        echo "export PATH=~/tools/$tool_dir/build/libexec/meme-5.5.2:\$PATH" >> ~/.zshrc
                        echo "export PATH=~/tools/$tool_dir/build/libexec/meme-5.5.2:\$PATH" >> ~/.bash_profile
                        sudo chown -R shiny:shiny-apps ~/tools/$tool_dir/build/
                        print_green "    MEME Suite successfully installed"
                    else
                        print_red "    Please install MEME suite later!"
                    fi # if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then
                fi # # if [ "$tool" == "fimo" ]; then
            done
        fi
    fi
fi

################################## 8. install R packages #################################
if [ true ]; then
    # 检查 R 是否已安装
    if command -v R >/dev/null 2>&1; then
        print_green_no_br "8. Would you like to install R packages? [y/N]: "
        # read -p "Would you like to install R packages? [y/N]: " answer
        read -p "" answer
        answer=${answer:-N} # Default to 'N' if no input provided

        if [ "$answer" == "Y" ] || [ "$answer" == "y" ]; then
            # 使用 R 脚本检查 rJava 是否已安装
            if ! Rscript -e "if (!requireNamespace('rJava', quietly = TRUE)) {quit(status = 1)}"; then
                # echo "rJava is not installed. Installing..."
                # 重新配置 Java 环境
                sudo R CMD javareconf                                         > /dev/null 2>&1
                curl -LO https://rforge.net/rJava/snapshot/rJava_1.0-6.tar.gz > /dev/null 2>&1
                tar fxz rJava_1.0-6.tar.gz                                    > /dev/null 2>&1
                R CMD INSTALL rJava                                           > /dev/null 2>&1 || print_red "    Failed to install rJava"
                rm rJava_1.0-6.tar.gz
                rm -rf rJava
            fi

            print_orange "Installing R packages... It takes minutes..."
            chmod a+x scripts/R_utils/install_packages.R
            # Rscript R/utils/install_packages.R
            Rscript scripts/R_utils/install_packages.R 2>&1 | tee log_R_packages_installation.log | grep -E "\* DONE \("
            awk '/The installed packages are as follows:/{flag=1} flag' log_R_packages_installation.log
        else
            print_orange "    No R packages installed"
        fi
    else
        print_red                "\n7. No R packages installed because there is no R installed."
        print_fluorescent_yellow "    Run 'scripts/R_utils/install_packages.R' to install R packages later. "
    fi # if command -v R >/dev/null 2>&1; then
fi

print_green "\nDONE"
