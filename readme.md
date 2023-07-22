This is a Shiny app developed for PMET.



## File tree

```shell
.
├── 01_CPP_debug_pmetindex_intervals.sh
├── 01_CPP_debug_pmetindex_promoters.sh
├── 02_CPP_debug_pmet_intervals_parallel.sh
├── 02_CPP_debug_pmet_promoters.sh
├── 02_CPP_debug_pmet_promoters_parallel.sh
├── data
├── performance
├── readme.md
├── results
├── scripts
├── src
└── visualize_pmet_php
```

![](https://raw.githubusercontent.com/duocang/images/master/PicGo/202307202339573.png)

## PMET index and PMET

Both are writen in C++, source code can be found in `src/indexing` and `src/pmetParallel`.

If necessary, it is possible to compile pmet index and pmet in different OS.

**PMET index**

```bash
# src/indexing
./build.sh
```

or

```bash
# src/indexing
g++  -g -Wall -std=c++11 cFimoFile.cpp cMotifHit.cpp fastFileReader.cpp main.cpp -o ../../scripts/pmetindex
```

**PMET**

```bash
# src/pmetParallel
./build.sh
```

or

```bash
# src/pmetParallel
g++ -g -Wall -std=c++11 Output.cpp motif.cpp motifComparison.cpp main.cpp -o ../../scripts/pmetParallel_linux -pthread
```

## Install GNU Parallel

GNU Parallel helps PMET index (FIMO and PMET index) to run in parallel mode.

```bash
sudo apt-get install parallel
```

Put GNU Parallel silent:

```bash
 # Run once
 parallel --citation
```

## Install The MEME Suite (FIMO and fasta-get-markov)

```bash
# cd a folder you want to put the software
wget https://meme-suite.org/meme/meme-software/5.5.2/meme-5.5.2.tar.gz

tar zxf meme-5.5.2.tar.gz
cd meme-5.5.2
./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
make
make test
make install
```

Add following into bash profile file.

```bash
# assuming you put meme folder under your home folder
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.2:$PATH
```

## Install samtools

Install from conda or mamba:

```bash
conda install -c bioconda samtools
```

Install from source:

> assuming you create a directory named `samtools` in home directory (~) and install samtools there.

```bash
wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2

cd samtools-1.17    # and similarly for bcftools and htslib
./configure --prefix=$HOME/samtools
make
make install

# Add following into bash profile file or .zshrc (if zsh used).

# assuming you put samtools-1.17 folder under your home folder
export PATH=$HOME/samtools/bin:$PATH
```

## Install bedtools

It is recommended to install bedtools via apt/yum or conda.

```bash
conda install -c bioconda bedtools
```

or

```bash
conda install bedtools

# Debian/Ubuntu
apt-get install bedtools

# Fedora/Centos
yum install BEDTools
```

It is possible to compile the bedtools by running the following commands.

```bash
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
```

## Python libraries

```bash
pip install numpy
pip install pandas
pip install scipy
pip install bio
pip install biopython
```