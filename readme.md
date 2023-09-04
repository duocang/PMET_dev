This is PMET.

## File tree

```shell
.
├── 01_homotypic_intervals.sh
├── 01_homotypic_promoters.sh
├── 02_heterotypic_intervals.sh
├── 02_heterotypic_promoters.sh
├── 02_heterotypic_promoters_single_CPU.sh
├── 03_test_new_fimo.sh
├── data
├── results
├── scripts
├── src
├── visualize_pmet_php
└── readme.md
```

![](https://raw.githubusercontent.com/duocang/images/master/PicGo/202307202339573.png)

## Compile PMET homotypic and PMET heterotypic

Both are writen in C++, source code can be found in `src/indexing` and `src/pmetParallel`.

If necessary, it is possible to compile `PMET index` and pmet in different OS.

**Compile in one bash**

```bash
chmod a+x 00_binary_compile.sh

bash 00_binary_compile.sh
```



After running the bash, all needed binary tools will be put in the `scripts` folder.

## 1. TEST PMET

```bash
chmod a+x 01_homotypic_promoters.sh
chmod a+x 02_heterotypic_promoters.sh
```

### 1.1 search and filter homotypic motifs matches in all promoters

```bash
bash 01_homotypic_promoters.sh
# This can take a long time.
```

### 1.2 search heterotypic motifs matching in all promoters

```bash
bash 02_heterotypic_promoters.sh
```


## 2. TEST new FIMO

Before running `PMET index`, we need to run FIMO to find all the homotypic motifs, and then `PMET index` will use the results from FIMO to run. This process will consume **IO** resources.

**IO is expensive.**

To mitigate the IO resource consumption associated with FIMO and the `PMET index`, we aim to integrate the capabilities of the `PMET index` directly into FIMO. Details of this integration can be found in the `src/meme-5.3.3`directory.

For instance, when querying 113 motif hits on the promoter of the Arabidopsis thaliana genome, the improved FIMO (referred to as NEW FIMO) can reduce write operations by 30GB and read operations by the same amount.

```bash
chmod a+x 03_test_new_fimo.sh
bash 03_test_new_fimo.sh
```

Using the hardware specifications listed below, the traditional combination of FIMO and the `PMET index` takes more than double the time compared to using NEW FIMO. On a mechanical hard drive, this time difference can be amplified, possibly reaching 5 to 10 times.

- Single core processor
- Intel i9-12900K
- Samsung 980 Pro SSD

An additional consideration is the potential to divide the meme files (motifs) into segments and employ `GNU Parallel` for concurrent processing. This approach would decrease run times. Moreover, it would amplify the efficiency of  `NEW FIMO`. Given that IO resources are finite, the IO resource usage of the combined `FIMO` and `PMET index` increases multiplicatively. **No less time with more threads**.

In contrast, `NEW FIMO` circumvents this issue entirely.


## Install GNU Parallel

GNU Parallel helps `PMET index` (FIMO and `PMET index`) to run in parallel mode.

```bash
sudo apt-get update
sudo apt-get install parallel
```

Put GNU Parallel silent:

```bash
 # Run once
 parallel --citation
```

## Install zentiy

```bash
sudo apt-get update
sudo apt-get install zenity
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