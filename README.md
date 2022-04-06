# **vsRNAfinder**

vsRNAfinder is a novel method for identifying viral small RNAs from small RNA-Seq data

## Dependencies

### Packages 

1. bowtie 1.3.0
2. samtools 1.9
3. bedtools 2.29.2
4. RNAfold 2.4.18 ( Viennarna 2.4.18 )

### R Dependencies

Base on R 3.6.1

pracma 2.3.3

### Python Dependencies

Base on python 3.6.8

1. Numpy 1.19.2
2. Pandas 1.1.5
3. Scipy 1.5.2
4. Pysam 0.15.3

## Download & Installation

Download vsRNAfinder with Git from GitHub

~~~
git clone https://github.com/ZenaCai/vsRNAfinder.git
~~~

Create a new conda environment

~~~
conda create -n vsRNAfinder python=3.6.8
conda activate vsRNAfinder
~~~

Install Packages

~~~
conda install -c bioconda bowtie=1.3.0
conda install -c bioconda samtools=1.9
conda install -c bioconda bedtools=2.29.2
conda install -c bioconda viennarna=2.4.18
~~~

Install R and Python dependencies

~~~
conda install r-base=3.6.1
conda install -c conda-forge r-pracma=2.3.3
conda install numpy=1.19.2
conda install pandas=1.1.5
conda install scipy=1.5.2
conda install -c bioconda pysam=0.15.3
~~~

## Usage

### Running vsRNAfinder pipeline

#### Step1: Preprocessing

Input fastq file should be trimmed for adapters and  too short reads (≤15 nt) should be discarded by cutadapt before running vsRNAfinder. 

~~~shell
usage: Preprocessing.py [-h] [--threads THREADS] --cleanfq CLEANFQ --genome
                        GENOME --outdir OUTDIR --prefix PREFIX

optional arguments:
  -h, --help         show this help message and exit
  --threads THREADS  running threads (default: 4)

required arguments:
  --cleanfq CLEANFQ  fastq file after trimming adapter
  --genome GENOME    genome fasta file
  --outdir OUTDIR    output path
  --prefix PREFIX    the prefix of output file

~~~

#### Step2: sRNA Identification and filtering sRNA by abundance, length,  significance of sRNA

~~~shell
usage: FindSmallRNA.py [-h] [--len LEN] [--pvalue PVALUE]
                       [--extendregion EXTENDREGION]
                       [--windows_size WINDOWS_SIZE] [--peakheight PEAKHEIGHT]
                       [--rpm RPM | --count COUNT] --data DATA --genomecovfile
                       GENOMECOVFILE --chromosome CHROMOSOME --outdir OUTDIR
                       --strand STRAND --readfile READFILE [--threads THREADS]
                       [--cleanread CLEANREAD | --mapInfor MAPINFOR]

optional arguments:
  -h, --help            show this help message and exit
  --len LEN             length threshold of sRNA (default: 15,34)
  --pvalue PVALUE       pvalue threshold based on the Poisson distribution
                        (default: 0.05)
  --extendregion EXTENDREGION
                        extending the upstream and downstream 10 nt of the
                        sRNA is considered as the background region (default:
                        10)
  --windows_size WINDOWS_SIZE
                        windows size (default: 11)
  --peakheight PEAKHEIGHT
                        minimum peak height (default: 5)
  --rpm RPM             RPM threshold of boundry sites (default: 5)
  --count COUNT         count threshold of boundry sites
  --threads THREADS     running threads (default: 4)
  --cleanread CLEANREAD
                        the number of clean reads
  --mapInfor MAPINFOR   bowtie mapInfor file

required arguments:
  --data DATA           output path of Preprocessing
  --genomecovfile GENOMECOVFILE
                        genome coverage file
  --chromosome CHROMOSOME
                        viral accession or chromosome, eg:NC_009334.1 or
                        chr1,chr2
  --outdir OUTDIR       output path
  --strand STRAND       strand of sRNA, eg: positive or negative
  --readfile READFILE   sequence alignments in BED file format

~~~

#### Step3: miRNA Identification for each chromosome (optional step)

~~~shell
usage: FindMiRNA.py [-h] [--rpm RPM | --count COUNT] --genome GENOME --speices
                    SPEICES --data DATA

optional arguments:
  -h, --help         show this help message and exit
  --rpm RPM          RPM threshold of boundry sites
  --count COUNT      count threshold of boundry sites

required arguments:
  --genome GENOME    genome fasta file
  --speices SPEICES  virus, animal or plant (default: virus)
  --data DATA        output path of FindSmallRNA
~~~

#### Step4: Quantification for each chromosome

~~~shell
usage: Quantification.py [-h] [--threads THREADS] [--threshold THRESHOLD]
                         [--format FORMAT] --data DATA --bam BAM
                         [--cleanread CLEANREAD | --mapInfor MAPINFOR]

optional arguments:
  -h, --help            show this help message and exit
  --threads THREADS     running threads (default: 4)
  --threshold THRESHOLD
                        minimum fraction of overlapping bases in a read that
                        is required for read assignment (default: 0.8)
  --format FORMAT       SAM or BAM format
  --cleanread CLEANREAD
                        the number of clean reads
  --mapInfor MAPINFOR   bowtie mapInfor file

required arguments:
  --data DATA           output path of FindSmallRNA
  --bam BAM             bam file
~~~

#### Output

- **outdir/chromosome/sRNA.xxx.poisson.txt**

  Result for small RNAs

- **outdir/chromosome/Result.txt**

  Result for miRNAs and other small RNAs

- **outdir/chromosome/sRNA.counts.txt**

  Quantification for small RNAs

Each column is described as follows:

| *Column*    | *Description*                                                |
| ----------- | ------------------------------------------------------------ |
| Site/sRNA   | Site of the sRNA in format Start-End-Strand-Chr              |
| Chr         | Chromosome of the sRNA                                       |
| Start       | Start position of the sRNA                                   |
| End         | End position of the sRNA                                     |
| Strand      | Strand of the sRNA                                           |
| Length      | Length of the sRNA                                           |
| Start_count | Number of reads starting at the start position of  the sRNA  |
| End_count   | Number of reads ending at the end position of  the sRNA      |
| Start_rpm   | The abundance of the start position of the sRNA is normalized using RPM (Reads Per Million) |
| End_rpm     | The abundance of the  end position of the sRNA is normalized using RPM (Reads Per Million) |
| Type        | Type of sRNA (miRNA or sRNA)                                 |
| Sequence    | Sequence of the sRNA                                         |
| Count       | Number of reads mapping to the sRNA                          |
| RPM         | The abundance of the sRNA is normalized using RPM (Reads Per Million) |

### How to run the test dataset using vsRNAfinder

```shell
# Step1
python Script/Preprocessing.py --cleanfq test_data/test.clean.fastq --genome test_data/NC_009334.1.fasta --outdir test_output --prefix test --threads 10
# Step2
python Script/FindSmallRNA.py --data test_output --genomecovfile test.sort.positive.bga.txt --chromosome NC_009334.1 --outdir test_output --strand positive --readfile test.sort.bed --mapInfor test.mapInfo.txt --threads 10
python Script/FindSmallRNA.py --data test_output --genomecovfile test.sort.negative.bga.txt --chromosome NC_009334.1 --outdir test_output --strand negative --readfile test.sort.bed --mapInfor test.mapInfo.txt --threads 10
# Step3
python Script/FindMiRNA.py --genome test_data/NC_009334.1.fasta --speices virus --data test_output/NC_009334.1
# Step4
python Script/Quantification.py --threads 10 --data test_output/NC_009334.1 --bam test_output/test.sort.bam --mapInfor test_output/test.mapInfo.txt

```





