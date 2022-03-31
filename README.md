# **vsRNAfinder**

vsRNAfinder is a novel method for identifying viral small RNAs from small RNA-Seq data

## Dependencies

### Packages 

1. bowtie 1.3.0
2. samtools 1.9
3. bedtools v2.29.2
4. RNAfold 2.4.18

### R Dependencies

Base on R 3.6.1

pracma 2.3.3

### Python Dependencies

Base on python 3.6.8

1. Numpy 1.19.2
2. Pandas 1.1.5
3. Scipy 1.5.2
4. Pysam 0.15.3

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

### How to run the test dataset using vsRNAfinder

```shell
# Step1
python Preprocessing.py --cleanfq test_data/test.clean.fastq --genome test_data/NC_009334.1.fasta --outdir test_output --prefix test --threads 10
# Step2
python FindSmallRNA.py --data test_output --genomecovfile test.sort.positive.bga.txt --chromosome NC_009334.1 --outdir test_output --strand positive --readfile test.sort.bed --mapInfor test.mapInfo.txt --threads 10
python FindSmallRNA.py --data test_output --genomecovfile test.sort.negative.bga.txt --chromosome NC_009334.1 --outdir test_output --strand negative --readfile test.sort.bed --mapInfor test.mapInfo.txt --threads 10
# Step3
python FindMiRNA.py --genome test_data/NC_009334.1.fasta --speices virus --data test_output/NC_009334.1
# Step4
python Quantification.py --threads 10 --data test_output/NC_009334.1 --bam test_output/test.sort.bam --mapInfor test_output/test.mapInfo.txt

```





