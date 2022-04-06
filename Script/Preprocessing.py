#!/usr/bin/env python
"""
@author: czn
"""
import os
import argparse
import time
parser = argparse.ArgumentParser()
parser.add_argument('--threads', help='running threads (default: 4)', default=4)
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--cleanfq', help='fastq file after trimming adapter', required=True)
requiredNamed.add_argument('--genome', help='genome fasta file', required=True)
requiredNamed.add_argument('--outdir', help='output path', required=True)
requiredNamed.add_argument('--prefix', help='the prefix of output file', required=True)
args = parser.parse_args()

threads = args.threads
cleanfq = args.cleanfq
genome = args.genome
outdir = args.outdir
prefix = args.prefix
if not os.path.exists(outdir):
    os.mkdir(outdir)
bowtie_index = outdir+'/'+prefix+'.bowtie_index.log'
sam = outdir+'/'+prefix+'.sam'
map_fastq = outdir+'/'+prefix+'.map.fastq'
mapInfo = outdir+'/'+prefix+'.mapInfo.txt'
unmap_fastq = outdir+'/'+prefix+'.unmap.fastq'
os.system('bowtie-build --threads '+str(threads)+' '+genome+' '+genome+'.bowtie'+' > ' + bowtie_index + ' 2>&1')
os.system('bowtie --threads '+str(threads)+' -v 0 -m 2 -a '+genome+'.bowtie '+cleanfq+' --un '+unmap_fastq+
          ' --al '+map_fastq+' -S '+sam+' > '+mapInfo + ' 2>&1')

bam = outdir+'/'+prefix+'.bam'
bam_sort = outdir+'/'+prefix+'.sort.bam'
bed = outdir+'/'+prefix+'.sort.bed'
os.system('samtools view -@ '+str(threads)+' -b -o '+bam+' '+sam)
os.system('samtools sort -@ '+str(threads)+' -o '+bam_sort+' '+bam)
os.system('samtools index '+bam_sort)
os.system('bedtools bamtobed -i '+bam_sort + ' > '+bed)
positive_depth = outdir+'/'+prefix+'.sort.positive.bga.txt'
negative_depth = outdir+'/'+prefix+'.sort.negative.bga.txt'
os.system('bedtools genomecov -ibam '+bam_sort+' -bga -split -strand + > '+positive_depth)
os.system('bedtools genomecov -ibam '+bam_sort+' -bga -split -strand - > '+negative_depth)
print('Preprocessing finish time:', time.asctime(time.localtime(time.time())))
