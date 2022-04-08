#!/usr/bin/env python
"""
@author: czn
"""
import os
import argparse
import time
import sys
parser = argparse.ArgumentParser()
parser.add_argument('--len', help='length threshold of sRNA (default: 15,34)', default='15,34')
parser.add_argument('--pvalue', help='pvalue threshold based on the Poisson distribution (default: 0.05)', default=0.05, type=float)
parser.add_argument('--extendregion', help='extending the upstream and downstream 10 nt of the sRNA is considered as '
                                           'the background region (default: 10)', default=10, type=int)
parser.add_argument('--windows_size', help='windows size (default: 11)', default=11, type=int)
parser.add_argument('--peakheight', help='minimum peak height (default: 5)', default=5, type=int)

group = parser.add_mutually_exclusive_group()
group.add_argument('--rpm', help='RPM threshold of boundry sites (default: 5)', type=int)
group.add_argument('--count', help='count threshold of boundry sites', type=int)

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--data', help='output path of Preprocessing', required=True)
requiredNamed.add_argument('--genomecovfile', help='genome coverage file', required=True)
requiredNamed.add_argument('--chromosome', help='viral accession or chromosome, eg:NC_009334.1 or chr1,chr2', required=True)
requiredNamed.add_argument('--outdir', help='output path', required=True)
requiredNamed.add_argument('--strand', help='strand of sRNA, eg: positive or negative', required=True)
requiredNamed.add_argument('--readfile', help='sequence alignments in BED file format', required=True)
parser.add_argument('--threads', help='running threads (default: 4)', default=4, type=int)

group2 = parser.add_mutually_exclusive_group()
group2.add_argument('--cleanread', help='the number of clean reads')
group2.add_argument('--mapInfor', help='bowtie mapInfor file')

args = parser.parse_args()
data = args.data
genomecovfile = args.genomecovfile.split('.txt')[0]
chromosome = args.chromosome.split(',')
windows_size = str(args.windows_size)
peakheight = str(args.peakheight)
outdir_ = args.outdir
strand = args.strand
readfile = args.readfile

softdir = sys.path[0]
peakInforfile = 'peak.y_smooth.'+strand+'.Infor.txt'
bedfile = 'sRNA.'+strand+'.bed'
outfile = 'sRNA.'+strand+'.txt'
ouputPoisson = 'sRNA.'+strand+'.poisson.txt'

minlen, maxlen = args.len.split(',')[0], args.len.split(',')[1]
rpm = str(args.rpm)
count = str(args.count)
mapInfor = args.mapInfor
pvalue = str(args.pvalue)
extendregion = str(args.extendregion)
threads = str(args.threads)
if args.rpm is None and args.count is None:
    rpm = '5'
if args.cleanread is not None:
    cleanread = args.cleanread
elif mapInfor is not None:
    file = open(data + '/' + mapInfor)
    cleanread = file.read().split('# reads processed: ')[1].split('\n')[0]
    file.close()
else:
    print('Please provide the number of clean reads or bowtie mapInfor')

for chr in chromosome:
    outdir = outdir_+'/'+chr
    outdir_tmp = outdir + '/_tmp'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        os.mkdir(outdir_tmp)
    os.system('python ' + softdir + '/Smooth.py --infile ' + data + '/' + genomecovfile + '.txt' + ' --outdir ' + outdir_tmp +
              ' --windows_size ' + windows_size + ' --chr ' + chr)
    if os.path.exists(outdir_tmp + '/' + genomecovfile + '.extend.txt.smooth.windows_size.' + windows_size + '.txt'):
        os.system('Rscript ' + softdir + '/FindPeak.R ' + genomecovfile + '.extend.txt.smooth.windows_size.' +
                  windows_size + '.txt ' + outdir_tmp + ' ' + strand + ' ' + peakheight)
        os.system('python ' + softdir + '/FindPeak.py --strand ' + strand +
                  ' --smoothfile ' + outdir_tmp + '/' + genomecovfile + '.extend.txt.smooth.windows_size.' + windows_size + '.txt'
                  ' --peakfile ' + outdir_tmp + '/peak.y_smooth.'+strand+'.txt'
                  ' --peakInforfile ' + outdir_tmp+ '/' + peakInforfile +
                  ' --readfile ' + data + '/' + readfile + ' --bed ' + outdir_tmp + '/' + bedfile + ' --chr ' + chr +
                  ' --outfile ' + outdir_tmp + '/' + outfile + ' --threads ' + threads)
        linecount = len(open(outdir_tmp + '/peak.y_smooth.'+strand+'.txt').readlines())
        if linecount != 0:
            if args.len == '15,34':
                if args.count is not None:
                    os.system('python ' + softdir + '/Poisson.py --method count --count ' + count + ' --cleanread ' + cleanread +
                              ' --genomecovfile ' + outdir_tmp + '/' + genomecovfile + '.extend.txt' +
                              ' --sRNAfile ' + outdir_tmp + '/' + outfile +
                              ' --pvalue ' + pvalue + ' --chr ' + chr + ' --extendregion ' + extendregion +
                              ' --output ' + ouputPoisson + ' --threads ' + threads + ' --outdir ' + outdir)
                else:
                    os.system('python ' + softdir + '/Poisson.py --method rpm --rpm ' + rpm + ' --cleanread ' + cleanread +
                              ' --genomecovfile ' + outdir_tmp + '/' + genomecovfile + '.extend.txt'
                              ' --sRNAfile ' + outdir_tmp + '/' + outfile +
                              ' --pvalue ' + pvalue + ' --chr ' + chr + ' --extendregion ' + extendregion +
                              ' --output ' + ouputPoisson + ' --threads ' + threads + ' --outdir ' + outdir)
            else:
                if args.count is not None:
                    os.system('python ' + softdir + '/Poisson.py --method count --count ' + count + ' --cleanread ' + cleanread +
                              ' --genomecovfile ' + outdir_tmp + '/' + genomecovfile + '.extend.txt' +
                              ' --sRNAfile ' + outdir_tmp + '/' + outfile +
                              ' --pvalue ' + pvalue + ' --chr ' + chr + ' --extendregion ' + extendregion +
                              ' --output ' + ouputPoisson + ' --minlen ' + minlen + ' --maxlen ' + maxlen +
                              ' --threads ' + threads + ' --outdir ' + outdir)
                else:
                    os.system('python ' + softdir + '/Poisson.py --method rpm --rpm ' + rpm + ' --cleanread ' + cleanread +
                              ' --genomecovfile ' + outdir_tmp + '/' + genomecovfile + '.extend.txt'
                              ' --sRNAfile ' + outdir_tmp + '/' + outfile +
                              ' --pvalue ' + pvalue + ' --chr ' + chr + ' --extendregion ' + extendregion +
                              ' --output ' + ouputPoisson + ' --minlen ' + minlen + ' --maxlen ' + maxlen +
                              ' --threads ' + threads + ' --outdir ' + outdir)

print('Finish time:', time.asctime(time.localtime(time.time())))
