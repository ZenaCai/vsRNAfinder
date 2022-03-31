#!/usr/bin/env python
"""
@author: czn
"""
import pandas as pd
import argparse
import os
import re
import pysam
parser = argparse.ArgumentParser()
parser.add_argument('--threads', help='running threads (default: 4)', default=4, type=int)
parser.add_argument('--threshold', help='minimum fraction of overlapping bases in a read that '
                                        'is required for read assignment (default: 0.8)', default=0.8)
parser.add_argument('--format', help='SAM or BAM format', default='BAM')
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--data', help='output path of FindSmallRNA', required=True)
requiredNamed.add_argument('--bam', help='bam file', required=True)
group = parser.add_mutually_exclusive_group()
group.add_argument('--cleanread', help='the number of clean reads')
group.add_argument('--mapInfor', help='bowtie mapInfor file')
args = parser.parse_args()

threads = args.threads
data = args.data
bam = args.bam
mapInfor = args.mapInfor
threshold = args.threshold
format = args.format

if args.cleanread is not None:
    cleanread = args.cleanread
elif mapInfor is not None:
    file = open(mapInfor)
    cleanread = file.read().split('# reads processed: ')[1].split('\n')[0]
    file.close()
else:
    print('Please provide the number of clean reads or bowtie mapInfor')

if format == 'BAM':
    mode = "rb"
elif format == 'SAM':
    mode = "r"
else:
    print('Please Input File Format')
def GetResult(samfile, chr, start, end, strand, threshold, cleanread, ID):  # start-1
    count = 0
    for read in samfile.fetch(chr, start, end):  # bam: base-1
        read_start, read_end, read_length = read.reference_start, read.reference_end, read.qlen
        overlap = read.get_overlap(start, end)
        percentage = overlap / read_length
        if percentage >= threshold:
            if (strand == '+' and not read.is_reverse) or (strand == '-' and read.is_reverse):
                count += 1
        else:
            if read_start == start or read_end == end:
                if (strand == '+' and not read.is_reverse) or (strand == '-' and read.is_reverse):
                    count += 1
    length = end - start
    output = ID + '\t' + chr + '\t' + str(start + 1) + '\t' + str(end) + '\t' + strand + '\t' + str(length) + '\t' + str(
        count) + '\t' + str(count * 1000000 / int(cleanread))
    return output

out = open(data + '/sRNA.counts.txt', 'w')
out.write("sRNA\tChr\tStart\tEnd\tStrand\tLength\tCounts\tRPM\n")
samfile = pysam.AlignmentFile(bam, mode, threads=threads)
for filename in os.listdir(data):
    if re.search('.saf', filename):
        linecount = len(open(data + '/' + filename).readlines())
        if linecount > 1:
            df_saf = pd.read_csv(data + '/' + filename, header=0, sep="\t")
            for row in df_saf.itertuples():
                chr, start, end, strand = str(getattr(row, 'Chr')), getattr(row, 'Start'), getattr(row, 'End'), getattr(row, 'Strand')
                ID = getattr(row, 'GeneID')
                output = GetResult(samfile, chr, start-1, end, strand, threshold, cleanread, ID)
                out.write(output+'\n')
        else:
            print('No data for Quantification')
out.close()
