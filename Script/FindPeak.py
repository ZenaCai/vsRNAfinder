#!/usr/bin/env python
"""
@author: czn
"""
import time
import pandas as pd
import linecache
from collections import Counter
import argparse
print('Start identify small RNAs:', time.asctime(time.localtime(time.time())))
parser = argparse.ArgumentParser()
parser.add_argument('--strand', help='positive or negative strand')
parser.add_argument('--smoothfile', help='genome coverage smooth file')
parser.add_argument('--peakfile', help='peak file')
parser.add_argument('--peakInforfile', help='output file of peak')
parser.add_argument('--readfile', help='sequence alignments in BED file format')
parser.add_argument('--outfile', help='output file of sRNA')
parser.add_argument('--bed', help='sRNA output in BED file format')
parser.add_argument('--chr', help='specific chromosome', type=str)
parser.add_argument('--threads', help='running threads', type=int)
args = parser.parse_args()

strand = args.strand
smoothfile = args.smoothfile
peakfile = args.peakfile
peakInforfile = args.peakInforfile
readfile = args.readfile
outfile = args.outfile
bedfile = args.bed
chr = args.chr
threads = args.threads

smooth = 'y_smooth'
df = pd.read_csv(smoothfile, header=0, sep='\t')
#####################################################################################
count = len(open(peakfile).readlines())
if count != 0:
    df_peak = pd.read_csv(peakfile, header=None, sep='\t')
    out = open(peakInforfile, 'w')
    linenum = 0
    out.write('Genome	Peak site	Peak start	Peak end	Peak abundance	Peak length	Peak gap\n')
    gap = 0
    with open(peakfile) as file:
        for line in file:
            linenum += 1
            line1 = line.strip().split('\t')
            index = []
            for i in range(int(line1[2]) - 1, int(line1[3])):
                index.append(i)
            peak_num = df.loc[index, ][smooth].tolist().count(float(line1[0]))
            if linenum != count:
                nexline = linecache.getline(peakfile, linenum + 1)
                gap = int(nexline.strip().split('\t')[5]) - int(line1[5])
                out.write(line1[4] + '\t' + line1[5] + '\t' + line1[8] + '\t' + line1[9] + '\t' + line1[7] + '\t' +
                          str(peak_num) + '\t' + str(gap) + '\n')
            else:
                out.write(line1[4] + '\t' + line1[5] + '\t' + line1[8] + '\t' + line1[9] + '\t' + line1[7] + '\t' +
                          str(peak_num) + '\t' + str(gap) + '\n')
    out.close()

    """Small RNA region"""
    def GetBoundary(start, end, add, cut, dict):
        site_lst = df['1'].tolist()
        list = []
        dict_site = {}
        site, num = 0, 0
        dict_site[num] = ['Nan']
        for i in range(start, end, add):
            if i in site_lst:
                abundance = float(df[df['1'] == i][smooth])
                if abundance >= cut:
                    list.append(i)
                    if i in dict:
                        value = dict[i]
                        if value > num:
                            site = i
                            num = value
                        if value not in dict_site:
                            dict_site[value] = [i]
                        else:
                            dict_site[value].append(i)
        site_ = dict_site[num]
        return site_, num

    df = pd.read_csv(smoothfile, header=0, sep='\t')
    dict_chr_start_positive = {}
    dict_chr_end_positive = {}
    dict_chr_start_negative = {}
    dict_chr_end_negative = {}
    with open(readfile) as file:
        for line in file:
            line1 = line.strip().split('\t')
            if line1[-1] == '+':
                if line1[0] not in dict_chr_start_positive:
                    dict_chr_start_positive[line1[0]] = [int(line1[1]) + 1]
                    dict_chr_end_positive[line1[0]] = [int(line1[2])]
                else:
                    dict_chr_start_positive[line1[0]].append(int(line1[1]) + 1)
                    dict_chr_end_positive[line1[0]].append(int(line1[2]))
            else:
                if line1[0] not in dict_chr_start_negative:
                    dict_chr_start_negative[line1[0]] = [int(line1[1]) + 1]
                    dict_chr_end_negative[line1[0]] = [int(line1[2])]
                else:
                    dict_chr_start_negative[line1[0]].append(int(line1[1]) + 1)
                    dict_chr_end_negative[line1[0]].append(int(line1[2]))
    positive_startlst, positive_endlst = [], []
    negative_startlst, negative_endlst = [], []
    if chr in dict_chr_start_positive:
        positive_startlst = dict_chr_start_positive[chr]
    if chr in dict_chr_end_positive:
        positive_endlst = dict_chr_end_positive[chr]
    if chr in dict_chr_start_negative:
        negative_startlst = dict_chr_start_negative[chr]
    if chr in dict_chr_end_negative:
        negative_endlst = dict_chr_end_negative[chr]
    def GetsRNA(line):
        output = 'Nan'
        percenta = 0
        line1 = line.strip().split('\t')
        peak_site, peak_abundance = int(line1[1]), float(line1[4])
        cut = peak_abundance * percenta
        if strand == 'positive':
            start_positive, end_positive = GetBoundary(peak_site, peak_site - 30, -1, cut,
                                                       Counter(positive_startlst)), \
                                           GetBoundary(peak_site + 1, peak_site + 30 + 1, 1, cut,
                                                       Counter(positive_endlst))
            if start_positive[0][0] != 'Nan' and end_positive[0][0] != 'Nan':
                output = (str(peak_site) + '\t' + str(min(start_positive[0])) + '\t' + str(max(end_positive[0])) +
                          '\t' + str(max(end_positive[0]) - min(start_positive[0]) + 1) + '\tpositive\t' +
                          str(peak_site) + '-' + str(peak_site - 30) + '\t' + str(peak_site + 1) + '-' +
                          str(peak_site + 30 + 1) + '\t' + str(start_positive) + '\t' + str(end_positive) + '\t' +
                          str(start_positive[1]) + ',' + str(end_positive[1]) + '\n')
            else:
                output = (str(peak_site) + '\tNan\tNan\tNan\tpositive\t' +
                          str(peak_site) + '-' + str(peak_site - 30) + '\t' + str(peak_site + 1) + '-' +
                          str(peak_site + 30 + 1) + '\t' + str(start_positive) + '\t' + str(end_positive) + '\tNan\n')
        elif strand == 'negative':
            start_negative, end_negative = GetBoundary(peak_site, peak_site - 30, -1, cut,
                                                       Counter(negative_startlst)), \
                                           GetBoundary(peak_site + 1, peak_site + 30 + 1, 1, cut,
                                                       Counter(negative_endlst))
            if start_negative[0][0] != 'Nan' and end_negative[0][0] != 'Nan':
                output = (str(peak_site) + '\t' + str(min(start_negative[0])) + '\t' + str(max(end_negative[0])) +
                          '\t' + str(max(end_negative[0]) - min(start_negative[0]) + 1) + '\tnegative\t' +
                          str(peak_site) + '-' + str(peak_site - 30) + '\t' + str(peak_site + 1) + '-' +
                          str(peak_site + 30 + 1) + '\t' + str(start_negative) + '\t' + str(end_negative) + '\t' +
                          str(start_negative[1]) + ',' + str(end_negative[1]) + '\n')
            else:
                output = (str(peak_site) + '\tNan\tNan\tNan\tnegative\t' +
                          str(peak_site) + '-' + str(peak_site - 30) + '\t' + str(peak_site + 1) + '-' +
                          str(peak_site + 30 + 1) + '\t' + str(start_negative) + '\t' + str(end_negative) + '\tNan\n')
        if output != 'Nan':
            return output

    from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
    p = ProcessPoolExecutor(threads)
    linenum = 0
    output_result = []
    with open(peakInforfile) as file:
        for line in file:
            linenum += 1
            if linenum != 1:
                output_result.append(p.submit(GetsRNA, line))
    p.shutdown()
    out = open(outfile, 'w')
    out.write('Peak site\tStart_min\tEnd_max\tLength\tStrand\tStartRange\tEndRange\tStart_lst\tEnd_lst\tcount\n')
    for i in output_result:
        out.write(i.result())
    out.close()

    df = pd.read_csv(outfile, sep='\t', header=0, index_col=None)
    df = df[(df['Start_min'] != 'Nan') & (df['End_max'] != 'Nan')]
    df.insert(0, 'site', (df['Start_min'].astype(str) + '-' + df['End_max'].astype(str) + '-' +
                          df['Strand']).tolist())
    length_lst = df['Length'].astype(int).tolist()
    site_lst = set(df['site'].tolist())
    print('small RNA region', len(site_lst))
    if len(site_lst) != 0:
        print('min length & max length', min(length_lst), max(length_lst))
        from itertools import groupby
        num = 5
        for k, g in groupby(sorted(length_lst), key=lambda x: x // num):
            print('{}-{}: {}'.format(k * num, (k + 1) * num - 1, len(list(g))))

        """ BED Format """
        df_bed = pd.DataFrame()
        df_bed['Genome'] = [chr]*df.shape[0]
        df_bed['Start'], df_bed['End'] = (df['Start_min'].astype(int) - 1).tolist(), \
                                         df['End_max'].astype(int).tolist()
        df_bed['Name'], df_bed['Score'] = df['site'].tolist(), ['.']*df.shape[0]
        if strand == 'positive':
            df_bed['Strand'] = df['Strand'].replace('positive', '+').tolist()
        elif strand == 'negative':
            df_bed['Strand'] = df['Strand'].replace('negative', '-').tolist()
        df_bed.to_csv(bedfile, sep='\t', index=None, header=None)
    else:
        print(chr, 'no sRNA')
else:
    print(chr, 'no sRNA')
