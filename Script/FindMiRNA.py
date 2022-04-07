#!/usr/bin/env python
"""
@author: czn
"""
import pandas as pd
from collections import Counter
import time
import os
import re
import argparse
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
group.add_argument('--rpm', help='RPM threshold of boundry sites', type=float)
group.add_argument('--count', help='count threshold of boundry sites', type=int)

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--genome', help='genome fasta file', required=True)
requiredNamed.add_argument('--speices', help='virus, animal or plant (default: virus)', required=True, default='virus')
requiredNamed.add_argument('--data', help='output path of FindSmallRNA', required=True)
args = parser.parse_args()
genome = args.genome
speices = args.speices
FindPeak_dir = args.data
RNAfold_dir = FindPeak_dir + '/RNAfold'
rpm = args.rpm
count = args.count

percent = 0.7
energy_cut = -25
##################################################################################
def GetmiRNA(filename, smiRNA_lst, energy_cut, percent, dict_seq):
    df = pd.read_csv(filename + '.bed', header=None, sep='\t')
    df['key'] = df[0].astype(str)+':'+df[1].astype(str)+'-'+df[2].astype(str)+'('+df[5]+')'
    with open(filename + '.str') as file:
        for line in file:
            if re.match('>', line):
                ID = line.strip().split('>')[1]
                sRNA = df[df['key'] == ID][3].tolist()[0]
                if re.search('\+', ID):
                    sRNA_start = int(sRNA.split('-')[0]) - int(df[df['key'] == ID][1].tolist()[0]) - 1
                    sRNA_end = int(sRNA.split('-')[1]) - int(df[df['key'] == ID][1].tolist()[0])
                    left = int(sRNA.split('-')[0]) - int(df[df['key'] == ID][1].tolist()[0]) - 1
                    rigth = int(df[df['key'] == ID][2].tolist()[0]) - int(sRNA.split('-')[1])
                else:
                    sRNA_start = -(int(sRNA.split('-')[1]) - int(df[df['key'] == ID][1].tolist()[0]))
                    sRNA_end = -(int(sRNA.split('-')[0]) - int(df[df['key'] == ID][1].tolist()[0]) - 1)
                    left = int(df[df['key'] == ID][2].tolist()[0]) - int(sRNA.split('-')[1])
                    rigth = int(sRNA.split('-')[1]) - int(df[df['key'] == ID][1].tolist()[0])
            elif not re.match('>', line) and not re.search('\(', line):
                sequence = line.strip()
                if sRNA_end == 0:
                    sRNA_seq = sequence[sRNA_start:]
                else:
                    sRNA_seq = sequence[sRNA_start:sRNA_end]
                dict_seq[sRNA] = [sRNA_seq]
            elif not re.match('>', line) and re.search('\(', line):
                loop = re.compile('(\(+\.+\)+)').findall(line)
                sRNA_infor = '*' * left + line.split(' ')[0][sRNA_start:sRNA_end] + '*' * rigth
                sRNA_len = len(line.split(' ')[0][sRNA_start:sRNA_end])
                if sRNA_len != 0:
                    if re.search(' ', line.split(' (')[1].split(')')[0]):
                        energy = float(line.split(' (')[1].split(')')[0].split(' ')[-1])
                    else:
                        energy = float(line.split(' (')[1].split(')')[0])

                    if len(loop) == 1:
                        ## single loop
                        for sRNA in df[df['key'] == ID][3].tolist():
                            if not re.search('\)', sRNA_infor):
                                if Counter(sRNA_infor)['('] / sRNA_len >= percent and energy <= energy_cut:
                                    smiRNA_lst.append(sRNA)
                            elif not re.search('\(', sRNA_infor):
                                if Counter(sRNA_infor)[')'] / sRNA_len >= percent and energy <= energy_cut:
                                    smiRNA_lst.append(sRNA)
                    else:
                        ## multiple loops
                        for sRNA in df[df['key'] == ID][3].tolist():
                            if not re.search('\)', sRNA_infor):
                                if Counter(sRNA_infor)['('] / sRNA_len >= percent and energy <= energy_cut:
                                    smiRNA_lst.append(sRNA)
                            elif not re.search('\(', sRNA_infor):
                                if Counter(sRNA_infor)[')'] / sRNA_len >= percent and energy <= energy_cut:
                                    smiRNA_lst.append(sRNA)

##################################################################################
genomeIndex = genome + '.fai'
if not os.path.exists(genomeIndex):
    os.system('samtools faidx ' + genome)
dict_chr = {}
with open(genomeIndex) as file:
    for line in file:
        line1 = line.strip().split('\t')
        dict_chr[line1[0]] = int(line1[1])
#################################################################################
if speices == 'plant':
    step = 7
    ii = 10
else:## animal & virus
    step = 5
    ii = 8
loop = 15
df_all = pd.DataFrame()
df_positive = pd.DataFrame()
df_negative = pd.DataFrame()
for filename in os.listdir(FindPeak_dir):
    if re.search('poisson', filename) and not re.search('.saf', filename) and not re.search('counts', filename):
        df_sub = pd.read_csv(FindPeak_dir + '/' + filename, header=0, sep='\t')
        if rpm is not None:
            df_sub = df_sub[(df_sub['Start_rpm'] >= rpm) & (df_sub['End_rpm'] >= rpm)]
        if count is not None:
            df_sub = df_sub[(df_sub['Start_count'] >= rpm) & (df_sub['End_count'] >= rpm)]
        df_sub['Gap'] = ['.']*df_sub.shape[0]
        df_sub['Strand'] = df_sub['Strand'].replace('negative', '-')
        df_sub['Strand'] = df_sub['Strand'].replace('positive', '+')
        df_sub_positive = df_sub[df_sub['Strand'] == '+']
        df_sub_negative = df_sub[df_sub['Strand'] == '-']
        colname = ["Site", "Chr", "Start", "End", "Strand", "Length", "Start_count", "End_count", "Start_rpm", "End_rpm", "Pvalue"]
        df_all = pd.concat([df_all, df_sub_positive[colname], df_sub_negative[colname]])
        df_positive = pd.concat([df_positive, df_sub_positive])
        df_negative = pd.concat([df_negative, df_sub_negative])
if not df_all.empty:
    if not os.path.exists(RNAfold_dir):
        os.mkdir(RNAfold_dir)
    smiRNA_lst = []
    dict_seq = {}
    for tag in [5, 3]:
        df = pd.DataFrame()
        if tag == 5:
            ### 5'
            df_positive['Start_extend'] = df_positive['Start'] - step*ii
            df_positive['End_extend'] = df_positive['End'] + df_positive['Length'] + loop + step*ii
            df_negative['Start_extend'] = df_negative['Start'] - (df_negative['Length'] + loop + step*ii)
            df_negative['End_extend'] = df_negative['End'] + step*ii
        else:
            ### 3'
            df_positive['Start_extend'] = df_positive['Start'] - (df_positive['Length'] + loop + step*ii)
            df_positive['End_extend'] = df_positive['End'] + step*ii
            df_negative['Start_extend'] = df_negative['Start'] - step*ii
            df_negative['End_extend'] = df_negative['End'] + (df_negative['Length'] + loop + step*ii)
        colname = ['Chr', 'Start_extend', 'End_extend', 'Site', 'Gap', 'Strand']
        df = pd.concat([df, df_positive[colname], df_negative[colname]])
        df.loc[df['Start_extend'] < 0, 'Start_extend'] = 0
        for i in dict_chr:
            df.loc[(df['Chr'].astype(str) == i) & (df['End_extend'] > dict_chr[str(i)]), 'End_extend'] = dict_chr[str(i)]
        filename = RNAfold_dir + '/miRNA.Precursor.'+str(tag)
        df.to_csv(filename + '.bed', sep='\t', header=None, index=None)
        os.system('bedtools getfasta -fi ' + genome +' -bed ' + filename + '.bed' + ' -s > ' + filename + '.fasta')
        os.system('RNAfold --noPS ' + filename +'.fasta' + ' > ' + filename + '.str')
        GetmiRNA(filename, smiRNA_lst, energy_cut, percent, dict_seq)
    df_seq = pd.DataFrame(dict_seq)

    sRNA_all = df_all['Site'].tolist()
    print('small RNA\t' + str(len(set(sRNA_all))))
    print('miRNA\t' + str(len(set(smiRNA_lst))))
    print('other small RNA\t' + str(len(set(sRNA_all) - set(smiRNA_lst))))
    df_sRNA = pd.DataFrame()
    df_sRNA_all = df_all[df_all['Site'].str.contains('|'.join(sRNA_all))]
    if len(set(smiRNA_lst)) != 0:
        df_smiRNA = df_all[df_all['Site'].str.contains('|'.join(list(set(smiRNA_lst))))].copy()
    else:
        df_smiRNA = df_all[df_all['Site'].str.contains('None')].copy()
    if len(set(sRNA_all) - set(smiRNA_lst)) != 0:
        df_sRNA_other = df_all[df_all['Site'].str.contains('|'.join(list(set(sRNA_all) - set(smiRNA_lst))))].copy()
    else:
        df_sRNA_other = df_all[df_all['Site'].str.contains('None')].copy()
    df_smiRNA['Type'] = ['miRNA'] * df_smiRNA.shape[0]
    df_sRNA_other['Type'] = ['sRNA'] * df_sRNA_other.shape[0]
    df_sRNA = pd.concat([df_sRNA, df_smiRNA, df_sRNA_other])
    df_sRNA['Sequence'] = df_seq[df_sRNA['Site'].tolist()].loc[0, ].tolist()
    df_sRNA.to_csv(FindPeak_dir + '/Result.txt', header=True, index=None, sep='\t')
else:
    print('No Data')
print('RNAfold finish time:', time.asctime(time.localtime(time.time())))
