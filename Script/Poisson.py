#!/usr/bin/env python
"""
@author: czn
"""
import re, os
import pandas as pd
import numpy as np
from scipy.stats import poisson
from collections import Counter
import time
import argparse
print('Filter by count or rpm, length and poisson:', time.asctime(time.localtime(time.time())))
parser = argparse.ArgumentParser()
parser.add_argument('--method', help='filter by count or rpm')
parser.add_argument('--rpm', help='RPM threshold of boundry sites', type=float)
parser.add_argument('--count', help='count threshold of boundry sites', type=int)
parser.add_argument('--minlen', help='minimum length of sRNA, default=15', default=15, type=int)
parser.add_argument('--maxlen', help='maximum length of sRNA, default=34', default=34, type=int)
parser.add_argument('--mapInfo', help='bowtie mapInfor file')
parser.add_argument('--cleanread', help='the number of clean reads', type=int)
parser.add_argument('--genomecovfile', help='genome coverage file')
parser.add_argument('--sRNAfile', help='sRNA file')
parser.add_argument('--pvalue', help='pvalue threshold based on the Poisson distribution', type=float)
parser.add_argument('--chr', help='specific chromosome', type=str)
parser.add_argument('--extendregion', help='extending the upstream and downstream 10 nt of the sRNA is considered as '
                                           'the background region, default: 10 nt', default=10, type=int)
parser.add_argument('--output', help='output file')
parser.add_argument('--threads', help='running threads', type=int)

args = parser.parse_args()
method = args.method
rpm = args.rpm
mincount = args.count
minlen, maxlen = args.minlen, args.maxlen
mapInfo = args.mapInfo
cleanread = args.cleanread
genomecovfile = args.genomecovfile
sRNAfile = args.sRNAfile
pvalue = args.pvalue
chr = args.chr
extendregion = args.extendregion
ouputfile = args.output
threads = args.threads

df_genomecov_all = pd.DataFrame()
df_sRNA = pd.DataFrame()
df_genomecov = pd.read_csv(genomecovfile, header=None, sep='\t')
df_genomecov_all = pd.concat([df_genomecov_all, df_genomecov], axis=0, ignore_index=True)
if os.path.exists(sRNAfile):
    df_sRNA = pd.read_csv(sRNAfile, header=0, sep='\t')
    df_sRNA['Chr'] = [chr]*df_sRNA.shape[0]
    df_sRNA = df_sRNA[df_sRNA['count'] != 'Nan']
    df_sRNA['Length'] = df_sRNA['Length'].astype(int)
    df_sRNA['Start_count'] = df_sRNA['count'].map(lambda x: x.split(',')[0])
    df_sRNA['End_count'] = df_sRNA['count'].map(lambda x: x.split(',')[1])
    df_sRNA['Start_rpm'] = (df_sRNA['count'].map(lambda x: x.split(',')[0]).astype(int) / cleanread) * 1000000
    df_sRNA['End_rpm'] = (df_sRNA['count'].map(lambda x: x.split(',')[1]).astype(int) / cleanread) * 1000000

    """ Filtered by count or rpm """
    if method == 'count':
        df_filter = (df_sRNA[(df_sRNA['Start_count'].astype(int) >= mincount) & (df_sRNA['End_count'].astype(int) >= mincount)])
    elif method == 'rpm':
        df_filter = (df_sRNA[(df_sRNA['Start_rpm'] >= rpm) & (df_sRNA['End_rpm'] >= rpm)])
    df_filter.insert(0, 'Site', (df_filter['Start_min'].astype(str) + '-' + df_filter['End_max'].astype(str) +
                                 '-' + df_filter['Strand'] + '-' + df_filter['Chr'].tolist()).tolist())
    df_filter = df_filter.drop_duplicates(subset='Site', keep='first')
    """ Filtered by length """
    df_filter = (df_filter[(df_filter['Length'].astype(int) >= minlen) & (df_filter['Length'].astype(int) <= maxlen)])

    def GetSigSite(site, df_filter, df_genomecov_all):
        count1, count2 = df_filter[df_filter['Site'] == site]['Start_count'].tolist()[0], \
                         df_filter[df_filter['Site'] == site]['End_count'].tolist()[0]
        start, end = int(site.split('-')[0]), int(site.split('-')[1])
        chr = site.split('-')[3]
        rawsite_lst = list(range(start, end + 1))
        df_genomecov = df_genomecov_all[df_genomecov_all[0].astype(str) == chr]
        extend_start, extend_end = start - extendregion, end + extendregion
        site_lst = list(range(extend_start, extend_end + 1))
        lam = np.average(df_genomecov[df_genomecov[1].isin(site_lst)][2])
        p_lst = []
        for j in df_genomecov[df_genomecov[1].isin(rawsite_lst)][2]:
            p_lst.append(poisson.pmf(j, int(lam)))
        p = max(p_lst)
        p_lst_filter = [i for i in p_lst if i <= 0.05]
        percentage = len(p_lst_filter) / len(p_lst)
        if p <= pvalue or (p > pvalue and percentage >= 0.8):
            return site

    from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
    p = ProcessPoolExecutor(threads)
    site_result = []
    for site in df_filter['Site'].tolist():
        site_result.append(p.submit(GetSigSite, site, df_filter, df_genomecov_all))
    p.shutdown()
    sig_site = []
    for i in site_result:
        sig_site.append(i.result())

    df_filter_sig = df_filter[df_filter['Site'].isin(sig_site)]
    column = ["Site", "Chr", "Start_min", "End_max", "Strand", "Length", "Start_count", "End_count", "Start_rpm", "End_rpm"]
    df_filter_sig_ouput = df_filter_sig[column].copy()
    df_filter_sig_ouput.columns = ["Site", "Chr", "Start", "End", "Strand", "Length", "Start_count", "End_count", "Start_rpm", "End_rpm"]
    df_filter_sig_ouput.to_csv(ouputfile, sep='\t', index=None)
    column = ["Site", "Chr", "Start_min", "End_max", "Strand"]
    df_saf = df_filter_sig[column].copy()
    df_saf.columns = ["GeneID", "Chr", "Start", "End", "Strand"]
    df_saf['Strand'].replace('positive', '+', inplace=True)
    df_saf['Strand'].replace('negative', '-', inplace=True)
    df_saf.to_csv(ouputfile+'.saf', sep='\t', index=None)
    print('small RNA region', df_filter[df_filter['Site'].isin(sig_site)].shape[0])
