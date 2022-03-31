#!/usr/bin/env python
"""
@author: czn
"""
import os
import argparse
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import time
print('Start smoothing:', time.asctime(time.localtime(time.time())))
parser = argparse.ArgumentParser()
parser.add_argument('--infile', help='depth file')
parser.add_argument('--outdir', help='output path')
parser.add_argument('--windows_size', help='windows size of smoothing', type=int)
parser.add_argument('--chr', help='specific chromosome', type=str)
args = parser.parse_args()
inputfile = args.infile
windows_size = args.windows_size
outdir = args.outdir+'/'
chr = args.chr

outfile = inputfile.split('/')[-1].split('.txt')[0]+'.extend.txt'
out = open(outdir+outfile, 'w')
gap = 100
tag = 0
with open(inputfile) as file:
    for line in file:
        line1 = line.strip().split('\t')
        if line1[0] == chr:
            if int(float(line1[3])) == 0:
                tag = 1
                zero_start, zero_end = int(line1[1]), int(line1[2])
                zero_length = zero_end - zero_start
            else:
                nonzero_start, nonzero_end = int(line1[1]), int(line1[2])
                if tag == 1:
                    if gap >= zero_length:
                        for i in range(zero_start+1, zero_end+1):
                            out.write(line1[0]+'\t'+str(i)+'\t0'+'\n')
                    else:
                        if zero_start != 0:
                            for i in range(zero_start + 1, zero_start + int(gap/2) + 1):
                                out.write(line1[0]+'\t'+str(i)+'\t0'+'\n')
                            for i in range(zero_end - 1 - int(gap/2), zero_end + 1):
                                out.write(line1[0]+'\t'+str(i)+'\t0'+'\n')
                        else:
                            for i in range(zero_end - int(gap/2) - 1, zero_end + 1):
                                out.write(line1[0]+'\t'+str(i)+'\t0'+'\n')
                for i in range(nonzero_start+1, nonzero_end+1):
                    out.write(line1[0]+'\t'+str(i)+'\t'+str(int(float(line1[3])))+'\n')
                tag = 0
out.close()
size = os.path.getsize(outdir+outfile)
if size != 0:
    df = pd.read_csv(outdir+outfile, sep='\t', header=None)
    sites = df[1].tolist()
    depth = df[2].tolist()
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth
    y_smooth = smooth(depth, int(windows_size))
    df['y_smooth'] = y_smooth
    print('Finish time:', time.asctime(time.localtime(time.time())))
    df['y_smooth'] = round(df['y_smooth'], 4)
    df.to_csv(outdir+outfile + '.smooth.windows_size.' + str(windows_size) + '.txt', sep='\t', index=None)
else:
    print('No smooth data')