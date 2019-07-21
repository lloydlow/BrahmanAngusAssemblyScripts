# -*- coding: utf-8 -*-
"""
Created on Mon May  6 15:12:20 2019

@author: dbickhart
This script is provided by Derek M Bickhart. It is free to use for academic/non-profit work in the hopes it will be useful!
PS: I'm very anti-pythonic, so approach the code with an open-mind
"""

import argparse
from collections import defaultdict
import numpy as np
import contextlib
import os

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Generate restriction enzyme maps for an assembly fasta and select fragments of a specific size range"
            )
    parser.add_argument('-c', '--gene', 
                        help="Input gene copy number list from annotation program",
                        type=str, required=True
                        )
    parser.add_argument('-p', '--pop', 
                        help="Population structure (sample(tab)popnumber) file",
                        type=str, required=True
                        )
    parser.add_argument('-o', '--output', 
                        help="Output file name",
                        type=str, required=True
                        )
    parser.add_argument('-s', '--skip', 
                        help="Skip this sample (can be specified more than once)",
                        action="append", default=[]
                        )
    parser.add_argument('-m', '--melt',
                        help="Write CN data to melted table for plotting in R",
                        type=str, default='none'
                        )
    parser.add_argument('-v', '--vfilt',
                        help="VST filter for melt file",
                        type=float, default=0.2
                        )
    parser.add_argument('-g', '--mag',
                        help="Minimum of difference between mean of CN values to report",
                        type=float, default=2.0
                        )
    return parser.parse_args()

def main(args):
    # Read in population structure and skip list
    skip = set()
    if len(args.skip) > 0:
        skip = {x for x in args.skip}
    
    popstruct = defaultdict(list)
    poplookup = {}
    with open(args.pop, 'r') as fh:
        for l in fh:
            l = l.rstrip()
            s = l.split()
            if s[0] in skip: continue
            poplookup[s[0]] = s[1]
            popstruct[s[1]].append(s[0])
            
    # Now to process the data on a line-by-line basis
    numericErrs = 0
    with open(args.gene, 'r') as fh, open(args.output, 'w') as out, smartFile(args.melt) as melt:
        # Process header to get the column indicies for the samples
        head = fh.readline()
        geno = genoLookup(poplookup, head, skip)
        if args.melt != 'none':
            melt.write('\t'.join(['Gene', 'Pop', 'Sample', 'CN']) + '\n')
        for l in fh:
            l = l.rstrip()
            s = l.split()
            if not isNumber(s[6]): 
                # Some of the entries are "null" because of the difference in chromosome names
                numericErrs += 1
                continue
            marker = s[1] + ':' + s[2] + '-' + s[3]
            worker = totalHolder(marker, s[6:], s[0])
            
            Vst = worker.Vst(geno, marker, skip)
            out.write(f'{s[1]}\t{s[2]}\t{s[3]}\t{Vst}\t{s[0]}\n')
            
            mtest = worker.magTest(args.mag)
            if args.melt != 'none' and Vst >= args.vfilt and mtest:
                print(f'Melt\t{worker.genes}\t{Vst}\t{worker.mag}')
                for l in worker.decorateCNData(geno, Vst, skip):
                    melt.write(l + '\n')
            
    print(f'Dealt with {numericErrs} null fields')

@contextlib.contextmanager
def smartFile(filename : str, mode : str = 'w'):
    if filename == 'none':
        fh = open(os.devnull, 'w')
    else:
        fh = open(filename, mode)
    try:
        yield fh
    finally:
        if filename != 'none':
            fh.close()
    
def isNumber(value : str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False

class totalHolder:
    
    def __init__(self, markername: str, values: list, gene: str):
        self.markername = markername
        self.values = []
        self.values.append([float(x) for x in values])
        self.genes = gene
        self.popavgs = defaultdict(float)
        self.mag = 0.0
    
    # Routines for overlapping gene regions
    def addValues(self, values: list):
        self.values.append(values)
        
    def addGene(self, gene: str):
        self.genes += ';' + gene
   
    def Vst(self, genoLookup, markers, skip) -> float:
        # Total variance calculation
        total = []
        numerator = 0
        denominator = 0
        # Data structures
        popnums = {}
        popvar = defaultdict(list)
        for i in range(0, len(self.values)):
            for j in range(0, len(self.values[i])):                
                sample, pop = genoLookup.getPopSamp(j)
                if sample == None:
                    continue
                
                total.append(self.values[i][j])
                if not pop in popnums:
                    popnums[pop] = 0
                popnums[pop] += 1
                popvar[pop].append(self.values[i][j])
                self.popavgs[pop] += self.values[i][j]
                
        for pop, num in popnums.items():
            pvar = np.var(popvar[pop])
            numerator += pvar * num
            denominator += num
            self.popavgs[pop] /= num
            
        Vt = np.var(total)
        Vs = numerator / denominator
        if Vt == 0: return 0.0
        
        return ((Vt - Vs) / Vt)     

    def magTest(self, mag):
        # Test if the difference between pop averages is above the difference
        diff = 0
        for k, v in self.popavgs.items():
            diff = abs(diff - v)

        self.mag = diff
        return diff > mag
    
    def decorateCNData(self, genoLookup, Vst, skip):
        # Print out the data in "melted form"
        geneid = '{} ({:0.2f})'.format(self.genes, Vst)
        content = []
        for i in range(0, len(self.values)):
            for j in range(0, len(self.values[i])):
                sample, pop = genoLookup.getPopSamp(j)
                if sample == None:
                    continue
                
                subc = [geneid, pop, sample, str(self.values[i][j])]
                content.append('\t'.join(subc))
        
        return content
        
        
class genoLookup:
    
    def __init__(self, poplookup, head, skip):
        self.pops = {}
        self.samples = {}
        head = head.rstrip()
        header = head.split()
        
        for i in range(6, len(header)):
            if header[i] in skip:
                self.samples[i - 6] = None
                self.pops[i - 6] = None
            else: 
                self.samples[i - 6] = header[i]
                self.pops[i - 6] = poplookup[header[i]]
            
    def getPopSamp(self, idx: int):
        return self.samples[idx], self.pops[idx]
    
    def getOrderedSamp(self):
        return self.samples
    
    def getOrderedPops(self):
        return self.pops

    
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
