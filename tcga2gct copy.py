#!/usr/bin/python3

# Converting TCGA-style text files of gene expression from FireBrowse to GCT files

import sys
#import csv
import pandas as pd
import numpy as np

def main(args):
    if(args[1]=='-h' or args[1]=='--help' or args[1]=='--h'):
        print("\nConverts TCGA-downloaded expression files (.txt) to GenePattern-compatible files (.gct)")
        print("\tUsage: tcga2gct /path/to/tcgafile.txt [/path/to/gctoutput.gct]\n")
        #return
    
    tissue=getTissue(args[1])
    fname=tissue+'_final.gct'
    
    df_i=pd.read_table(args[1], low_memory=False)
    df=df_i.drop(df_i.index[[0]])
    sample_df=df.iloc[:,1:]
    samples=df.columns[1:]

    genes=df.iloc[:,0]
    cgenes=cleanGenes(genes)
    
    sample_df.insert(0,'genes',pd.Series(cgenes, sample_df.index))
    sample_df.insert(1,'description',pd.Series("NONE", sample_df.index))
    
    with open(fname, 'w') as outfile:
        outfile.write('#1.2\n')
        outfile.write(str(len(cgenes))+'\t'+str(len(samples))+'\n')
        
    with open(fname, 'a') as outfile:
        sample_df.to_csv(outfile, sep='\t', header=False, index=False, mode='a')
        


def cleanGenes(glist):
    glist=list(glist)
    cglist=[]
    iq=0
    for g in glist:
        g=g.split(sep='|')[0]
        if(g=='?'):
           g='?'+'.'+str(iq)
           iq+=1
        cglist.append(g.upper())
    return(cglist)

def getTissue(fname):
    fns=fname.split(sep='/')
    ttype=fns[len(fns)-1].split(sep='.')
    return(ttype[0])

main(sys.argv)