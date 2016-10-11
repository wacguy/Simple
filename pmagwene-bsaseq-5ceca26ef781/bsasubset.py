
import sys
import csv
import argparse
import os, os.path




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Utility for finding common SNPs from a set of input files.")

    parser.add_argument('-i', '--input', nargs='+', type=str, required=True,
            help='input files to use to find common subset')

    parser.add_argument('-s','--suffix', type=str, default="-filtered.out", 
            help='string to add as suffix to the output file name') 

    args = parser.parse_args()

    outsuffix = args.suffix
    if outsuffix == "":
        outsuffix = "-filtered.out"

    infiles = args.input

    fname0 = infiles[0]
    ccset = set()
    with open(fname0, 'rU') as f:
        for line in f:
            words = line.split()
            ccset.add( (words[0], words[1]) )


    for fname in infiles[1:]:
        cc = set()
        with open(fname, 'rU') as f:
            for line in f:
                words = line.split()
                cc.add( (words[0], words[1]) )
            ccset &= cc

    for fname in infiles:
        with open(fname, 'rU') as f:
            lines = []
            for line in f:
                words = line.split()
                if (words[0], words[1]) in ccset:
                    lines.append(line)
            inprefix, inext = os.path.splitext(fname)
            outname = inprefix + outsuffix
            outfile = open(outname, 'w')
            outfile.writelines(lines)
            outfile.close() 
                            