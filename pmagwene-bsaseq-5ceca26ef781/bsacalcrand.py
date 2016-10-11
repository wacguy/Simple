# -*- coding: utf-8 -*-

import sys
import argparse
import csv

import numpy as np
#from scipy import signal

import kernelsmoothing as ks
from kernelsmoothing import  uniform, triangle, tricube, triweight, quartic, connected_maxminfilter
kernel_dict = {'uniform':uniform, 'triangle':triangle, 'tricube':tricube, 
                'triweight':triweight, 'quartic':quartic}




def read_alleles(fname):
    reader = csv.reader(open(fname,'rU'),delimiter="\t")
    lines = [row for row in reader]
    cset = set(((i[0],i[1]) for i in lines))  # in python 2.7 could use set comprehensions
    return lines, cset


def filter_lines(lines, ccset):
    return [i for i in lines if tuple(i[:2]) in ccset]

def write_lines(lines, fname):
    writer = csv.writer(open(fname, 'w'), delimiter="\t")
    writer.writerows(lines)



def load_allele_cts(fname):
    """ Parse a tab delimited file of allele counts at biallelic loci.

    Returns:
    
    1. a list of arrays, each array representing the parsed data for a 
    given chromosome.
    2. a list with the original chromsome names

    Expects a tab delimited file consisting of 4 columns, where the 
    columns are:

    1. chromosome names -- treated as a string
    2. chromosomal coordinate -- integer
    3. observed count of allele A -- integer
    4. observed count of allele a -- integer
    """
    reader = csv.reader(open(fname,'rU'),delimiter="\t")
    
    chroms = []
    alleles = []
    thechrom = []
    for row in reader:
        if row[0] not in chroms:  # if is a new chromosome
            if len(thechrom): # append the old chromosome to the alleles list
                alleles.append(np.array(thechrom[:]))
            chroms.append(row[0]) 
            thechrom = [] # reset the current chromosome
        thechrom.append((chroms.index(row[0]), int(row[1]),int(row[2]),int(row[3])))
    alleles.append(np.array(thechrom)) # append the last chromosome
    
    return alleles



def join_bulks(bulk1, bulk2):
    """ A convenience function to join two bulks into one.
    """
    chroms = []
    nchroms = len(bulk1)
    for i in range(nchroms):
        b1, b2 = bulk1[i], bulk2[i]
        thechrom = np.column_stack((b1[:,:2],b1[:,2:]+b2[:,2:]))
        chroms.append(thechrom)
    return chroms



def Gtest_indep(a, b, c, d):
    """Calculates G-test for independence of 2 x 2 tables.
    
    see Sokal and Rohlf 1994, eqn 17.11 and Box 17.6 p. 731
    """
    a = np.array(a,dtype=np.float)
    b = np.array(b,dtype=np.float)
    c = np.array(c,dtype=np.float)
    d = np.array(d,dtype=np.float)
    log = np.log
    n = a + b + c + d
    q1 = a*log(a) + b*log(b) + c*log(c) + d*log(d) + n*log(n)
    ab = a + b
    cd = c + d
    ac = a + c
    bd = b + d
    q2 = ab*log(ab) +  ac*log(ac) + bd*log(bd) + cd*log(cd)  
    G = 2.0 * (q1 - q2)
    try:
        G[G < 0] = 0 # possible due to rounding error
    except TypeError: # can happen if input are not arrays
        pass
    return G

def Gtest_onepool(a, b, pi=0.5):
    """Calculates G-statistic based on observed counts a and b, relative to expected proportion pi.
    """
    a = np.array(a,dtype=np.float)
    b = np.array(b,dtype=np.float)
    n = a + b
    expected_a = n * pi
    expected_b = n * (1.0-pi)
    return 2 *  (a * np.log(a/expected_a) + b * np.log(b/expected_b))



def abs_diff_allelefreq(low, high):
    reflow = low[:,2].astype(np.float)/np.sum(low[:,2:], axis=1)
    refhigh = high[:,2].astype(np.float)/np.sum(high[:,2:], axis=1)
    return np.abs(reflow - refhigh)
    

def bulkG(low, high, const=0.01):
    """Given arrays representing allele counts in two bulks calculate G at each locus.
    
    low and high are (n x 4) numpy arrays with allele counts in the columns
    """
    a,b = low[:,2] + const, low[:,3] + const
    c,d = high[:,2] + const, high[:,3] + const  
    G =  Gtest_indep(a,b,c,d)
    return G

def bulkG_randomized(low, high, const=0.01):
    """Given arrays representing allele counts in two bulks calculate G at each locus.
    
    low and high are (n x 4) numpy arrays with allele counts in the columns
    """
    Low = low.copy()
    High = high.copy()
    np.random.shuffle(Low)
    np.random.shuffle(High)
    a,b = Low[:,2] + const, Low[:,3] + const
    c,d = High[:,2] + const, High[:,3] + const  
    G =  Gtest_indep(a,b,c,d)
    return G


def bulkG_onepool(low, pi=0.5, const=0.01):
    """Given arrays representing allele counts in two bulks calculate G at each locus.
    
    low and high are (n x 4) numpy arrays with allele counts in the columns
    """
    a,b = low[:,2] + const, low[:,3] + const
    G =  Gtest_onepool(a,b,pi)
    return G    



def remove_jointly_fixed(low, high):
    """ Remove sites that are fixed for the same allele in both the low and high samples.
    """
    notfixed = []
    for i in range(len(low)):
        lowfixed = is_fixed(low[i,2], low[i,3])
        highfixed = is_fixed(high[i,2], high[i,3])
        if (lowfixed[0] == True) and (lowfixed == highfixed):
            continue
        else:
            notfixed.append(i)
    return low[notfixed,:], high[notfixed,:]

def remove_all_jointly_fixed(lowchroms, highchroms):
    newlow, newhigh = [], []
    for i in range(len(lowchroms)):
        l, h = remove_jointly_fixed(lowchroms[i], highchroms[i])
        newlow.append(l)
        newhigh.append(h)
    return newlow, newhigh


def is_fixed(refct, altct):
    if (refct == 0) and (altct > 0):
        return True, 0
    elif (altct == 0) and (refct > 0):
        return True, 1
    else:
        return False, None



class Results:
    pass
    
def Gstats(low, high, halfwidth=10000, kernel=tricube, maxgap=10000, maxminsize=7):
    """ Given low and high bulks, calculate G and G'.
    """
    nchroms = len(low)
    G = [bulkG(l, h) for l,h in zip(low,high)]
    r = Results()
    r.chroms, r.coords, r.G, r.smoothG, r.maxminG = [],[],[],[],[]
    for i in range(nchroms):
        chromi = low[i][:,0]
        coordsi = low[i][:,1]
        Gi = G[i]
        smoothGi = np.array(ks.kernel_smooth(coordsi, Gi, halfwidth, kernel=kernel))
        maxmincoords, maxminG = connected_maxminfilter(coordsi, Gi, size=maxminsize, maxgap=maxgap)
        r.chroms.append(chromi)
        r.coords.append(coordsi)
        r.G.append(Gi)
        r.smoothG.append(smoothGi)
        r.maxminG.append(maxminG)
    return r

def Gstats_randomized(low, high, halfwidth=10000, kernel=tricube, maxgap=10000, maxminsize=7):
    """ Given low and high bulks, calculate G and G'.
    """
    nchroms = len(low)
    G = [bulkG_randomized(l, h) for l,h in zip(low,high)]
    r = Results()
    r.chroms, r.coords, r.G, r.smoothG, r.maxminG = [],[],[],[],[]
    for i in range(nchroms):
        chromi = low[i][:,0]
        coordsi = low[i][:,1]
        Gi = G[i]
        smoothGi = np.array(ks.kernel_smooth(coordsi, Gi, halfwidth, kernel=kernel))
        maxmincoords, maxminG = connected_maxminfilter(coordsi, Gi, size=maxminsize, maxgap=maxgap)
        r.chroms.append(chromi)
        r.coords.append(coordsi)
        r.G.append(Gi)
        r.smoothG.append(smoothGi)
        r.maxminG.append(maxminG)
    return r

def Gstats_onepool(low, pi=0.5, halfwidth=10000, kernel=tricube, maxgap=10000, maxminsize=7):
    """ Given one bulk, calculate G and G'.
    """
    nchroms = len(low)
    G = [bulkG_onepool(l, pi) for l in low]
    r = Results()
    r.chroms, r.coords, r.G, r.smoothG, r.maxminG = [],[],[],[],[]
    for i in range(nchroms):
        chromi = low[i][:,0]
        coordsi = low[i][:,1]
        Gi = G[i]
        smoothGi = np.array(ks.kernel_smooth(coordsi, Gi, halfwidth, kernel=kernel))
        maxmincoords, maxminG = connected_maxminfilter(coordsi, Gi, size=maxminsize, maxgap=maxgap)
        r.chroms.append(chromi)
        r.coords.append(coordsi)
        r.G.append(Gi)
        r.smoothG.append(smoothGi)
        r.maxminG.append(maxminG)
    return r






def save_Gstats(fname, results, withmaxmin=False):
    """ Save the results of the analysis to a text file.
    """
    chroms = np.hstack(results.chroms)
    coords = np.hstack(results.coords)
    G = np.hstack(results.G)
    smoothG = np.hstack(results.smoothG)
    if withmaxmin:
        maxminG = np.hstack(results.maxminG)
        sresults = np.column_stack((chroms,coords,G,smoothG,maxminG))
        fmt = "%2.1d\t%5.1d\t%5.4f\t%5.4f\t%5.4f"
    else:
        sresults = np.column_stack((chroms,coords,G,smoothG))
        fmt = "%2.1d\t%5.1d\t%5.4f\t%5.4f"

    np.savetxt(fname, sresults, fmt=fmt) 




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BSA mapping from high-throughput sequencing")

    parser.add_argument('-L','--low', nargs='*', required=True, 
            help="one or more files giving allele freqs in low bulk")
    parser.add_argument('-H','--high',nargs='*', required=True,
            help="one or more files giving allele freqs in high bulk")
    parser.add_argument('-k','--kernel', type=str, default="tricube", 
            choices=['uniform','tricube','triweight','quartic'], 
            help='smoothing kernel (default=tricube)')             
    parser.add_argument('-w','--width', type=int, default=20000, 
            help='smoothing kernel half-width (default=20000)')    
    parser.add_argument('-o', '--outfile',nargs='?', type=argparse.FileType('w'), default=sys.stdout)

    parser.add_argument('--filterfixed', action='store_true', dest='jointfixed', default=False,
                help="Whether to filter out sites that are fixed for same allele in both low and high bulks")

    parser.add_argument('--maxmin', action='store_true', dest='maxmin', default=False,
                help="Include the moving maxmin values of the G statistic in the output")

    parser.add_argument('--onepool', action='store_true', dest='onepool', default=False,
                help="Whether to do a one pool G statistic (using Low bulk)")

    args = parser.parse_args()

    lbulks, hbulks = [],[]
    
    for l in args.low:
        bulk = load_allele_cts(l)
        lbulks.append(bulk)
    
    lowbulk = lbulks[0]    
    if len(lbulks) > 1:
        for bulk in lbulks[1:]:
            lowbulk = join_bulks(lowbulk, bulk)    

    for h in args.high:
        bulk = load_allele_cts(h)
        hbulks.append(bulk)        

    highbulk = hbulks[0]    
    if len(hbulks) > 1:
        for bulk in hbulks[1:]:
            highbulk = join_bulks(highbulk, bulk) 

    if args.jointfixed and not args.onepool:
        newlow, newhigh = [], []
        for i in range(len(lowbulk)):
            l, h = remove_jointly_fixed(lowbulk[i], highbulk[i])
            newlow.append(l)
            newhigh.append(h)
        lowbulk = newlow
        highbulk = newhigh
    
    if args.onepool:
        results = Gstats_onepool(lowbulk, 0.5, halfwidth=args.width, kernel=kernel_dict[args.kernel]) 
    else:
        results = Gstats(lowbulk, highbulk, halfwidth=args.width, kernel=kernel_dict[args.kernel]) 

    save_Gstats(args.outfile, results, withmaxmin=args.maxmin)


