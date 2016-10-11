import math
import webbrowser

import numpy as np
from scipy import stats

import kernelsmoothing as ks
from kernelsmoothing import  uniform, triangle, tricube, triweight, quartic
kernel_dict = {'uniform':uniform, 'triangle':triangle, 'tricube':tricube,
                'triweight':triweight, 'quartic':quartic}

import bsacalc





def haldane(x):
    """ Haldane's function.
    """
    z = (0.5 * (1.0 - np.exp(-2 * (np.abs(x)))))
    return z

def linkageF(x, fA=1):
    r = haldane(x)
    return fA*(1-2*r)+r


def linkageG(x, fL=0.5, fH=0.5):
    log2 = np.log(2)
    fBl = linkageF(x, fL)
    notBl = 1 - fBl

    fBh = linkageF(x, fH)
    notBh = 1 - fBh

    gL = fBl * (np.log(fBl) - np.log(notBl)) + np.log(notBl) +  log2
    gH = fBh * (np.log(fBh) - np.log(notBh)) + np.log(notBh) +  log2
    return gL + gH




def theoryG_exact(ns, C, v, maphalfdist, kernel=tricube):
    """ Calculate E[G], Var[G], and Var[G'] based on theoretical considerations.

    ns = number of segregants
    C = average coverage per variable site
    v = effective number of sites per window
    maphalfdist = smoothing window half-width in centimorgans

    See Magwene, Willis, Kelly 2011.
    """
    a = 2.0 + (1./(2.*C)) + ((1.0 + 2*C)/ns)
    b =  (C**2 * (4.0*ns-1.0))/(8.0 * ns**3)

    E_G = 1.0 + (C/(2.0*ns))
    Var_G = a + b

    k = get_k(v, kernel)
    k2sum = np.sum(k*k)

    r = get_r(maphalfdist, v)
    rksum = sum_over_kandr(k, r)

    Var_Gbar = (a+b)*k2sum + (b * rksum)

    return E_G, Var_G, Var_Gbar


def get_k(npts, kernel=tricube):
    q =  kernel(np.linspace(-1,1,npts))
    return q/np.sum(q)

def get_r(maphalfdist,npts):
    m = np.linspace(-maphalfdist,maphalfdist,npts)
    dists = np.array([np.abs(i-m) for i in m])
    return haldane(dists)

def sum_over_kandr(k,r):
    inners = []
    for i in range(len(k)):
        ri = r[i]
        ksub = np.concatenate((k[:i],k[i+1:]))
        rsub = np.concatenate((ri[:i],ri[i+1:]))
        innersum = np.sum(k[i] * ksub * (1.0-2*rsub)**2)
        inners.append(innersum)
    return np.sum(inners)

def theoryG_exact_alt(ns, C, v, maphalfdist, kernel=tricube):
    a = 2.0 + (1./(2.*C)) + ((1.0 + 2*C)/ns)
    b =  (C**2 * (4.0*ns-1.0))/(8.0 * ns**3)

    E_G = 1.0 + (C/(2.0*ns))
    Var_G = a + b

    k = get_k(v, kernel)
    k2sum = np.sum(k*k)

    r = get_r(maphalfdist, v)
    #m = np.abs(np.linspace(-maphalfdist,maphalfdist,v))
    #r = np.array([np.abs(i-m) for i in m])

    #rksum = sum_over_kandr(k, r)

    inners = []
    for i in range(len(k)):
        ri = r[i]
        ksub = np.concatenate((k[:i],k[i+1:]))
        rsub = np.concatenate((ri[:i],ri[i+1:]))
        innersum = np.sum(b * k[i] * ksub * (1.0-2*rsub)**2)
        inners.append(innersum)
    rksum = np.sum(inners)

    Var_Gbar = (a+b)*k2sum + rksum

    return E_G, Var_G, Var_Gbar


def theoryG(C,ns,v):
    a = 2.0 + (1./(2.*C)) + ((1.0 + 2*C)/ns)
    b = 350./(247.*v)  # for tri-cubed
    c =  (C**2 * (4.0*ns-1.0))/(8.0 * ns**3)

    E_G = 1.0 + (C/(2.0*ns))
    Var_G = a + c

    Var_Gbar = (a*b)+c
    return E_G, Var_G, Var_Gbar

def theoryG_diploidhomozygous(C,ns,v):
    a = 2.0 + (1./(2.*C)) + ((1.0 + 2*C)/(0.5*ns))
    #b = (8.*(64-112*v**2 + 27*v**6))/((21. * v**3) * ((4.-3*v**2)**2))  ## for cubed
    b = 350./(247.*v)  # for tri-cubed
    c =  (C**2 * (2.0*ns-1.0))/(ns**3)

    E_G = 1.0 + (C/float(ns))
    Var_G = a + c

    Var_Gbar = (a*b)+c
    return E_G, Var_G, Var_Gbar

def theoryG_old(C,ns,W,meanR):
    p3 =  (C**2 * (4.0*ns-1.0))/(8.0 * ns**3)
    a = 2.0 + (1./(2.*C)) + ((1.0 + 2*C)/ns) + p3
    E_G = 1.0 + (C/(2.0*ns))
    Var_G = a
    b =  p3 * (1.0-meanR)**2
    Var_Gbar = (a/W) + b
    return E_G, Var_G, Var_Gbar

def theoryG_haploid(C,ns,W,meanR):
    p3 =  (C**2 * (2*ns-1.0))/(2.0 * ns**3)
    a = 2.0 + (1./(2.*C)) + ((1. + 2*C)/ns) + p3
    E_G = 1.0 + C/float(ns)
    Var_G = a
    b =  p3 * (1.0-meanR)**2
    Var_Gbar = (a/W) + b
    return E_G, Var_G, Var_Gbar



def MAD(x,axis=None):
    """ Calculate the Median Absolute Deviation.
    """
    x = np.array(x)
    medx = np.median(x,axis=axis)
    return np.median(np.abs(x-medx), axis=axis)


def leftMAD(x,axis=None):
    x = np.array(x)
    medx = np.median(x, axis=axis)
    left = x[x < medx]
    return np.median(np.abs(left-medx))


def compute_half_sample(x):
    if len(x) <= 3:
        return x
    half = len(x)/2
    n = len(x)
    s = sorted(x)

    minsample = s[0:half]
    minrange = s[-1] - s[0]

    for i in range(0, n - half + 1):
        samp = s[i:i+half]
        srng = samp[-1] - samp[0]
        if srng < minrange:
            minsample = samp
            minrange = srng
    return minsample


def half_sample_mode(x):
    """Robust estimator of the mode of a continuous distn.

    See D. R. Bickel and R. Fruhwirth  On a Fast, Robust Estimator of the Mode:
    Comparisons to Other Robust Estimators with Applications,
    Computational Statistics and Data Analysis 50, 3500-3530 (2006)
    """
    x = sorted(x)
    if len(x) == 1:
        return x[0]
    if len(x) == 2:
        return sum(x)/2.0
    if len(x) == 3:
        return min( [(x[2] + x[1])/2., (x[1] + x[0])/2.] )

    big = True
    while big:
        x = compute_half_sample(x)
        if len(x) > 3:
            continue
        else:
            if len(x) == 1:
                return x[0]
            if len(x) == 2:
                return sum(x)/2.0
            if len(x) == 3:
                return min( [(x[2] + x[1])/2., (x[1] + x[0])/2.] )



def lognormal_mean_var(Ex, Vx):
    """ Calculate mu and sigma^2 for a log-normal distn based on E[X] and Var[X].
    """
    mu = math.log(Ex) - 0.5*math.log(1+Vx/Ex**2)
    s2 = math.log(1+Vx/Ex**2)
    return mu,s2




def identify_outliers(X, madmult=5.2):
    med = np.median(X)
    mode = half_sample_mode(X)
    ctr = (med + mode)/2.
    left = X[X < ctr]
    MAD = np.median(np.abs(left-ctr))
    non = np.less_equal(X, mode + madmult * MAD)
    return non

def robust_lognormal_mean_var(X, modefunc=half_sample_mode, madmult=5.2):
    """ Robust, non-parametric estimator of a 'contaminated' log-normal distn.

    We assume the contaminating distns represent less than 50% of the observations.

    We assume the mean of the contaminating distributions lies to the right
    of the null distribution.
    """
    logX = np.log(X)
    # logMED = np.median(logX)
    # logMAD = leftMAD(logX)
    #
    # nonoutlier = np.less_equal(logX, logMED + madmult*logMAD)
    # trimX = X[nonoutlier]
    trimX = X[identify_outliers(logX)]

    mu = math.log(np.median(trimX))
    mode = modefunc(trimX)
    s2 = mu - math.log(mode)
    return mu,s2



def estimate_smoothG_cutoff_lognorm_theory(smoothG, tmean, tvar, fdr=0.05, independent=False):
    tmean, tvar = float(tmean), float(tvar)
    lmu, ls2 = lognormal_mean_var(tmean, tvar)
    nulldist = stats.lognorm(math.sqrt(ls2),scale=math.exp(lmu))
    pvals = [nulldist.sf(i) for i in smoothG]
    pcutoff = BH_fdr_procedure(pvals, fdr, independent=independent)
    if pcutoff is None:
        return None, None
    gcutoff = nulldist.isf(pcutoff)
    return pcutoff, gcutoff

def estimate_smoothG_cutoff_lognorm_robust(smoothG, fdr=0.05, independent=False, modefunc=half_sample_mode, madmult=5.2):
    lmu, ls2 = robust_lognormal_mean_var(smoothG, modefunc=modefunc, madmult=madmult)
    nulldist = stats.lognorm(math.sqrt(ls2),scale=math.exp(lmu))
    pvals = [nulldist.sf(i) for i in smoothG]
    pcutoff = BH_fdr_procedure(pvals, fdr, independent=independent)
    if pcutoff is None:
        return None, None
    gcutoff = nulldist.isf(pcutoff)
    return pcutoff, gcutoff


def simulate_bulking_negbinom(freqsLow, freqsHigh, ns, C, coords=None):
    """Simulate generation of bulk samples.

    Under assumption that chromsomes are samples as Binomial and read depth is sampled
    as Negative Binomial
    """
    twoNs = 2.0 * ns
    nsnps = len(freqsLow)

    if coords is None:
        coords = range(nsnps)
    else:
        coords = sorted(coords)

    lowpop = []
    for i in range(2*ns):
        p = np.array([np.random.binomial(1,f) for f in freqsLow])
        lowpop.append(p)
    lowpop = np.array(lowpop)

    highpop = []
    for i in range(2*ns):
        p = np.array([np.random.binomial(1,f) for f in freqsHigh])
        highpop.append(p)
    highpop = np.array(highpop)

    lowP = np.sum(lowpop,axis=0)/twoNs
    highP = np.sum(highpop,axis=0)/twoNs

    lowref = []
    lowalt = []
    for i in range(nsnps):
        p = lowP[i]
        ref = int(C * p)
        alt = stats.nbinom.rvs(ref, 1-p)
        lowref.append(ref)
        lowalt.append(alt)

    highref = []
    highalt = []
    for i in range(nsnps):
        p = highP[i]
        ref = int(C * p)
        alt = stats.nbinom.rvs(ref, 1-p)
        highref.append(ref)
        highalt.append(alt)

    l = np.column_stack((lowref,lowalt))
    h = np.column_stack((highref,highalt))
    L, H, idx = simulated_bulks_to_samples(coords, l, h)
    #G = bsacalc.bulkG(L,H)

    return L, H, idx




def truth_runs(tflist, minlen=2):
    runs = []
    rstart, rend = None, None
    n = len(tflist)
    for i, e in enumerate(tflist):
        if e:
            if rstart is None:
                rstart, rend = i, i+1
            else:
                rend += 1
            if i == (n-1):
                if rstart != (n-1):
                    runs.append((rstart,n-1))
        else:
            if rstart is None:
                continue
            if (rend-rstart) >= minlen:
                runs.append((rstart,rend))
            rstart, rend = None, None
    return runs

def regions_above_threshold(gstats, threshold, minrun=10, minlen=10000):
    rdict = {}
    for i in range(len(gstats.chroms)):
        rdict[i] = []
        coords = gstats.coords[i]
        smoothG = gstats.smoothG[i]
        above = np.greater_equal(smoothG, threshold)
        runs = truth_runs(above, minrun)
        if len(runs):
            for (a,b) in runs:
                c1, c2 = coords[a], coords[b]
                if (c2-c1) < minlen:
                    continue
                maxidx = smoothG[a:b].argmax()
                maxG = max(smoothG[a:b])
                maxcoord = coords[a:b][maxidx]
                rdict[i].append(((coords[a],coords[b],maxcoord,maxG)))

    return rdict

class Peak:
    def __init__(self):
        self.regionleft = None
        self.regionright = None
        self.apex = None
        self.height = None
        self.left = None
        self.right = None
        self.chrom = None

    def __repr__(self):
        return str((self.apex, (self.left, self.right)) )

    def __str__(self):
        return __repr__(self)



def peaks(gstats, threshold, minrun=10, minlen=10000, perc=0.95):
    rdict = {}
    for i in range(len(gstats.chroms)):
        rdict[i] = []
        coords = gstats.coords[i]
        smoothG = gstats.smoothG[i]
        above = np.greater_equal(smoothG, threshold)
        runs = truth_runs(above, minrun)
        if len(runs):
            for (a,b) in runs:
                c1, c2 = coords[a], coords[b]
                if (c2-c1) < minlen:
                    continue
                maxidx = smoothG[a:b].argmax()
                maxG = max(smoothG[a:b])
                maxcoord = coords[a:b][maxidx]

                #right margin
                rightidx = maxidx
                for c in range(maxidx, len(smoothG[a:b])):
                    if smoothG[a:b][c] > (maxG*perc):
                        rightidx = c
                    else:
                        break
                rightmarg = coords[a:b][rightidx]

                #left margin
                leftidx = maxidx
                for c in range(maxidx, -1, -1):
                    if smoothG[a:b][c] > (maxG*perc):
                        leftidx = c
                    else:
                        break
                leftmarg = coords[a:b][leftidx]

                if leftmarg == rightmarg:
                    continue

                p = Peak()
                p.regionleft = coords[a]
                p.regionright = coords[b]
                p.apex = maxcoord
                p.height = maxG
                p.left = leftmarg
                p.right = rightmarg
                p.chrom = i+1
                p.snplocs = coords[a:b]

                rdict[i].append(p)
    return rdict


def view_sgdpeak(p):
    url = r"""http://fasolt.stanford.edu/cgi-bin/ORFMAP/ORFmap?chr=%d&beg=%d&end=%d"""
    webbrowser.open(url % (p.chrom, p.left, p.right))


def BH_fdr_procedure(pvals, q, independent=True):
    """ p-vals, FDR (q) -> p-value threshold to give desired false discovery rate

    Arguments:
        * pvals - seq of p-values for each test
        * q - desired False Discovery Rate (e.g. 0.05)
        * independent - specifies whether tests are independent

    if independent==True, use the standard Benjamini and Hochberg procedure.
    in independent==False, use the modified Benjamini and Yekutieli procedure.

    References:
    BENJAMINI Y, HOCHBERG Y
    CONTROLLING THE FALSE DISCOVERY RATE - A PRACTICAL AND POWERFUL APPROACH TO MULTIPLE TESTING
    J ROY STAT SOC B MET 57 (1): 289-300 1995

    Benjamini Y, Yekutieli D
    The control of the false discovery rate in multiple testing under dependency
    ANN STAT 29 (4): 1165-1188 AUG 2001

    """
    pvals = list(pvals)
    pvals.sort()
    if independent:
        jalpha = np.arange(1,len(pvals)+1) * (q/len(pvals))
    else:
        csum = sum(1.0/np.arange(1,len(pvals)+1))
        jalpha = np.arange(1,len(pvals)+1) * (q/(len(pvals)*csum))

    diff = pvals - jalpha
    neg = np.less(diff, 0)
    nz = [i for i in range(len(neg)) if neg[i] == True]
    if not len(nz):
        return None
    cutidx = max(nz)
    return pvals[cutidx]








def permute_counts(alleles):
    """ permute lists of allele counts as returned by load_allele_cts in bsacalc
    """
    allchroms = np.vstack([chrom[:,2:] for chrom in alleles])
    permuted = np.random.permutation(allchroms)
    newalleles = []
    idx = 0
    for chrom in alleles:
        n = len(chrom)
        newchrom = np.hstack([chrom[:,:2],permuted[idx:idx+n]])
        newalleles.append(newchrom)
        idx = idx + n
    return newalleles

def counts_to_frequencies(alleles):
    falleles = [i[:,2:].astype('f') for i in alleles]  # grab columns 3 and 4 which represent counts
    freqA0 = [i[:,0]/np.sum(i, axis=1) for i in falleles] # calculate freq of allele A0 
    return freqA0


def color_on_thresholds(vals, threshold, c1='r', c2='b'):
    clrs = np.where(vals > threshold, c1, c2)
    return list(clrs)


def Gstats_randomized_distributions(low, high, nperms=10, halfwidth=10000):
    G, smoothG = [],[]
    for i in range(nperms):
        lowperm = permute_counts(low)
        highperm = permute_counts(high)
        r = bsacalc.Gstats(lowperm, highperm, halfwidth=halfwidth)
        G += list(np.hstack(r.G))
        smoothG += list(np.hstack(r.smoothG))
    G.sort()
    smoothG.sort()
    return G, smoothG


def Gstats_randomized_distributions_onepool(low, p=0.5, nperms=10, halfwidth=10000):
    G, smoothG = [],[]
    for i in range(nperms):
        lowperm = permute_counts(low)
        r = bsacalc.Gstats_onepool(lowperm, p, halfwidth=halfwidth)
        G += list(np.hstack(r.G))
        smoothG += list(np.hstack(r.smoothG))
    G.sort()
    smoothG.sort()
    return G, smoothG




