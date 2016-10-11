import math
import numpy as np
import scipy as sp
import scipy.stats as stats
import scipy.signal as signal




def map2recomb(d, func="haldane"):
    """ map distance (in Morgans) -> recombination prob.

    Possible options to func are "haldane" and "kosambi". Defaults
    to haldane.
    """
    d = np.abs(d)
    if func == "kosambi":
        r = 0.5* (1 - np.exp(-4.0 * d))/(1 + np.exp(-4.0 * d))
    else:  
        r = 0.5 * (1.0 - np.exp(-2.0 * d)) # default to haldane's function
    return r


def simulate_recombinant_chromosomes(nchroms, nsites=1000, interval=0.01):
    """ simulate the generation of recombinant chromosomes from a cross btw inbred strains
    
    # chromosome, # sites, interval between sites in Morgans -> simulated recombinant chromosomes.

    Return an array of dim (nchroms, nsites), so that rows represent chromosomes and
    columns represent sites.

    The two parental alleles are designated -1 (ref) and 1 (alt).
    """
    r = map2recomb(interval)

    # poisson sampling to pick intervals where crossovers occur
    xover = []
    while len(xover) < nchroms:
        chrom = stats.poisson.rvs(r, size=nsites)
        if np.sum(chrom) >= 1:  # require at least one crossover
            xover.append(chrom)
    xover = np.array(xover)

    xover[xover % 2 == 1] = 1  # odd number of crossovers are the only ones that change state
    xover = -xover             # designate intervals where crossovers occur with value -1
    xover[xover % 2 == 0] = 1  # even number of crossovers make for no change of allelic state

    # randomly assign an initial parental state at the first site
    startg = np.random.randint(2, size=(nchroms,)) 
    startg[startg == 0] = -1
    xover[:,0] = startg

    # this is a cute trick. The elems of the rows of xprod are either -1 or 1. 
    # -1 indicates a cross-over occurs, 1 indicates no (effective) crossover occurs
    # given initial alleles in column 0, the cumulative product (columnwise) gives
    # us the corresponding genotypes generated from those crossover events
    return np.cumprod(xover,axis=1)


def to_01genotypes(genotypes):
    """ Two value numeric encoding of genotype -> to (0,1) encoding

    Assumes minimum value indicates ancestral allele which is coded 0.  Any other
    values is coded as 1.
    """
    genotypes = np.asarray(genotypes)
    ref = genotypes.min()
    G = np.where(genotypes == ref, 
                 np.zeros_like(genotypes), 
                 np.ones_like(genotypes))
    return G


def chrompair_to_diploid_genotypes(c1, c2):
    """Pair of arrays representing genotypes on each chrom -> (-1,0,1) diploid genotypes
    """
    C1 = to_01genotypes(c1)
    C2 = to_01genotypes(c2)
    return (C1+C2) - 1



def haploid_phenotypes(genotypes, effectsize):
    """effect size, genotypes -> simulated phenotype (column vector)

    genotypes = matrix of biallelic genotypes at causal sites. Each row corresponds to an individual,
                 columns to sites. Genotypes are coded as -1 = ancestral, 1 = derived

    effectsizes = floating point number in range (0,1) representing effect size of causal variant

    """
    pi = effectsize
    g = to_01genotypes(genotypes)
    n = len(g)
    p = np.sum(g)/float(n)  # frequency of causal polymorphisms
    q = 1.0 - p
    wtPi = pi / (p * q)
    genetic_effect = g * math.sqrt(wtPi)
    random_effect = math.sqrt(1.0 - pi) * np.random.normal(0,1,size=n)
    return genetic_effect + random_effect


def diploid_phenotypes(genotypes, effectsize):
    """effect size, genotypes -> simulated phenotype (column vector)

    genotypes = matrix of biallelic genotypes at causal sites. Each row corresponds to an individual,
                 columns to sites. Genotypes are coded as -1 = ancestral, 1 = derived

    effectsizes = floating point number in range (0,1) representing effect size of causal variant

    """
    pi = effectsize
    g = np.asarray(genotypes)
    n = len(g)
    p = (np.sum(g == 1) + 0.5*np.sum(g == 0))/float(n)  # frequency of causal polymorphisms
    q = 1.0 - p
    wtPi = pi / (2 * p * q)
    genetic_effect = g * math.sqrt(wtPi)
    random_effect = math.sqrt(1.0 - pi) * np.random.normal(0,1,size=n)
    return genetic_effect + random_effect


def genotypes_to_ref_freqs(genotypes):
    return (genotypes == -1).sum(axis=0)/float(len(genotypes))


def simulate_BSA_bulking(phenos, genos, bulkpct=20):
    phenos = np.asarray(phenos)
    pbulks, gbulks = [], []
    lowsegs = phenos <= np.percentile(phenos, bulkpct)
    highsegs = phenos >= np.percentile(phenos, 100 - bulkpct)
    for bulk in (lowsegs, highsegs):
        pbulks.append(phenos[bulk])
        gbulks.append(genos[bulk])
    return pbulks, gbulks
 

def simulate_BSA_sequencing(C, ref_freqs, model='poisson'):
    if model == 'nbinom':
        ctref = stats.nbinom.rvs(C*ref_freqs, 0.5)
        ctalt = stats.poisson.rvs(C*(1.0-ref_freqs), 0.5)
    else:
        ctref = stats.poisson.rvs(C*ref_freqs)
        ctalt = stats.poisson.rvs(C*(1.0-ref_freqs))
    return ctref, ctalt


def simulate_BSA_experiment(popsize, nsites = 1001, interval = 0.01,
            causalsite=500, effectsize=0, bulkpct=20, diploid=False, 
            coverage=100, coverage_model="poisson"):
    """ Simulate a bulk segregant analysis experiment based on F2 offspring.
    """
    if diploid:
        c1 = simulate_recombinant_chromosomes(popsize, nsites, interval)
        c2 = simulate_recombinant_chromosomes(popsize, nsites, interval)
        genos = chrompair_to_diploid_genotypes(c1, c2)
        phenos = diploid_phenotypes(genos[:,causalsite], effectsize)
    else:
        genos = simulate_recombinant_chromosomes(popsize, nsites, interval)
        phenos = haploid_phenotypes(genos[:,causalsite], effectsize)
    pbulks, gbulks = simulate_BSA_bulking(phenos, genos, bulkpct)
    freqbulks = [genotypes_to_ref_freqs(i) for i in gbulks]
    cvgbulks = [np.column_stack(simulate_BSA_sequencing(coverage, f, coverage_model)) for f in freqbulks]

    class Results:
        pass

    r = Results
    r.genotypes = genos
    r.causal_genotypes = genos[:,causalsite]
    r.phenotypes = phenos
    r.phenotypic_bulks = pbulks
    r.genotypic_bulks = gbulks
    r.freq_bulks = freqbulks
    r.observed_counts = cvgbulks
    return r

def Gstats_from_BSA_simulation(bsaresults, halfwidth=50, sgdegree=6):
    lowbulk = bsaresults.observed_counts[0]
    highbulk = bsaresults.observed_counts[1]
    n1, n2 = lowbulk[:,0], lowbulk[:,1]
    n3, n4 = highbulk[:,0], highbulk[:,1]
    G = Gtest_indep(n1, n2, n3, n4)
    smoothG = signal.savgol_filter(G, 2*halfwidth+1, sgdegree, mode='mirror')

    bsaresults.G = G
    bsaresults.smoothG = smoothG
    return bsaresults


def Gtest_indep(a, b, c, d, const=0.01):
    """Calculates G-test for independence of 2 x 2 tables.
    
    see Sokal and Rohlf 1994, eqn 17.11 and Box 17.6 p. 731
    """
    a = np.array(a,dtype=np.float) + const
    b = np.array(b,dtype=np.float) + const
    c = np.array(c,dtype=np.float) + const
    d = np.array(d,dtype=np.float) + const
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


def extreme_value_cutoff_from_BSA_simulations(bsasims, percentile=95):
    M = [max(i.smoothG) for i in bsasims]
    fxn = stats.frechet_r
    params = fxn.fit(M)
    f = fxn(*params)
    simcutoff = np.percentile(M, percentile)
    fitcutoff = f.ppf(percentile/100.)

    class Results:
        pass
    
    r = Results()
    r.simulation_cutoff = simcutoff
    r.evt_cutoff = fitcutoff
    r.maxima = M
    r.fit_distn = f
    return r


