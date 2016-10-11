# -*- coding: utf-8 -*-

import math

import numpy as np
from scipy import stats, signal


def uniform(x):
    """Uniform(square) kernel."""
    z = np.ones(len(x),np.float)
    z[np.abs(x) >= 1] = 0
    return z         
    
def triangle(x):
    """Triangular kernel."""
    z = (1.0 - np.abs(x))
    z[np.abs(x) >= 1] = 0
    return z        

def tricube(x):
    """Tricube kernel."""
    z = (1.0 - np.abs(x)**3)**3
    z[np.abs(x) >= 1] = 0
    return z
    
def cubed(x):
    """Cubic kernel."""
    z = (1.0 - np.abs(x)**3)
    z[np.abs(x) >= 1] = 0
    return z        

def triweight(x):
    """Triweight kernel."""
    z = 1.09375 * ((1.0 - np.array(x)**2)**3)
    z[np.abs(x) >= 1] = 0
    return z
    
def epanechnikov(x):
    """Epanechnikov kernel."""
    z = 0.75 * (1.0 - np.array(x)**2)
    z[np.abs(x) >= 1] = 0
    return z   
        
def quartic(x):
    """Quartic kernel."""
    z = 0.9375 * ((1.0 - np.array(x)**2)**2)
    z[np.abs(x) >= 1] = 0
    return z  
    
biweight = quartic 


rsqrt2pi = 1./math.sqrt(2.0*math.pi)

def gaussian(x):
    """Gaussian kernel."""
    z = rsqrt2pi * np.exp(-0.5 * x**2)
    return z    
    
    
pi4 = math.pi/4.
pi2 = math.pi/2.

def cosine(x):
    """Cosine kernel."""
    z = pi4 * np.cos(pi2*x)
    z[np.abs(x) >= 1] = 0
    return z    


def reflect_ends(y, h):
    ystart = y[:h]
    yend = y[-h:]
    Y = np.hstack((ystart[::-1], y, yend[::-1]))
    return Y


def fast_smooth(y, h, kernel=epanechnikov, kinterval=(-1,1)):
    """Uses convolution to smooth y, using a kernel of half-bandwidth h, where h is # of points.

    Like kernel_smooth but assumes uniform interpoint distances so can
    use numpy.convolve for fast results.
    """
    wts = kernel(np.linspace(kinterval[0], kinterval[1], 2*h + 1))
    sumwts = np.sum(wts)
    Y = reflect_ends(y, h)
    smoothY = np.convolve(Y, wts/sumwts, mode='valid')
    return smoothY

def fast_kernel_smooth(x, y, h, kernel=epanechnikov, maxgap=10000):
    """ A faster version of kernel smoothing.

    Here the half-width h is defined in terms of the number of points in the 
    filter window.  Thus we're effectively assuming approximately uniform sampling
    over the window.  To avoid "jumps" in smoothing, you can use the maxgap parameter 
    to specify the maximum distant in x, between adjacent points, for which you want
    to treat the points as neighbors.  This allows you to to break up the signal into subparts, 
    each of which is filtered independently.        
    """
    x = np.asarray(x)
    y = np.asarray(y)

    xparts, yparts = connected_intervals(x,y,maxgap)
    smoothedparts = []
    for part in yparts:
        npart = len(part)
        if (npart % 2) == 0:  # even
            npart -= 1
        width = min(npart,h)
        spart = fast_smooth(part, width, kernel)
        smoothedparts.append(spart)
    return np.concatenate(smoothedparts)


def kernel_smooth(x,y,h, kernel=uniform, ctinband=False):
    """Calculate a `kernel smoothed' moving average of y, at the given
    x-coords and with half-width, h. Assumes the x's are already in sorted order. 
    The half-width here is defined in terms of values of x NOT number of points, so
    the number of points in the smoothing window varies over the range of x.
    
    This is a function to calculate moving averages given uneven sampling. 
    Calculates moving average of y at coordinate x_i by weighted averaging 
    over all points in range (x_i-h, x_i+h). Weights are given by the kernel
    function. This is equivalent to the Nadaraya-Watson regression estimate.
    
    To deal with beginning/end of intervals, reflects left/right half-bandwiths
    worth of data. These reflected half bandwidths are trimmed off before
    the smoothed data is returned.
    
    If ctinband==True will return the number of points in each smoothing window.
    
    """
    x = np.array(x)
    y = np.array(y)
    
    olen = len(x)
    
    # at the head and tail of the sequence reflect the right and left (respectively)
    # half-bandwidths
    xfirst = x[0]
    startband = x <= xfirst + h
    xstart = xfirst - (x[startband] - xfirst)
    ystart = y[startband]
    nstart = len(xstart)-1
    
    xlast = x[-1]
    endband = x >= xlast - h
    xend = xlast + (xlast - x[endband])
    yend = y[endband]
    nend = len(xend) - 1    
    x = np.hstack((xstart[::-1][:-1], x, xend[::-1][1:]))
    y = np.hstack((ystart[::-1][:-1], y, yend[::-1][1:]))

    
    lx, ly = len(x), len(y)
    if lx != ly:
        raise Exception("x and y must be same length.")
    z = []
    ninband = []
    for i in range(lx):
        c = x[i]
        inband =  np.logical_and(x >= (c-h), x <= (c+h))
        if ctinband:
            ninband.append(len(np.flatnonzero(inband)))
        xfrac = (np.abs(x[inband] - c))/float(h)
        xwt = kernel(xfrac)
        ywin = np.sum(y[inband]*xwt)/np.sum(xwt)
        z.append(ywin)
    
    # trim off the reflected right/left half-bandwidths
    if ctinband:
        return z[nstart:olen+nstart], ninband[nstart:olen+nstart]        
    return z[nstart:olen+nstart]   





def connected_intervals(x, y, maxgap=10000):
    """Determine the connected intervals over a set of x,y observations by 
    identifying the 'gaps'.
    
    Gaps are defined as adjacent x values where x[i+1]-x[i] > maxgap
    
    This is useful when you want to draw a plot over a set of data but you 
    don't want to connect points that span intervals where there is no data.
    """
    x = np.array(x)
    y = np.array(y)
    diff = np.diff(x)
    gaps = [i+1 for i in np.flatnonzero(diff > maxgap)]
    if not len(gaps):
        return [x],[y]
    newx, newy = [], []
    idx = 0
    for i,j in enumerate(gaps):
        jx = x[idx:j]
        jy = y[idx:j]
        newx.append(jx)
        newy.append(jy)
        if i == len(gaps)-1:
            newx.append(x[j:])
            newy.append(y[j:])
        idx = j
    return newx, newy
 

def savgol_smooth(x,y, h, maxgap=10000, degree=6, mode='mirror'):
    """Apply Savitsky-Golay filtering of a signal y, over the domain of x.

    Here the half-width h is defined in terms of the number of points in the 
    filter window.  Thus we're effectively assuming approximately uniform sampling
    over the window.  To avoid "jumps" in smoothing, you can use the maxgap parameter 
    to specify the maximum distant in x, between adjacent points, for which you want
    to treat the points as neighbors.  This allows you to to break up the signal into subparts, 
    each of which is filtered independently.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    xparts, yparts = connected_intervals(x,y,maxgap)
    smoothedparts = []
    for part in yparts:
        npart = len(part)
        if npart == 1:
            smoothedparts.append(part)
            continue
        if (npart % 2) == 0:  # even
            npart -= 1
        width = min(npart,2*h+1)
        polydeg = min(width-1, degree)
        spart = signal.savgol_filter(part, width, polydeg, mode=mode)
        smoothedparts.append(spart)
    return np.concatenate(smoothedparts)

def minfilter(x, size=5):
    return signal.order_filter(x, np.ones(size), 0)

def maxfilter(x, size=5):
    return signal.order_filter(x, np.ones(size), size-1)
    #return signal.order_filter(x, np.ones(size), size-1)

def minmaxfilter(x, size=5):
    return minfilter(maxfilter(x, size), size)

def maxminfilter(x, size=5):
    return maxfilter(minfilter(x, size), size)

def connected_maxminfilter(x, y, size=5, maxgap=10000):
    gapx, gapy = connected_intervals(x, y, maxgap)
    maxminy = [maxminfilter(i, size) for i in gapy]
    return np.concatenate(gapx), np.concatenate(maxminy)
