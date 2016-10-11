import sys
import math, csv
from itertools import cycle

import numpy as np
from matplotlib import pylab
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter, MaxNLocator



def kilobases(x,pos):
    return "%1.0fKb" % (x * 1e-3)    
kbaseformatter = FuncFormatter(kilobases)  


def chrom_offset(chromlens):
    return [sum(chromlens[:i]) for i in range(len(chromlens))]

    
def connected_intervals(x, y, maxgap=10000):
    """Determine the connected intervals over a set of x,y observations by 
    identifying the 'gaps'.
    
    Gaps are defined as adjacent x values where x[i+1]-x[i] > maxgap
    
    This is useful when you want to draw a plot over a set of data but you 
    don't want to connect points that span intervals where there is no data.
    """
    x = np.array(x)
    y = np.array(y)
    diff = x[1:] - x[:-1]
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
            
    

def setup_chromfig(chromlen, ymax=None, fig=None, ax=None):
    """Configure a matplotlib figure object for drawing a single chromosome.
    """    
    
    if fig is None:
        fig = pyplot.figure(figsize=(8,4),dpi=150)

    if ax is None:
        ax = fig.add_subplot(111) 
    
    ax.set_xlim(0,chromlen)
    if ymax is not None:
        ax.set_ylim(0, ymax)
        
    ax.spines['left'].set_position(('outward',10))
    ax.spines['bottom'].set_position(('outward',10))    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.yaxis.set_ticks_position('left')      
    ax.xaxis.set_ticks_position('bottom')      
    ax.set_ylabel('G-statistic')
    return fig 
    
def setup_many_chromfig(chromlens, ymax, chrbounds=True, 
        boundstyle='dotted', chrfontsize=8, fig=None, ax=None, labelchroms=True,
        chromlabeloffset=0.025):
    """Configure a matplotlib figure object for drawing multiple chromosomes end to end.
    """
    runsum = [sum(chromlens[:i]) for i in range(len(chromlens))]
    nchroms = len(chromlens)
    xmax = (runsum[-1] + chromlens[-1])
    xborder = xmax * 0.01

    if fig is None:
        fig = pyplot.figure(figsize=(10,6),dpi=150)

    if ax is None:
        ax = fig.add_subplot(111)

    chrlabelY = -(round(ymax * chromlabeloffset))
        
    if chrbounds:
        ax.vlines(runsum,-2,ymax,linestyles=boundstyle,linewidth=0.25,color='0.2')    
        for i in range(len(runsum)):
            xpt = runsum[i]+(chromlens[i]/2.0)
            if labelchroms:
                ax.text(xpt, chrlabelY, str(i+1),fontsize=chrfontsize)

    ax.set_ylim(0, ymax)
    ax.set_xlim(0, xmax)   
    ax.set_xticks([]) 
    ax.spines['left'].set_position(('outward',10))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    ax.set_ylabel('G-statistic')    
    return fig
            
    

def plot_chrom(x, y, maxgap=10000, color='r', ax=None, **kw):
    """ Plot a connected curve over a set of chromosomal coordinates.
    """
    newx,newy = connected_intervals(x, y, maxgap)
    for j in range(len(newx)):
        nx = newx[j]
        ny = newy[j]
        if ax is None:
            pyplot.plot(nx, ny, color=color, **kw)  
        else:
            ax.plot(nx, ny, color=color, **kw)


def fillbetween_chrom(x, y1, y2=0, maxgap=10000, color='r', **kw):
    """ Plot a filled curve over a set of chromosomal coordinates.
    """
    newx,newy1 = connected_intervals(x, y1, maxgap)
    newy2 = [0] * len(newx)
    if y2 is not 0:
        newx, newy2 = connected_intervals(x, y2, maxgap)
    for j in range(len(newx)):
        nx = newx[j]
        ny1 = newy1[j]
        ny2 = newy2[j]
        pyplot.fill_between(nx, ny1, ny2, color=color, **kw)


        
def as_colorcycle(c):
    if isinstance(c, str):
        return cycle([c])
    else:
        return cycle(list(c))
        

def plot_many_chroms(X, Y, chromlens, maxgap=10000, color='r', **kw):
    """ Given multiple chromosomal coords and y-values, plot connected curves
    over each chromosome.
    
    Chromosomes are drawn end to end.
    """
    offset = [sum(chromlens[:i]) for i in range(len(chromlens))]
    nchroms = len(X)
    clrs = as_colorcycle(color)    
    for i in range(nchroms):
        x = X[i] + offset[i]
        y = Y[i]
        plot_chrom(x, y, maxgap, color=clrs.next(), **kw)

def fillbetween_many_chroms(X, Y,  chromlens, maxgap=10000, color='r', **kw):
    offset = [sum(chromlens[:i]) for i in range(len(chromlens))]
    nchroms = len(X)
    clrs = as_colorcycle(color)    
    for i in range(nchroms):
        x = X[i] + offset[i]
        y = Y[i]
        fillbetween_chrom(x, y, y2=0, maxgap=maxgap, color=clrs.next(), **kw)    


def plot_raw(x,y,color='k',size=1, marker='.', ax=None, **kw):
    """ Draw y-values over chromosomal coordinates as points.
    """
    if ax is None:
        pyplot.plot(x, y, marker=marker, ls='None', markersize=size, color=color, **kw)
    else:
        ax.plot(x, y, marker=marker, ls='None', markersize=size, color=color, **kw)


def plot_many_raw(X, Y, chromlens, color='k', size=1, marker='.', **kw):
    """ Draw y-values over chromosomal coordinates as points for mulitple chromosomes.
    """
    offset = [sum(chromlens[:i]) for i in range(len(chromlens))]
    nchroms = len(X)
    clrs = as_colorcycle(color)  
    for i in range(nchroms):
        x = X[i] + offset[i]
        y = Y[i]
        plot_raw(x, y, color=clrs.next(), size=size, marker=marker, **kw)    


def draw_region(start, end,  height, **kwargs):
    """ Draw a rectangular patch that spans coordinates start to end, with given height.
    """
    r = pylab.Rectangle((start, 0), abs(start-end), height, **kwargs)
    pylab.gca().add_patch(r)
    pylab.draw()

def draw_ftr(ftr,  height, **kwargs):
    """ Convenience function for draw_region, for objects with
    start and end attributes.
    """
    draw_region(ftr.start, ftr.end, height, **kwargs)



def plotG_one(i, gresults, chromlens, drawraw=True,
        rawcolor='k',  rawsize=1, rawmarker='.', rawalpha=0.5,
        smoothcolor='r', smoothwidth=1.5, smoothalpha=1, 
        ymax=None, maxgap=10000, fig=None):  
    """ A convenience function for drawing G-stats over a single chromosome.
    """
    if ymax is None:
        if drawraw:
            ymax = max(gresults.G[i]) * 1.1
        else:
            ymax = max(gresults.smoothG[i]) * 1.1    
    x = gresults.coords[i]
    if fig is None:
        fig = setup_chromfig(chromlens[i], ymax)
    if drawraw:
        plot_raw(x, gresults.G[i], color=rawcolor, size=rawsize, marker=rawmarker, alpha=rawalpha)
    plot_chrom(x, gresults.smoothG[i], color=smoothcolor, linewidth=smoothwidth, 
        maxgap=maxgap, alpha=smoothalpha)
    return fig
    


def plotG_all(gresults, chromlens,  drawraw=True,
        rawcolor='k',  rawsize=1, rawmarker='.', rawalpha=0.25,
        smoothcolor='r', smoothwidth=1.0,  smoothalpha=1,
        chrbounds=True, boundstyle='dotted', chrfontsize=8,
        ymax=None, maxgap=10000, fig=None):  
    """ A convenience function for drawing G-stats over a genome.
    """
    if ymax is None:
        if drawraw:
            ymax = max([max(i) for i in gresults.G]) * 1.1
        else:
            ymax = max([max(i) for i in gresults.smoothG]) * 1.1
    
    if fig is None:
        fig = setup_many_chromfig(chromlens, ymax, 
            chrbounds=chrbounds, boundstyle=boundstyle, chrfontsize=chrfontsize)
    X = gresults.coords
    if drawraw:
        plot_many_raw(X, gresults.G, chromlens, 
                    color=rawcolor, size=rawsize, marker=rawmarker, alpha=rawalpha)
        
    plot_many_chroms(X, gresults.smoothG, chromlens,
                    color=smoothcolor, linewidth=smoothwidth,
                    maxgap=maxgap, alpha=smoothalpha)
    return fig


class Results:
    pass
    
def load_Gstats(fname, withmaxmin=False):
    """ Load Gstats output file from bsacalc.py
    """
    if withmaxmin:
        chroms, coords, G, smoothG, maxminG = np.loadtxt(fname, dtype=('i8,i8,f8,f8,f8'), unpack=True)
    else:
        chroms, coords, G, smoothG = np.loadtxt(fname, dtype=('i8,i8,f8,f8'), unpack=True)        

    chromi = sorted(list(set(chroms)))
    r = Results()
    r.chroms = [chroms[chroms == i] for i in chromi]
    r.coords = [coords[chroms == i] for i in chromi]
    r.G = [G[chroms == i] for i in chromi]
    r.smoothG = [smoothG[chroms == i] for i in chromi]
    if withmaxmin:
        r.maxminG = [maxminG[chroms == i] for i in chromi]
    return r

    
def subset_Gstats(gresults, totake):
    """ Subset a Gstats result object from the chromosomes indicated in totake.
    """
    r = Results()
    r.chroms = [gresults.chroms[i] for i in totake]
    r.coords = [gresults.coords[i] for i in totake]
    r.G = [gresults.G[i] for i in totake]
    r.smoothG = [gresults.smoothG[i] for i in totake]
    return r
    


if __name__ == "__main__":
    import os.path, argparse
    parser = argparse.ArgumentParser(description="Draw BSA mapping from high-throughput sequencing")
 
    parser.add_argument('-g','--gstats', nargs='?', help="G stats, output of bsacalc.py",
                 type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-c','--chromlens', required=True, 
                 type=argparse.FileType('r'), help="File giving chrom lens")
    parser.add_argument('-o', '--outfile',required=True,
                 type=str, help="Output graphics file")
    parser.add_argument('-n', type=int, default=0, dest='chromn',
                help="Specify chromosome number to draw (1-indexed)")
    parser.add_argument('--coords', nargs=2, type=int, dest='coords',
                default=argparse.SUPPRESS, required=False,
                help="Specify coordinate range to draw (1-indexed)")
    parser.add_argument('--ylim', nargs=2, type=int, dest='ylim',
                default=argparse.SUPPRESS, required=False,
                help="Specify y-axis limits for figure")
    parser.add_argument('--figsize', nargs=2, type=float, dest='figsize',
                default=argparse.SUPPRESS, required=False,
                help="Specify width, height of figure in inches")
    parser.add_argument('--nticks', nargs=1, type=int, dest='nticks',
                default=argparse.SUPPRESS, required=False,
                help="Specify number of ticks on x-axis when drawing single chromosome") 
    parser.add_argument('--maxgap', nargs=1, type=int, dest='maxgap',
                default=10000, required=False,
                help="Specify maximum interval between connected sites for curve drawing")   

    parser.add_argument('--maxmin', action='store_true', dest='maxmin', default=False,
                help="Include the moving maxmin values of the G statistic in the output")                                                               
    
    parser.add_argument('--threshold', nargs=1, type=float, dest='threshold',
                default=0.0, required=False,
                help="Specify value along the ordinate at which to draw significance threshold") 

    parser.add_argument('--noraw', action='store_false', dest='drawraw',
                help="Don't draw raw G values")
    parser.add_argument('--rawclr',type=str, default='gray')
    parser.add_argument('--smoothclr',type=str, default='r')
    parser.add_argument('--thresholdclr',type=str, default='k')
    parser.add_argument('--smoothwidth',type=float,default=1)
    args = parser.parse_args()
    
    r = load_Gstats(args.gstats, withmaxmin=args.maxmin)
    clens = [int(i) for i in args.chromlens.readlines()] 
    tlen = sum(clens)

    smoothclr = args.smoothclr
    rawclr = args.rawclr
    smoothwidth = args.smoothwidth

    n = args.chromn
    if n <= 0:
        fig = plotG_all(r, clens, drawraw=args.drawraw,rawcolor=rawclr,
                        smoothcolor=smoothclr, smoothwidth=smoothwidth,
                        maxgap=args.maxgap) 
        if args.maxmin:
            fillbetween_many_chroms(r.coords, r.maxminG, clens, maxgap=args.maxgap, alpha=0.5)

    else:
        fig = plotG_one(n-1, r, clens, drawraw=args.drawraw,rawcolor=rawclr,
                        smoothcolor=smoothclr, smoothwidth=smoothwidth,
                        maxgap=args.maxgap) 
        if args.maxmin:
            fillbetween_chrom(r.coords[n-1], r.maxminG[n-1], maxgap=args.maxgap, alpha=0.5)

    ax = fig.gca()
    if 'coords' in vars(args):
        ax.set_xlim(*args.coords)

    if n > 0:
        if 'nticks' in vars(args):
            m = MaxNLocator(*args.nticks)
            m.view_limits(*ax.get_xlim())
            ax.xaxis.set_major_locator(m)
        else:
            ax.xaxis.set_major_locator(MaxNLocator(4))

    if 'ylim' in vars(args):
        ax.set_ylim(*args.ylim)
    if 'figsize' in vars(args):
        fig.set_size_inches(*args.figsize)

    if args.threshold > 0:
        ax.hlines(args.threshold, 0, tlen,  color=args.thresholdclr, linestyle='dashed')
    
    name,ext = os.path.splitext(args.outfile)
    # check extension of outfile
    if ext in ('.png','.pdf','.jpg','.ps','.eps','.svg'):
        ftype = ext[1:]
    else:   # if extension not recognized output png as default
        ftype = 'png'
    fig.savefig(args.outfile, format=ftype)
    
                
