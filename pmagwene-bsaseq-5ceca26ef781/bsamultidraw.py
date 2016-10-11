
import sys
import math, csv
from itertools import cycle

import numpy as np
from matplotlib import pylab
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter, MaxNLocator

import bsadraw






if __name__ == "__main__":
    import os.path, argparse
    parser = argparse.ArgumentParser(description="Draw BSA figure comparing multiple analysis files on same genome.")
 
    parser.add_argument('-g','--gstats', nargs='+', help="G stats, output of bsacalc.py",
                 type=str, required=True)
    
    parser.add_argument('-c','--chromlens', required=True, 
                 type=argparse.FileType('r'), help="File giving chrom lens")
    
    parser.add_argument('-o', '--outfile',required=True,
                 type=argparse.FileType('w'), help="Output graphics file")

    parser.add_argument('--colors',type=str, default='r', nargs='+')
    parser.add_argument('--smoothwidth',type=float,default=1.5)    
    parser.add_argument('--smoothalpha',type=float,default=0.75)


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
    
    parser.add_argument('--threshold', nargs=1, type=float, dest='threshold',
                default=0.0, required=False,
                help="Specify value along the ordinate at which to draw significance threshold") 

    
    parser.add_argument('--thresholdclr',type=str, default='k')

    args = parser.parse_args()

    clens = [int(i) for i in args.chromlens.readlines()] 
    n = args.chromn
    tlen = sum(clens)
    smoothclr = bsadraw.as_colorcycle(args.colors)
    smoothwidth = args.smoothwidth
    smoothalpha = args.smoothalpha


    rs = []
    maxys = []
    for fname in args.gstats:
        r = bsadraw.load_Gstats(fname)
        rs.append(r)
        rmax = max([max(i) for i in r.smoothG])
        maxys.append(rmax)

    ymax = max(maxys) * 1.1

    if n <= 0:
        fig = bsadraw.setup_many_chromfig(clens, ymax)
    else:
        fig = bsadraw.setup_chromfig(clens[n-1], ymax)

    for r in rs:
        if n <= 0:
            bsadraw.plotG_all(r, clens, fig=fig, drawraw=False, 
                         smoothcolor=smoothclr.next(), smoothwidth=smoothwidth,
                         maxgap=args.maxgap,smoothalpha=smoothalpha) 
        else:
            bsadraw.plotG_one(n-1, r, clens, fig=fig,drawraw=False,
                            smoothcolor=smoothclr.next(), smoothwidth=smoothwidth,
                            maxgap=args.maxgap) 

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

    name,ext = os.path.splitext(args.outfile.name)
    # check extension of outfile
    if ext in ('.png','.pdf','.jpg','.ps','.eps','.svg'):
        ftype = ext[1:]
    else:   # if extension not recognized output png as default
        ftype = 'png'
    fig.savefig(args.outfile, format=ftype)

