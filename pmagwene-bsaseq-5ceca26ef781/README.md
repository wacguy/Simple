# bsaseq -- A python package for doing BSA-seq analyses.

IMPORTANT NOTE: As of September 2016, the dafault kernel smoothing approach has changed.  The default is now based on Savitsky-Golay filtering, which I have found to produce better results than the default tricube smoothing previously recommended in our 2011 paper.


This softare implements the analysis framework described in:

*  Magwene, P. M., J. H. Willis, and J. K. Kelly. 2011. The statistics of bulk segregant analysis using next-generation sequencing. PLoS Computational Biology, 7(11):e1002255. 


# Python Environment

The code in this repository requires packages from the standard "SciPy" stack -- numpy, matplotlib, etc.  If you are new to Python, the easiest way to get up and going quickly is to install the [Anaconda Scientific Python Distribution](https://store.continuum.io/cshop/anaconda/).

WINDOWS NOTE:  PDF output in Matplotlib 1.2.1 appears to be broken at Windows 64-bit (Win32 hasn't been tested). I will investigate but until further notice avoid the use of the PDF backend under Windows.


# bsacalc.py

`bsacalc.py` calculates raw and smoothed G-statistics on a per site basis.

## Required arguments 

- -L, -H : Files containing allele counts for one or more 'low' and one ore more 'high' bulks. Multiple files can be specified as so:

        bsacalc.py -L low1 low2 -H high1 high2 
    

### Structure of allele count files

The allele count files are tab-delimited text files with 4 columns, where columns are:

1. Chromosome name
2. Chromosomal coordinate of site
3. Counts of allele A0
4. Counts of allele A1    


## Optional arguments

- -o, --outfile: file for output
- -k, --kernel: type of smoothing kernel
- -w, --width: width of the smoothing kernel (in bp)


## Output

By default output is written to stdout. Output is tab-delimited with 4 columns, where columns are:

1. Chromosomal number as given by order of appearance in the input (0-indexed)
2. Chromosomal coordinates
3. Raw G statistics
4. Smoothed G-statistics



## Example Usage

From the examples directory:

    $ python ../bsacalc.py -L het1cts-simple.txt -H het3cts-complex.txt -w 33750 -k tricube -o ccm-example.out


# bsadraw.py

`bsadraw.py` draws figures illustrating the results of a BSA-seq analysis using output from `bsacalc.py`

## Required arguments

- -g, --gstats: G-statistics output from bsacalc.py. If -g is not specified that will look for input from stdin.
- -c, --chromlens: A file with a list of chromosome lengths (in bp). Order of chromosome should be same as order given in input to bsacalc.py.  Plain text file, one chromosome per line.
- -o, --outfile: The name of the file where the output figure is to be written to. The figure format is inferred from the filename extension. Supported extensions -- .pnf, .pdf, .jpg, .ps, .eps, .svg.

## Optional arguments

- -n: Chromosome number to draw (1-indexed). If n=0 (default) the output is the whole genome.
- --coords: start and stop coordinates to include in figure (1-indexed)
- --ylim: y-axis limits of drawing
- --figsize: width, and height of figure in inches
- --nticks: number of ticks on x-axis when drawing single chromosome
- --maxgap: maximum interval between sites for which curve lines are connected
- --noraw: Suppresses drawing of raw G-values 
- --rawclr: Color of points used to depict raw G-values
- --smoothclr: Color of lines used to depict smoothed G-values
- --smoothwidth: Width of lines used to depict smoothed G-values

## Output

File with format specified by `-o` option.

## Example usage

### Smoothed G-stats only, genome wide

    $ python ../bsadraw.py -g ccm-example.out -c scer_chromlen.txt -o genomic-smooth.png --noraw
    
### Smooth and raw G-stats, chromosome 10

    $ python ../bsadraw.py -g ccm-example.out -c scer_chromlen.txt -o chrom10.png --rawclr='0.5' -n 10 --figsize 5 4 --ylim 0 40 --smoothclr='green'

This changes the raw points to a medium gray (`--rawclr='0.5'`), the smoothed curve to green, and cuts of the y-axis at 40