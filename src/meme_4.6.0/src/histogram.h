/*
 * $Id: histogram.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:40:04  nadya
 * Initial revision
 *
 */


#ifndef HISTOGRAM_H
#  define HISTOGRAM_H

#define N_BINS 20
#define N_TICKS 40
#define BAR_SYMBOL '*'
#define CHUNK 1000

/* create and print a histogram by reading from a file made by getscores */
extern void file_histogram(
  char *fname,          /* name of file; NULL means stdin */
  int n_bins,           /* number of bins (X-axis, plotted vertically) */
  int n_ticks,          /* number of ticks (Y-axis, plotted horizontally) */
  char bar_symbol       /* symbol to build bars from */
);

/* create and print a histogram from a file or array of values */
extern void histogram(
  FILE *file,           /* file containing values to plot */
  double *data_array,   /* array containing values to plot */
  int n,                /* size of array */
  double min,           /* minimum value on X-axis; ignored if array */
  double max,           /* maximum value on X-axis; ignored if array */
  int n_bins,           /* number of bins (X-axis, plotted vertically) */
  int n_ticks,          /* number of ticks (Y-axis, plotted horizontally) */
  char bar_symbol       /* symbol to build bars from */
);

/* set up to make a histogram */
extern void init_histogram(
  BOOLEAN stream,	/* reading from a file */
  double *data_array,	/* array containing values to plot */
  int n,		/* size of array */
  int *bins,		/* bins of histogram */
  int n_bins,		/* number of bins (X-axis, plotted vertically) */
  double *bin_size,	/* size of a bin */
  double *minimum,	/* minimum value on X-axis; ignored if array */
  double *maximum 	/* maximum value on X-axis; ignored if array */
);

/* fill up the bins */
extern void fill_histogram(
  FILE *file,		/* file containing values to plot */
  double *data_array,	/* array containing values to plot */
  int n,		/* size of array */
  double min,		/* minimum value on X-axis; ignored if array */
  double bin_size,	/* size of a bin */
  int *bins,		/* bins of histogram */
  int n_bins		/* number of bins (X-axis, plotted vertically) */
);

/* print the histogram set up with init_ and fill_ */ 
extern void finish_histogram(
  double min,		/* minimum value on X-axis; ignored if array */
  int n_bins,		/* number of bins (X-axis, plotted vertically) */
  int *bins,		/* bins of histogram */
  double bin_size,	/* size of a bin */
  int n_ticks,		/* number of ticks (Y-axis, plotted horizontally) */
  char bar_symbol,	/* symbol to build bars from */
  double *avg,          /* mean of data */
  double *sd            /* standard deviation of data */
);
extern void finish_blastogram(
  double min,		/* minimum value on X-axis; ignored if array */
  int n_bins,		/* number of bins (X-axis, plotted vertically) */
  int *bins,		/* bins of histogram */
  double bin_size,	/* size of a bin */
  int n_ticks,		/* number of ticks (Y-axis, plotted horizontally) */
  char bar_symbol,	/* symbol to build bars from */
  double *avg,		/* mean of data */
  double *sd		/* standard deviation of data */
);

#endif

