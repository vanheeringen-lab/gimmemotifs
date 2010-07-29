/* Convolve vectors using O(n log n) FFT algorithm */
/* Uses FFTW library */
#ifndef GLAM2_CONVOLVE_H
#define GLAM2_CONVOLVE_H

#include "glam2_fftw3.h"

typedef struct {
  int size;

  double *x;  /* 1st input */
  double *y;  /* 2nd input */
  double *z;  /* output */

  fftw_plan x_plan;
  fftw_plan y_plan;
  fftw_plan z_plan;
} fft_convolver;

/* Setup to convolve vectors of size <= max_size */
void fft_init(fft_convolver *f, const int max_size);

/* z = convolution of x and y */
void fft_convolve(double *z, const double *x, const double *y, const int size, fft_convolver *f);

#endif
