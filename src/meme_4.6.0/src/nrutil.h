/*
 * $Id: nrutil.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:46:17  nadya
 * Initial revision
 *
 */

#ifndef NRUTIL_H
#define NRUTIL_H
void nrerror(char *error_text);
float *vector(int nl, int nh);
int *ivector(int nl, int nh);
double *dvector(int nl, int nh);
float **matrix(int nrl, int nrh, int ncl, int nch);
float **matrix(int nrl, int nrh, int ncl, int nch);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float **submatrix(float a, int oldrl, int oldrh, int oldcl,
  int oldch, int newrl, int newcl);
void free_vector(float *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_submatrix(float **b, int nrl, int nrh, int ncl, int nch);
float **convert_matrix(float *a, int nrl, int nrh, int ncl, int nch);
void free_convert_matrix(float **b, int nrl, int nrh, int ncl, int nch);
#endif
