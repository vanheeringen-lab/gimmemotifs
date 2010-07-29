/*
 * $Id: regress.c 67 2005-07-29 17:21:41Z nadya $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 17:25:42  nadya
 * Initial revision
 *
 */

#include "macros.h"

/*
	regress

	Least squares regression on points (x,y) to give
		y = mx + b

	Returns the root mean squared error of the fit.
*/
extern double regress(
  int n,			/* number of points */
  double *x,			/* x values */
  double *y,			/* y values */
  double *m,			/* slope */
  double *b 			/* y intercept */
)
{
  int i;
  double sx=0, sy=0, sxx=0, sxy=0;
  double mse=0;

  for (i=0; i<n; i++) {
    sx += x[i];
    sy += y[i];
    sxx += x[i]*x[i];
    sxy += x[i]*y[i];
  }

  *m = (n*sxy - sy*sx) / (n*sxx - sx*sx);
  *b = (sy - *m*sx)/n;

  for (i=0; i<n; i++) {
    double err = y[i] - (*m*x[i] + *b);
    mse += err * err;
  }
  mse = sqrt(mse);
  mse /= n;

  return mse;
}
