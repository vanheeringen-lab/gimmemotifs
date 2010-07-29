/* MoAn is a motif discovery tool for DNA sequences */
/* Copyright (C) 2006 Eivind Valen */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

#ifndef _PSSMFUNC_H_
#define _PSSMFUNC_H_

#define WM(ind)         search->wm[2 * ind + search->cur[ind]]    /* The selected WM */
#define NWM(ind)        search->wm[2 * ind + !search->cur[ind]]   /* The unselected WM  */
#define SWITCH_WM(ind)  search->cur[ind] = !search->cur[ind]      /* Switch between the two WMs */


PSSM load_log_matrices(char *filename, int alphlen);
PSSM load_count_matrices(char *filename, int alphlen);
void print_pssm(PSSM pssm);
PSSM *random_pssms(unsigned int count, const unsigned char order,  const unsigned char min_length, const unsigned char max_length, const unsigned char alphlen, const int seq);
PSSM  *random_pssm(unsigned char order,  const unsigned char min_length, const unsigned char max_length, unsigned char alphlen, int seq);

/* inline void normalize_col(PSSM pssm, unsigned int pos, unsigned int sum); */
/* inline void normalize_all(PSSM pssm, unsigned int sum); */


#endif
