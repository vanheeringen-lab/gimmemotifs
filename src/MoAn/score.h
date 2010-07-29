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

#ifndef _SCOREL_H_
#define _SCORE_H_

#include "anneal.h"


/* Scoring functions */
float score_cooc_likelihood(Annealing *search);
float score_single_likelihood(Annealing *search);
float score_likelihood_lennorm(Annealing *search);
float score_single_aic(Annealing *search);
float score_single_bayesian(Annealing *search);
float score_bayesian(Annealing *search);
float score_single_aic(Annealing *search);
float score_aic(Annealing *search);


#endif
