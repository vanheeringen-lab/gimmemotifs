/****************************************************************************
 * FILE: evomodel_with_loss.h
 * AUTHOR: JOHN HAWKINS
 * CREATE DATE: 07/02/2008
 * PROJECT: Phylogenetic Motif Scanning Power
 * DESCRIPTION: Evolutionary models for nucleotide substitutions 
 *		allowing motif loss
 * COPYRIGHT: 2008, UQ
 ****************************************************************************/
#ifndef EVOMODEL_WITH_LOSS_H
#define EVOMODEL_WITH_LOSS_H

#include "motif.h"
#include "evomodel.h"

/****************************************************************************
 *  Return the MATRIX of substitution probabilities under the loss model
 *****************************************************************************/ 
MATRIX_T* get_lossmodel_prob_matrix(
	double t, 
	const EVOMODEL_T* motif, 
	const EVOMODEL_T* background, 
	double lambda
);

#endif
