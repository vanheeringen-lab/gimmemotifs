/****************************************************************************
 * FILE: evomodel_with_loss.c
 * AUTHOR: John C. Hawkins
 * CREATE DATE: 05/02/2008
 * PROJECT: PHYLOGENETIC MOTIF SCAN
 * DESCRIPTION: Evolutionary models for nucleotide substitutions that allow for
		motifs to be lost with a specified probability
 * COPYRIGHT: 2008, UQ
 ****************************************************************************/
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "alphabet.h"
#include "evomodel.h"
#include "seq.h"
#include "utils.h"
#include "evomodel_with_loss.h"

// Nucleotide alphabest order as in motif.h
extern char alphabet[];

static double evol_with_loss_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *motif,
  const EVOMODEL_T *background,
  const double lambda
);

MATRIX_T* get_lossmodel_prob_matrix(double t, const EVOMODEL_T* motif, const EVOMODEL_T* background, double lambda) {

  assert(motif != NULL);
  assert(background != NULL);

  MATRIX_T* prob_matrix = NULL;

  // I am going to ignore the halpern bruno option for the moment
  // if (model->use_halpern_bruno == TRUE) {
  //  prob_matrix = matrix_exponential(t, model->rate_matrix);

  int alph_size = get_alph_size(ALPH_SIZE);
  prob_matrix = allocate_matrix(4, 4);
  int i;
  for (i = 0; i < alph_size; i++) {
    char old_base = alphabet[i];
    int j;
    for (j = 0; j < alph_size; j++) {
      char new_base = alphabet[j];
      // Calculate the probability of the old_base mutating
      // to the new_base in time t.
      // we use the new evaluator
      double value = evol_with_loss_evaluator(old_base, new_base, t, motif, background, lambda);
      set_matrix_cell(i, j, value, prob_matrix);
    }
  }
  return prob_matrix;
}

// I HAVE COPIED THIS STUFF TO GET IT TO COMPILE, SOMEHOW IT NEEDS TO BE 
// REFERENCED FROM THE FILE: evomodel.c
// 
/****************************************************************************
  * Data structure representing an evolutionary model for
  * nucleotide substitution
  ****************************************************************************/
struct evomodel {
  MODEL_TYPE_T   type;
  BOOLEAN_T      use_halpern_bruno;
  int            num_params;   // Number of model parameters
  char**         param_names;  // Names of the parameters, useful for debugging
  double*        param_values;
  ARRAY_T*       equil_freqs;
  ARRAY_T*       site_spec_freqs;
  MATRIX_T*      rate_matrix;
  // Function evaluating probability of nucleotide new being substituted
  // for nucleotide old in time t
  double         (*evaluator)(
                   const char old_base, 
                   const char new_base, 
                   const double t,
                   const EVOMODEL_T* model);
    // FIXME: Free string-list multiple times.
};
/****************************************************************************
  * Calculate probability for nucleotide substitution using Tamura-Nei
  * model. See for example Felsenstein, "Inferring Phylogenies",
  * page 203.
  ****************************************************************************/

static double inner_tamura_nei_evaluator(
    const char old_base, 
    const char new_base, 
    const double t,
    const double Pa,
    const double Pc,
    const double Pg,
    const double Pt,
    const double non_specific_rate, 
    const double pur_pur_rate,
    const double pyr_pyr_rate
)
{

  // Check that both bases are in the alphabet.
  assert(strchr(alphabet, old_base));
  assert(strchr(alphabet, new_base));

  // Derive intermediate values used by model.
  const double Pr = Pa + Pg;  // Pyrimidine freq.
  const double Py = Pc + Pt;  // Purine freq
  double Pj = 0.0;            // Equilibrium freq of new base
  switch (new_base) {
    case 'A':
      Pj = Pa;
      break;
    case 'G':
      Pj = Pg;
      break;
    case 'C':
      Pj = Pc;
      break;
    case 'T':
      Pj = Pt;
      break;
  }

  assert(Pr > 0.0);
  assert(Py > 0.0);

  // Evaluate Kronecker delta
  const int Dij = (old_base == new_base) ? 1 : 0;

  // Classify old_base
  double Ai = 0.0;
  if (old_base == 'A' || old_base == 'G' ) {
    // old_base is purine
    Ai = pur_pur_rate;
  } else {
    // old_base is pyrimidine
    Ai = pyr_pyr_rate;
  }

  // Classify new_base
  double SkEjkPk = 0.0; // Value of sum over k of Ejk Pk = (Pr | Py)
  if (new_base == 'A' || new_base == 'G') {
    // new_base is purine
    SkEjkPk = Pr;
  } else {
    // new_base is pyrimidine
    SkEjkPk = Py;
  }

  int Eij = 0;      // Value of Watson-Kronecker given old_base, new_base
  // Are both bases R or Y?
  if ((old_base == 'A' || old_base == 'G') && 
      (new_base == 'A' || new_base =='G')) {
    // Both purines
    Eij = 1;
  } else if ((old_base == 'C' || old_base == 'T') && 
             (new_base == 'C' || new_base == 'T')) {
    // Both pyrimidines
    Eij = 1;
  } else {
    Eij = 0;
  }

  const double first_term =
    (Dij == 0) ? 0 : exp(-((Ai + non_specific_rate) * t));
  const double second_term =
    (Eij == 0) ? 0 : (exp(-(non_specific_rate * t)) 
  		      * (1.0 - exp(-(Ai * t))) * (Pj)/SkEjkPk);
  const double third_term = (1.0 - exp(-(non_specific_rate * t))) * Pj;
  const double value = first_term + second_term + third_term;

  return value;
}

// END COPY



/****************************************************************************
  * Evaluate the inner sum of integral of the loss model
  *
  * The integral of the loss version of the Tamura-Nei model
  * consists of a a sum over all possible bases {a} that occur at
  * the point of loss. The function inside this sum is the
  * behemouth printed below. Each of the parameters must be
  * calaculated for {a} and then fed into this function.
  * See the supplementary material for the paper for an explanation
  * of each of the parameters.
  ****************************************************************************/
static double evaluate_integral_inner_sum(
    const double t,
    const double b,
    const double c,
    const double D,
    const double L,
    //const double m,
    //const double n,
    const double p, 
    const double q,
    const double r,
    const double s,
    const double u,
    const double v,
    const double w,
    const double y
)
{
/*
 Lbc/(-L-s-u+r+v) e^{-Lt -st - ut - rD + rt - vD + vt} \\
+ Lcq/(-L-s-u+v) e^{-Lt - st - ut - vD + vt} \\
- Lcq/(-L-s-u+v+r) e^{-Lt -st -ut -vD +vt -rD +rt} \\
+ Lcy/(-L-s-u) e^{-Lt -st -ut} \\
- Lcy/(-L-s-u+v) e^{-Lt -st - ut -vD + vt}\\
+ Lbp/(-L-u+r+v) e^{-Lt -ut - rD + rt  - vD + vt } \\
- Lbp/(-L-u+r+v-s) e^{-Lt -ut - rD + rt  - vD + vt -st} \\
+ Lpq/(-L-u+v) e^{-Lt - ut -vD +vt} \\
- Lpq/(-L-u+v+r) e^{-Lt - ut -vD +vt -rD +rt}\\
- Lpq/(-L-u+v-s) e^{-Lt - ut -vD +vt -st}\\
+ Lpq/(-L-u+v-s+r) e^{-Lt - ut -vD +vt -st-rD +rt}	\\
+ Lpy/(-L-u) e^{-Lt -ut} \\
- Lpy/(-L-u+v) e^{-Lt -ut -vD +vt} \\
- Lpy/(-L-u-s) e^{-Lt -ut-st} \\
+ Lpy/(-L-u-s+v) e^{-Lt -ut -st -vD +vt}\\
+ Lbw/(-L+r+v) e^{-Lt - rD + rt - vD  + vt} \\
- Lbw/(-L+r+v-u) e^{ -Lt - rD + rt - vD  + vt -ut}\\
+ Lqw/(-L+v) e^{-Lt-vD +vt} \\
- Lqw/(-L+v+r) e^{-Lt-vD +vt -rD + rt} \\
- Lqw/(-L+v-u) e^{-Lt-vD +vt -ut} \\
+ Lqw/(-L+v-u+r) e^{-Lt-vD +vt -ut- rD + rt}\\
+ Lwy/(-L) e^{-Lt} \\
- Lwy/(-L-u) e^{-Lt -ut} \\
- Lwy/(-L+v) e^{-Lt -vD +vt} \\
+ Lwy/(-L-u+v) e^{-Lt -ut -vD +vt}
*/
	return ((L*b*c)/(-L-s-u+r+v)) *  exp( -(L*t) -(s*t) - (u*t) - (r*D) + (r*t) - (v*D) + (v*t) )
		+ ((L*c*q)/(-L-s-u+v)) * exp( -(L*t) - (s*t) - (u*t) - (v*D) + (v*t) )
		- ((L*c*q)/(-L-s-u+v+r)) * exp(-(L*t) - (s*t) - (u*t) - (v*D) + (v*t) -(r*D) +(r*t) )
		+ ((L*c*y)/(-L-s-u)) * exp(-(L*t) - (s*t) - (u*t) )
		- ((L*c*y)/(-L-s-u+v)) * exp( -(L*t) - (s*t) - (u*t) -(v*D) + (v*t) )
		+ ((L*b*p)/(-L-u+r+v)) * exp( -(L*t) - (u*t) - (r*D) + (r*t)  - (v*D) + (v*t) ) 
		- ((L*b*p)/(-L-u+r+v-s)) * exp( -(L*t) - (u*t) - (r*D) + (r*t)  - (v*D) + (v*t) - (s*t) )
		+ ((L*p*q)/(-L-u+v)) * exp( -(L*t) - (u*t) -(v*D ) +(v*t) ) 
		- ((L*p*q)/(-L-u+v+r)) *  exp( -(L*t) - (u*t) -(v*D ) +(v*t) -(r*D) +(r*t) )
		- ((L*p*q)/(-L-u+v-s)) * exp( -(L*t) - (u*t) - (v*D) + (v*t) -(s*t) )
		+ ((L*p*q)/(-L-u+v-s+r)) * exp(-(L*t) - (u*t) - (v*D) + (v*t) -(s*t) - (r*D) + (r*t) )
		+ ((L*p*y)/(-L-u)) * exp( -(L*t) - (u*t) ) 
		- ((L*p*y)/(-L-u+v)) * exp( -(L*t) - (u*t) - (v*D) + (v*t) )
		- ((L*p*y)/(-L-u-s)) * exp( -(L*t) - (u*t) -(s*t) )
		+ ((L*p*y)/(-L-u-s+v)) * exp( -(L*t) - (u*t) -(s*t) -(v*D) + (v*t) )
		+ ((L*b*w)/(-L+r+v)) * exp( -(L*t) - (r*D) + (r*t) - (v*D)  + (v*t) )
		- ((L*b*w)/(-L+r+v-u)) * exp( -(L*t) - (r*D) + (r*t) - (v*D)  + (v*t) -(u*t) )
		+ ((L*q*w)/(-L+v)) * exp( -(L*t) -(v*D) +(v*t) )
		- ((L*q*w)/(-L+v+r)) * exp( -(L*t) -(v*D) +(v*t) -(r*D) + (r*t) )
		- ((L*q*w)/(-L+v-u)) * exp( -(L*t) -(v*D) +(v*t) -(u*t) )
		+ ((L*q*w)/(-L+v-u+r)) * exp( -(L*t) -(v*D) +(v*t) -(u*t) -(r*D) + (r*t) )
		+ ((L*w*y)/(-L)) * exp(-(L*t) )
		- ((L*w*y)/(-L-u)) * exp( -(L*t) -(u*t) )
		- ((L*w*y)/(-L+v)) * exp( -(L*t) -(v*D) +(v*t) )
		+ ((L*w*y)/(-L-u+v)) * exp( -(L*t) -(u*t) -(v*D) +(v*t) );
}

/****************************************************************************
  * Calculate probability for nucleotide substitution using a model
  * that allows the specified motif model to be lost with a certain
  * probability over the specified evolutionary time.
  * The calculation is based on an integration of the Tamura-Nei
  * model (Felsenstein, "Inferring Phylogenies", page 203).
  * 
  ****************************************************************************/
static double inner_tamura_nei_evaluator_with_loss(
    const char old_base, 
    const char new_base, 
    const double t,
    const double MPa,
    const double MPc,
    const double MPg,
    const double MPt,
    const double Mnon_specific_rate, 
    const double Mpur_pur_rate,
    const double Mpyr_pyr_rate,
    const double BPa,
    const double BPc,
    const double BPg,
    const double BPt,
    const double Bnon_specific_rate, 
    const double Bpur_pur_rate,
    const double Bpyr_pyr_rate,
    const double lambda
)
{

  // Check that both bases are in the alphabet.
  assert(strchr(alphabet, old_base));
  assert(strchr(alphabet, new_base));

  // We need to derive all the intermediate values used by model.
  // First for the motif
  const double MPr = MPa + MPg;  // Pyrimidine freq.
  const double MPy = MPc + MPt;  // Purine freq
  double MPj = 0.0;              // Equilibrium freq of new base
  switch (new_base) {
    case 'A':
      MPj = MPa;
      break;
    case 'G':
      MPj = MPg;
      break;
    case 'C':
      MPj = MPc;
      break;
    case 'T':
      MPj = MPt;
      break;
  }
  assert(MPr > 0.0);
  assert(MPy > 0.0);

  // Now for the background model
  const double BPr = BPa + BPg;  // Pyrimidine freq.
  const double BPy = BPc + BPt;  // Purine freq
  double BPj = 0.0;            // Equilibrium freq of new base
  switch (new_base) {
    case 'A':
      BPj = BPa;
      break;
    case 'G':
      BPj = BPg;
      break;
    case 'C':
      BPj = BPc;
      break;
    case 'T':
      BPj = BPt;
      break;
  }
  assert(BPr > 0.0);
  assert(BPy > 0.0);

  double noLossValue = inner_tamura_nei_evaluator( old_base, new_base, t, MPa, MPc, MPg, MPt, Mnon_specific_rate, Mpur_pur_rate, Mpyr_pyr_rate);

  // The first part of the solution is the value of the probability of no loss occuring
  double first_term = noLossValue * exp( -(lambda * t) );

  double second_term = 0.0;

  double totalProbOfLoss = 1;//1 - exp( -(lambda * t) );
  if(totalProbOfLoss > 0) {
  // The second part is the integral over all possible ways the site could be lost.
  // This is a sum over all possible residues that could occur at the loss point.
  int z=0;
  for(z=0; z<4; z++) {
	char middle_base = alphabet[z];
  	double IPj = 0.0;            // Equilibrium freq of intermediate base under motif
  	switch (middle_base) {
    	case 'A':
      		IPj = MPa;
      		break;
    	case 'G':
      		IPj = MPg;
      		break;
    	case 'C':
      		IPj = MPc;
      		break;
   	case 'T':
      		IPj = MPt;
      		break;
  	}

	// Now we generate all the intermediate values needed by the inner sum 
	// of the solution to the integral

        // First the delta function which indicates whether the initial base is equal to the intermediate
	// c = \delta_{\sigma_1,a}
	const int c = (old_base == middle_base) ? 1 : 0;
	// Now the delta function for whether the middle base is equal to the final
	// b = \delta_{a,\sigma_i}
	const int b = (middle_base == new_base) ? 1 : 0;

	// The next value provides the proportion of purine/pyrimidine frequency
 	// contributed by the middle base if both the initial and middle are of the 
	// same type, zero otherwise.
	// p = \Big( \frac{\pi^M_a\epsilon_{{\sigma_1},a}}{\sum_k \epsilon_{a,k} \pi^M_k} \Big
	double SkEjkPk = 0.0; // Value of sum over k of Ejk Pk = (Pr | Py)
  	if (middle_base == 'A' || middle_base == 'G') {
    		// middle is purine
    		SkEjkPk = MPr;
  	} else {
    		// middle is pyrimidine
    		SkEjkPk = MPy;
  	}
	int Eij = 0;      // Value of Watson-Kronecker given old_base, middle_base
  	// Are both bases R or Y?
  	if ((old_base == 'A' || old_base == 'G') && (middle_base == 'A' || middle_base =='G')) {
    		// Both purines
    		Eij = 1;
  	} else if ((old_base == 'C' || old_base == 'T') &&  (middle_base == 'C' || middle_base == 'T')) {
    		// Both pyrimidines
    		Eij = 1;
  	} else {
    		Eij = 0;
  	}
	const double p = Eij* IPj/SkEjkPk;

	// The next value provides the proportion of purine/pyrimidine frequency
	// under the background model
 	// contributed by the final base if both the middle and final are of the 
	// same type, zero otherwise.
	// q = \Big( \frac{\pi^B_{\sigma_i}\epsilon_{a,\sigma_i}}{\sum_k \epsilon_{\sigma_i,k} \pi^B_k} \Big) 
	SkEjkPk = 0.0; // Value of sum over k of Ejk Pk = (Pr | Py)
  	if (new_base == 'A' || new_base == 'G') {
    		// new_base is purine
    		SkEjkPk = BPr;
  	} else {
    		// new_base is pyrimidine
    		SkEjkPk = BPy;
  	}
	// Value of Watson-Kronecker given middle_base, new_base
	// Are both bases R or Y?
  	if ((middle_base == 'A' || middle_base == 'G') && (new_base == 'A' || new_base =='G')) {
    		// Both purines
    		Eij = 1;
  	} else if ((middle_base == 'C' || middle_base == 'T') && (new_base == 'C' || new_base == 'T')) {
    		// Both pyrimidines
    		Eij = 1;
  	} else {
    		Eij = 0;
  	}
	const double q = Eij* BPj/SkEjkPk;

	// The next value is the pur-pur or pyr-pyr rate of the intermediate base
	// under the background model
	// r = \alpha^B_a
	double r = 0;
	if (middle_base == 'A' || middle_base == 'G' ) {
    		// middle_base is purine
    		r = Bpur_pur_rate;
  	} else {
    		// middle_base is pyrimidine
    		r = Bpyr_pyr_rate;
  	}
	// The next value is the pur-pur or pyr-pyr rate of the initial base
	// under the motif model
	// s = \alpha^M_{\sigma_1}
	double s = 0;
	if (old_base == 'A' || old_base == 'G' ) {
    		// old_base is purine
    		s = Mpur_pur_rate;
  	} else {
    		// old_base is pyrimidine
    		s = Mpyr_pyr_rate;
  	}

	// The next value is the non specific rate of the motif model 
	// u = \beta^M
	const double u = Mnon_specific_rate;
	// and the non specific rate of the background model 
	// v = \beta^B
	const double v = Bnon_specific_rate;

	// The next value is the equilibrium frequency of the middle amino acid
	// under the motif model
	// w = \pi^M_a
	const double w = IPj;

	// The next value is the equilibrium frequency of the final amino acid
	// under the background model
	// y = \pi^B_{\sigma_i}
	const double y = BPj;

	const int D = t;
	const double L = lambda;

	// Now that all the intermediates are calculated
	// we plug them in to evaluate the integral at the two values for the 
	// region being integrated over
	double upper = evaluate_integral_inner_sum( t, b, c, D, L, p, q, r, s, u, v, w, y );
	double lower = evaluate_integral_inner_sum( 0, b, c, D, L, p, q, r, s, u, v, w, y );
	//printf("%f - %f\n", upper, lower);
	second_term += (upper - lower);
  }
     second_term=second_term*totalProbOfLoss;
  }	
  const double value = first_term + second_term;

  return value;
}




/****************************************************************************
  * FUNCTION FOR PERFORMING THE CALCULATIONS USING A NUMERICAL APPROXIMATION 
  ****************************************************************************/


/****************************************************************************
  * Evaluate the inner part of the function being integrated numerically
  *
  * The function of the loss version of the Tamura-Nei model
  * consists of a a sum over all possible bases {a} that occur at
  * the point of loss. The function inside this sum is the
  * one below, it is integrated in the other version. 
  * Each of the parameters must be
  * calaculated for {a} and then fed into this function.
  * See the supplementary material for the paper for an explanation
  * of each of the parameters.
  ****************************************************************************/
static double evaluate_function_inner_sum(
    const double t,
    const double b,
    const double c,
    const double D,
    const double L,
    const double p, 
    const double q,
    const double r,
    const double s,
    const double u,
    const double v,
    const double w,
    const double y
)
{
		return ( (exp(-(s + u)*t) * c + exp(-u*t) * (1 - exp(-s*t)) * p + (1 - exp(-u*t)) * w) * 
			 (exp(-(r + v)*(D-t)) * b + exp(-v*(D-t)) * (1 - exp(-r*(D-t))) * q + (1 - exp(-v*(D-t))) * y ) );
}

/****************************************************************************
  * Calculate probability for nucleotide substitution using a numerical
  * approximation to the required integral
  ****************************************************************************/
static double inner_tamura_nei_evaluator_with_loss_numerical(
    const char old_base, 
    const char new_base, 
    const double t,
    const double MPa,
    const double MPc,
    const double MPg,
    const double MPt,
    const double Mnon_specific_rate, 
    const double Mpur_pur_rate,
    const double Mpyr_pyr_rate,
    const double BPa,
    const double BPc,
    const double BPg,
    const double BPt,
    const double Bnon_specific_rate, 
    const double Bpur_pur_rate,
    const double Bpyr_pyr_rate,
    const double lambda
)
{

  // Check that both bases are in the alphabet.
  assert(strchr(alphabet, old_base));
  assert(strchr(alphabet, new_base));
  // We need to derive all the intermediate values used by model.
  // First for the motif
  const double MPr = MPa + MPg;  // Pyrimidine freq.
  const double MPy = MPc + MPt;  // Purine freq
  double MPj = 0.0;              // Equilibrium freq of new base
  switch (new_base) {
    case 'A':
      MPj = MPa;
      break;
    case 'G':
      MPj = MPg;
      break;
    case 'C':
      MPj = MPc;
      break;
    case 'T':
      MPj = MPt;
      break;
  }
  assert(MPr > 0.0);
  assert(MPy > 0.0);

  // Now for the background model
  const double BPr = BPa + BPg;  // Pyrimidine freq.
  const double BPy = BPc + BPt;  // Purine freq
  double BPj = 0.0;            // Equilibrium freq of new base
  switch (new_base) {
    case 'A':
      BPj = BPa;
      break;
    case 'G':
      BPj = BPg;
      break;
    case 'C':
      BPj = BPc;
      break;
    case 'T':
      BPj = BPt;
      break;
  }
  assert(BPr > 0.0);
  assert(BPy > 0.0);

  double noLossValue = inner_tamura_nei_evaluator( old_base, new_base, t, MPa, MPc, MPg, MPt, Mnon_specific_rate, Mpur_pur_rate, Mpyr_pyr_rate);

  // The first part of the solution is the value of the probability of no loss occuring
  double first_term = noLossValue * exp( -(lambda * t) );
  double second_term = 0.0;
  double totalProbOfLoss = 1;//1 - exp( -(lambda * t) );

  if(totalProbOfLoss > 0) {

  // The second part is the integral over all possible ways the site could be lost.
  // This is a sum over all possible residues that could occur at the loss point.

  int numSlices = 1000;
  double sliceWidth = t/numSlices;
  int i=0;
  for(i=0; i<numSlices; i++) {
	double evalAt = (sliceWidth*i) + (sliceWidth/2);
        int z=0;
  	for(z=0; z<4; z++) {
		char middle_base = alphabet[z];
  		double IPj = 0.0;            // Equilibrium freq of intermediate base under motif
  		switch (middle_base) {
    		case 'A':
      			IPj = MPa;
      			break;
    		case 'G':
      			IPj = MPg;
      			break;
    		case 'C':
      			IPj = MPc;
      			break;
   		case 'T':
      			IPj = MPt;
      			break;
  		}

		// Now we generate all the intermediate values needed by the inner sum 
		// of the solution to the integral

        	// First the delta function which indicates whether the initial base is equal to the intermediate
		// c = \delta_{\sigma_1,a}
		const int c = (old_base == middle_base) ? 1 : 0;
		// Now the delta function for whether the middle base is equal to the final
		// b = \delta_{a,\sigma_i}
		const int b = (middle_base == new_base) ? 1 : 0;

		// The next value provides the proportion of purine/pyrimidine frequency
 		// contributed by the middle base if both the initial and middle are of the 
		// same type, zero otherwise.
		// p = \Big( \frac{\pi^M_a\epsilon_{{\sigma_1},a}}{\sum_k \epsilon_{a,k} \pi^M_k} \Big
		double SkEjkPk = 0.0; // Value of sum over k of Ejk Pk = (Pr | Py)
  		if (middle_base == 'A' || middle_base == 'G') {
    			// middle is purine
    			SkEjkPk = MPr;
  		} else {
    			// middle is pyrimidine
    			SkEjkPk = MPy;
  		}
		int Eij = 0;      // Value of Watson-Kronecker given old_base, middle_base
  		// Are both bases R or Y?
  		if ((old_base == 'A' || old_base == 'G') && (middle_base == 'A' || middle_base =='G')) {
    			// Both purines
    			Eij = 1;
  		} else if ((old_base == 'C' || old_base == 'T') &&  (middle_base == 'C' || middle_base == 'T')) {
    			// Both pyrimidines
    			Eij = 1;
  		} else {
    			Eij = 0;
  		}
		const double p = Eij* IPj/SkEjkPk;

		// The next value provides the proportion of purine/pyrimidine frequency
		// under the background model
 		// contributed by the final base if both the middle and final are of the 
		// same type, zero otherwise.
		// q = \Big( \frac{\pi^B_{\sigma_i}\epsilon_{a,\sigma_i}}{\sum_k \epsilon_{\sigma_i,k} \pi^B_k} \Big) 
		SkEjkPk = 0.0; // Value of sum over k of Ejk Pk = (Pr | Py)
  		if (new_base == 'A' || new_base == 'G') {
    			// new_base is purine
    			SkEjkPk = BPr;
  		} else {
    			// new_base is pyrimidine
    			SkEjkPk = BPy;
  		}
		// Value of Watson-Kronecker given middle_base, new_base
		// Are both bases R or Y?
  		if ((middle_base == 'A' || middle_base == 'G') && (new_base == 'A' || new_base =='G')) {
    			// Both purines
    			Eij = 1;
  		} else if ((middle_base == 'C' || middle_base == 'T') && (new_base == 'C' || new_base == 'T')) {
    			// Both pyrimidines
    			Eij = 1;
  		} else {
    			Eij = 0;
  		}
		const double q = Eij* BPj/SkEjkPk;

		// The next value is the pur-pur or pyr-pyr rate of the intermediate base
		// under the background model
		// r = \alpha^B_a
		double r = 0;
		if (middle_base == 'A' || middle_base == 'G' ) {
    			// middle_base is purine
    			r = Bpur_pur_rate;
  		} else {
    			// middle_base is pyrimidine
    			r = Bpyr_pyr_rate;
  		}
		// The next value is the pur-pur or pyr-pyr rate of the initial base
		// under the motif model
		// s = \alpha^M_{\sigma_1}
		double s = 0;
		if (old_base == 'A' || old_base == 'G' ) {
    			// old_base is purine
    			s = Mpur_pur_rate;
  		} else {
    			// old_base is pyrimidine
    			s = Mpyr_pyr_rate;
  		}

		// The next value is the non specific rate of the motif model 
		// u = \beta^M
		const double u = Mnon_specific_rate;
		// and the non specific rate of the background model 
		// v = \beta^B
		const double v = Bnon_specific_rate;

		// The next value is the equilibrium frequency of the middle amino acid
		// under the motif model
		// w = \pi^M_a
		const double w = IPj;

		// The next value is the equilibrium frequency of the final amino acid
		// under the background model
		// y = \pi^B_{\sigma_i}
		const double y = BPj;

		const int D = t;
		const double L = lambda;

		// Now that all the intermediates are calculated
		// we plug them in to evaluate the integral at the two values for the 
		// region being integrated over
		second_term += evaluate_function_inner_sum( evalAt, b, c, D, L, p, q, r, s, u, v, w, y ) * sliceWidth;
  	}
     }
  }
  const double value = first_term + (totalProbOfLoss*second_term);

  return value;
}





/****************************************************************************
  * Calculate probability for nucleotide substitution using a numerical
  * approximation to the required integral
  ****************************************************************************/
static double inner_tamura_nei_evaluator_with_loss_numerical2(
    const char old_base, 
    const char new_base, 
    const double t,
    const double MPa,
    const double MPc,
    const double MPg,
    const double MPt,
    const double Mnon_specific_rate, 
    const double Mpur_pur_rate,
    const double Mpyr_pyr_rate,
    const double BPa,
    const double BPc,
    const double BPg,
    const double BPt,
    const double Bnon_specific_rate, 
    const double Bpur_pur_rate,
    const double Bpyr_pyr_rate,
    const double lambda
)
{

  // Check that both bases are in the alphabet.
  assert(strchr(alphabet, old_base));
  assert(strchr(alphabet, new_base));

  double noLossValue = inner_tamura_nei_evaluator( old_base, new_base, t, MPa, MPc, MPg, MPt, Mnon_specific_rate, Mpur_pur_rate, Mpyr_pyr_rate);

  // The first part of the solution is the value of the probability of no loss occuring
  double first_term = noLossValue * exp( -(lambda * t) );
  double second_term = 0.0;
  double totalProbOfLoss = 1;//1 - exp( -(lambda * t) );

  // The second part is the integral over all possible ways the site could be lost.
  // This is a sum over all possible residues that could occur at the loss point.
  if(totalProbOfLoss > 0) {
	// Now we perform the numerical approximation
	int numSlices = 100000;
	double sliceWidth = t/numSlices;
	second_term = 0.0;
        int i=0;
  	for(i=0; i<numSlices; i++) {
	   double evalAt = (sliceWidth*i) + (sliceWidth/2);
  	   double expPart = ( lambda * exp( -(lambda * evalAt) ));

	   double innerSum = 0.0;
	   int z=0;
  	   for(z=0; z<4; z++) {
	      char middle_base = alphabet[z];
	      double motifPart = inner_tamura_nei_evaluator( old_base, middle_base, evalAt, MPa, MPc, MPg, MPt, Mnon_specific_rate, Mpur_pur_rate, Mpyr_pyr_rate);
	      double basePart = inner_tamura_nei_evaluator( middle_base, new_base, (t-evalAt), BPa, BPc, BPg, BPt, Bnon_specific_rate, Bpur_pur_rate, Bpyr_pyr_rate);
	      innerSum += motifPart * basePart;

	   }
           double functValue = expPart * innerSum;
	   second_term += functValue * sliceWidth;
  	}
     
  }
  const double value = first_term + (totalProbOfLoss*second_term);

  return value;
}




/****************************************************************************
  * VISIBLE FUNCTION FOR OBTAINING THE SUB PROBS WITH LOSS
  ****************************************************************************/


/****************************************************************************
  * Calculate probability for nucleotide substitution using 
  * the general felsenstein model allowing for motif site loss.
  ****************************************************************************/
static double evol_with_loss_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *motif,
  const EVOMODEL_T *background, 
  const double lambda
)
{
  assert(motif != NULL);
  assert(background != NULL);
  assert(motif->param_names != NULL);
  assert(motif->param_values != NULL);
  assert(background->param_names != NULL);
  assert(background->param_values != NULL);


  // OK now we need to work out the parameters to feed into the 
  // the evaluator function: inner_tamura_nei_evaluator_with_loss
  // These will depend on the type of each of the models

  // Start with the motif model: We assume Jukes-Cantor
  // and then change the values depending on the model type
  double rate = (motif->param_values)[0];
  double MPa = 0.25;
  double MPc = 0.25;
  double MPg = 0.25;
  double MPt = 0.25;
  double Mnon_specific_rate = 4.0 * rate / 3.0;
  double Mpur_pur_rate = 0.0;
  double Mpyr_pyr_rate = 0.0;
  // Some variable for continuous reuse 
  double Mtransition_transversion = 0.0;
  double MPr = 0;
  double MPy = 0;
  //fprintf( stderr, "FG Model is type %d from options K2 [%d] F81 [%d] F84 [%d] HKY [%d] \n", motif->type, K2_MODEL, F81_MODEL, F84_MODEL, HKY_MODEL);
  if(K2_MODEL == motif->type) {
	Mtransition_transversion = (motif->param_values)[1];
  	MPr = MPa + MPg;
  	MPy = MPc + MPt;

  	Mpyr_pyr_rate = rate * ((MPr * MPy * Mtransition_transversion) 
    		- (MPa * MPg) - (MPc * MPt)) / (2.0 * (1.0 + Mtransition_transversion)
	  	* ((MPy * MPa * MPg) + (MPr * MPc * MPt)));
  	Mpur_pur_rate = Mpyr_pyr_rate;
  	Mnon_specific_rate = rate / ( 2.0 * MPr * MPy * (1.0 + Mtransition_transversion));

  } else if (F81_MODEL == motif->type) {
  	MPa = motif->param_values[1];
  	MPc = motif->param_values[2];
  	MPg = motif->param_values[3];
  	MPt = motif->param_values[4];
  	MPr = MPa + MPg;
  	MPy = MPc + MPt;

  	Mpyr_pyr_rate = 0;
  	Mpur_pur_rate = 0;
  	Mnon_specific_rate = rate / (3.0 * MPr * MPy);

  } else if (F84_MODEL == motif->type) {
  	MPa = motif->param_values[1];
  	MPc = motif->param_values[2];
  	MPg = motif->param_values[3];
  	MPt = motif->param_values[4];
  	MPr = MPa + MPg;
  	MPy = MPc + MPt;
  	Mtransition_transversion = (motif->param_values)[5];

  	assert(MPy > 0.0);
  	assert(MPr > 0.0);

  	Mpyr_pyr_rate = rate * ((MPr * MPy * Mtransition_transversion) 
			- (MPa * MPg) - (MPc * MPt)) / (2.0 * (1.0 + Mtransition_transversion)
	  		* ((MPy * MPa * MPg) + (MPr * MPc * MPt)));
  	Mpur_pur_rate = Mpyr_pyr_rate;
  	Mnon_specific_rate =  rate / ( 2.0 * MPr * MPy * (1.0 + Mtransition_transversion));

  } else if (HKY_MODEL == motif->type) {
	MPa = motif->param_values[1];
  	MPc = motif->param_values[2];
  	MPg = motif->param_values[3];
  	MPt = motif->param_values[4];
  	MPr = MPa + MPg;
  	MPy = MPc + MPt;
	Mtransition_transversion = (motif->param_values)[5];
	assert(MPy > 0.0);
  	assert(MPr > 0.0);
	Mpyr_pyr_rate = rate * ((MPr * MPy * Mtransition_transversion) 
    			- (MPa * MPg) - (MPc * MPt)) / (2.0 * (1.0 + Mtransition_transversion)
	  		* ((MPy * MPa * MPg * MPr / MPy) + (MPr * MPc * MPt)));
  	Mpur_pur_rate = Mpyr_pyr_rate * MPr / MPy;
	Mnon_specific_rate = rate / ( 2.0 * MPr * MPy * (1.0 + Mtransition_transversion));
  }


  // Now with the background model: We assume Jukes-Cantor
  // and then change the values depending on the model type
  rate = (background->param_values)[0];
  double BPa = 0.25;
  double BPc = 0.25;
  double BPg = 0.25;
  double BPt = 0.25;
  double Bnon_specific_rate = 4.0 * rate / 3.0;
  double Bpur_pur_rate = 0.0;
  double Bpyr_pyr_rate = 0.0;
  // Some variable for continuous reuse 
  double Btransition_transversion = 0.0;
  double BPr = 0;
  double BPy = 0;
  if(K2_MODEL == background->type) {
	Btransition_transversion = (background->param_values)[1];
  	BPr = BPa + BPg;
  	BPy = BPc + BPt;

  	Bpyr_pyr_rate = rate * ((BPr * BPy * Btransition_transversion) 
    		- (BPa * BPg) - (BPc * BPt)) / (2.0 * (1.0 + Btransition_transversion)
	  	* ((BPy * BPa * BPg) + (BPr * BPc * BPt)));
  	Bpur_pur_rate = Bpyr_pyr_rate;
  	Bnon_specific_rate = rate / ( 2.0 * BPr * BPy * (1.0 + Btransition_transversion));

  } else if (F81_MODEL == background->type) {
  	BPa = background->param_values[1];
  	BPc = background->param_values[2];
  	BPg = background->param_values[3];
  	BPt = background->param_values[4];
  	BPr = BPa + BPg;
  	BPy = BPc + BPt;

  	Bpyr_pyr_rate = 0;
  	Bpur_pur_rate = 0;
  	Bnon_specific_rate = rate / (3.0 * BPr * BPy);

  } else if (F84_MODEL == background->type) {
  	BPa = background->param_values[1];
  	BPc = background->param_values[2];
  	BPg = background->param_values[3];
  	BPt = background->param_values[4];
  	BPr = BPa + BPg;
  	BPy = BPc + BPt;
  	Btransition_transversion = (background->param_values)[5];

  	assert(BPy > 0.0);
  	assert(BPr > 0.0);

  	Bpyr_pyr_rate = rate * ((BPr * BPy * Btransition_transversion) 
			- (BPa * BPg) - (BPc * BPt)) / (2.0 * (1.0 + Btransition_transversion)
	  		* ((BPy * BPa * BPg) + (BPr * BPc * BPt)));
  	Bpur_pur_rate = Bpyr_pyr_rate;
  	Bnon_specific_rate =  rate / ( 2.0 * BPr * BPy * (1.0 + Btransition_transversion));

  } else if (HKY_MODEL == motif->type) {
  	BPa = background->param_values[1];
  	BPc = background->param_values[2];
  	BPg = background->param_values[3];
  	BPt = background->param_values[4];
  	BPr = BPa + BPg;
  	BPy = BPc + BPt;
  	Btransition_transversion = (background->param_values)[5];

  	assert(BPy > 0.0);
  	assert(BPr > 0.0);

  	Bpyr_pyr_rate = rate * ((BPr * BPy * Btransition_transversion) 
			- (BPa * BPg) - (BPc * BPt)) / (2.0 * (1.0 + Btransition_transversion)
	  		* ((BPy * BPa * BPg * BPr / BPy) + (BPr * BPc * BPt)));
  	Bpur_pur_rate = Bpyr_pyr_rate * BPr / BPy;
	Bnon_specific_rate = rate / ( 2.0 * BPr * BPy * (1.0 + Btransition_transversion));
  }


  return inner_tamura_nei_evaluator_with_loss_numerical2(
    old_base, 
    new_base, 
    t,
    MPa,
    MPc,
    MPg,
    MPt,
    Mnon_specific_rate, 
    Mpur_pur_rate,
    Mpyr_pyr_rate,
    BPa,
    BPc,
    BPg,
    BPt,
    Bnon_specific_rate, 
    Bpur_pur_rate,
    Bpyr_pyr_rate,
    lambda
  );

}
