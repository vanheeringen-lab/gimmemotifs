/****************************************************************************
 * FILE: evomodel.h
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 11/18/2004
 * PROJECT: EVOMCAST
 * DESCRIPTION: Evolutionary models for nucleotide substitutions
 * COPYRIGHT: 2004, UW
 ****************************************************************************/
#ifndef EVOMODEL_H
#define EVOMODEL_H

#include "motif.h"

typedef enum {
  SKIP_GAPS, 
  FIXED_GAP_COST, 
  WILDCARD_GAP, 
  MIN_GAPS 
} GAP_SUPPORT_T;

typedef enum {
  SINGLE_MODEL, 
  AVERAGE_MODEL, 
  JC_MODEL, 
  K2_MODEL, 
  F81_MODEL,
  F84_MODEL, 
  HKY_MODEL, 
  TAMURA_NEI_MODEL 
} MODEL_TYPE_T;

typedef enum {
  SIMULATED_COL_FREQS,
  EMPIRICAL_COL_FREQS,
  PRECOMPUTED_COL_FREQS
} COLUMN_FREQS_TYPE_T;

// An evomodel object.
typedef struct evomodel EVOMODEL_T;

// Get methods for an EVOMODEL_T
const char* get_model_name(const EVOMODEL_T* model);
int get_num_model_params(const EVOMODEL_T* model);
char** get_model_param_names(const EVOMODEL_T* model);
double* get_model_param_values(const EVOMODEL_T* model);
double get_model_equil_freq(char base, const EVOMODEL_T* model);
ARRAY_T* get_model_equil_freqs(const EVOMODEL_T* model);
MODEL_TYPE_T get_model_type(const EVOMODEL_T* model);
MATRIX_T* get_model_prob_matrix(const double t, const EVOMODEL_T* model);
MATRIX_T* old_get_model_prob_matrix(const EVOMODEL_T* model, double t);
MATRIX_T* get_model_rate_matrix(EVOMODEL_T* model);
MATRIX_T* matrix_exponential(const double t, MATRIX_T *a);
BOOLEAN_T uses_halpern_bruno(const EVOMODEL_T* model);

/*******************************************************************
 * Build a set of evolutionary models for the positions of a motif. 
 ********************************************************************/
EVOMODEL_T** make_motif_models(
  MOTIF_T* motif, 
  ARRAY_T* bg_freqs,
  MODEL_TYPE_T model_type,
  double fg_rate,
  double bg_rate,
  double purine_pyrimidine,
  double transition_transversion,
  BOOLEAN_T use_halpern_bruno
); 

/*******************************************************************
 * Build an evolutionary model. 
 * The type of model to be built is specified by the
 * model_type parameter.  The supported model types are:
 *    psfm = Position specific freq. matrix
 *    jc = Jukes-Cantor
 *    k2 = Kimura 2 parameter
 *    f81 = Felsenstein 81
 *    f84 = Felsenstein 84
 *    hky = HKY
 *    tn = Tamura-Nei
 * Each of these models requires a different set of parameters.
 * In order to support the most general model all of the parameters
 * will be required even if the model chosen doesn't make use of them.
 *
 * If an array of site-specific frequencies is provided the Halpern-Brun
 * modification of the model will be used.
 ********************************************************************/
EVOMODEL_T* make_model(
  const MODEL_TYPE_T model_type, 
  const double rate,
  const double transition_transversion,
  const double purine_pyrimidine,
  ARRAY_T* equil_freqs,
  ARRAY_T* position_specific_freqs
);

/****************************************************************************
 *  Allocate a data structure for a PSFM model
 *    equil_freqs contains the equilibrium frequencies for the model.
 ****************************************************************************/ 
EVOMODEL_T* make_psfm_model(ARRAY_T* equil_freqs);

/****************************************************************************
 *  Allocate a data structure for average PSFM model
 *    equil_freqs contains the equilibrium frequencies for the model.
 ****************************************************************************/ 
EVOMODEL_T* make_avg_psfm_model(ARRAY_T* equil_freqs);

/****************************************************************************
 *  Allocate a data structure for a Jukes-Cantor model
 ****************************************************************************/
EVOMODEL_T* make_jc_model(
  const double rate, 
  ARRAY_T* equil_freqs
);

/****************************************************************************
 *  Allocate a data structure a Kimura 2-parameter model
 ****************************************************************************/
EVOMODEL_T* make_k2_model(
  const double rate, 
  ARRAY_T* equil_freqs, 
  const double transition_transversion
);

/****************************************************************************
 *  Allocate a data structure for an F81 model
 *    rate is overall mutation rate
 *    equil_freqs contains the equilibrium frequencies for the model.
 ****************************************************************************/ 
EVOMODEL_T* make_f81_model(
  const double rate,
  ARRAY_T* equil_freqs
);

/****************************************************************************
 *  Allocate a data structure for an F84 model
 *    rate is overall mutation rate
 *    equil_freqs contains the equilibrium frequencies for the model.
 *    transition_transversion is transition/transversion rate ratio
 ****************************************************************************/ 
EVOMODEL_T* make_f84_model(
  const double rate,
  ARRAY_T* equil_freqs,
  const double transition_transversion
);

/****************************************************************************
 *  Allocate a data structure for an HKY model
 *    rate is overall mutation rate
 *    equil_freqs contains the equilibrium frequencies for the model.
 *    transition_transversion transition/transversion rate ratio
 ****************************************************************************/ 
EVOMODEL_T* make_hky_model(
  const double rate,
  ARRAY_T* equil_freqs,
  const double transition_transversion
);

/****************************************************************************
 *  Allocate a data structure for a Tamura-Nei model
 *    rate is overall mutation rate
 *    equil_freqs contains the equilibrium frequencies for the model.
 *    transition_transversion transition/transversion rate ratio
 *    pyrimidine_purine R/R to Y/Y transitions rate ratio.
 ****************************************************************************/ 
EVOMODEL_T* make_tamura_nei_model(
  const double rate,
  ARRAY_T* equil_freqs,
  const double transition_transversion,
  const double pyrimidine_purine
);

/****************************************************************************
 *  Turn model into Halpern-Bruno model
 ****************************************************************************/
void make_halpern_bruno_model(EVOMODEL_T* model, ARRAY_T* site_specific_freqs);

/****************************************************************************
 *  Free the data structure for an evolutionary model
 ****************************************************************************/
void free_model(EVOMODEL_T* model);

/****************************************************************************
 *  Return the equilibrium frequency for a base in a model.
 *****************************************************************************/ 
double get_model_equilibrium_freq(
  const char base, 
  const EVOMODEL_T* model
);
/****************************************************************************
   Print a evolutionary model.
  ****************************************************************************/
void print_evomodel(EVOMODEL_T* model);
#endif
