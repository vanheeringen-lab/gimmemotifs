/****************************************************************************
 * FILE: evomodel.c
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 11/18/2004
 * PROJECT: EVOMCAST
 * DESCRIPTION: Evolutionary models for nucleotide substitutions
 * COPYRIGHT: 2004, UW
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

// Nucleotide alphabest order as in motif.h
extern char alphabet[];

static char * model_names[] = { 
  "single",
  "average",
  "Jukes-Cantor",
  "Kimura 2 parameter",
  "F81",
  "F84",
  "HKY",
  "Tamura-Nei"
};

static  MATRIX_T* make_rate_matrix(
  ARRAY_T* equil_freqs,
  const double non_specific_rate,
  const double pur_pur_rate,
  const double py_py_rate
);

static double jc_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *model
);

static double k2_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *model
);

static double f81_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *model
);

static double f84_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *model
);

static double hky_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *model
);

static double tamura_nei_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *model
);

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
   Print a evolutionary model.
  ****************************************************************************/
void print_evomodel(EVOMODEL_T* model)
{
  print_array(model->equil_freqs, 8, 4, TRUE, stdout);
  print_matrix(model->rate_matrix, 8, 4, TRUE, stdout);
}

/*******************************************************************
 * Build an array of evolutionary models for the positions of a motif. 
 * The background model is stored in the first element of the array.
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
) 
{
  int alph_size = get_alph_size(ALPH_SIZE);
  int num_models = motif->length + 1;
  EVOMODEL_T** models = mm_malloc(num_models * sizeof(EVOMODEL_T*));

  // Make the background model and save it as models[0].
  // We never use HB for the background model.
  // AVG_MODEL uses JC_MODEL to compute column frequencies (for p-values).
  EVOMODEL_T* model = make_model(
    (model_type == AVERAGE_MODEL) ? JC_MODEL : model_type,
    bg_rate,
    transition_transversion,
    purine_pyrimidine,
    bg_freqs,
    NULL
  );
  int i = 0;
  models[i] = model;
  
  // print_evomodel(model);

  for (i = 1; i < num_models; i++) {
    // Use motif PSFM for equilibrium freqs. for model.
    ARRAY_T* site_specific_freqs = allocate_array(alph_size);
    int j = 0;
    for(j = 0; j < alph_size; j++) {
      double value = get_matrix_cell(i - 1, j, motif->freqs);
      set_array_item(j, value, site_specific_freqs);
    }
    if (use_halpern_bruno == FALSE) {
      // If not using Halpern-Bruno use site specific freqs for 
      // equilibrium freqs.
      model = make_model(
        model_type, 
        fg_rate, 
        transition_transversion,
        purine_pyrimidine,
        site_specific_freqs,
        NULL
      );
    }
    else {
      // If using Halpern-Bruno use background freqs for 
      // equilibrium freqs and pass on site specific freqs.
      model = make_model(
        model_type, 
        fg_rate, 
        transition_transversion,
        purine_pyrimidine,
        bg_freqs,
        site_specific_freqs
      );
    }
    models[i] = model;
    //FIXME: print_evomodel(model);
    free_array(site_specific_freqs);
  }

  return models;
}

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
 * If an array of site-specific frequencies is provided the Halpern-Bruno
 * modification of the model will be used.
 ********************************************************************/
EVOMODEL_T* make_model(
  const MODEL_TYPE_T model_type, 
  const double rate,
  const double transition_transversion,
  const double purine_pyrimidine,
  ARRAY_T* equil_freqs,
  ARRAY_T* site_specific_freqs
) 
{

  EVOMODEL_T* model = NULL;

  switch(model_type) {
    case SINGLE_MODEL:
      if (site_specific_freqs != NULL) {
        die("Single model doesn't support Halpern-Bruno\n");
      }
      model = make_psfm_model(equil_freqs);
      break;
    case AVERAGE_MODEL:
      if (site_specific_freqs != NULL) {
        die("Average model doesn't support Halpern-Bruno\n");
      }
      model = make_avg_psfm_model(equil_freqs);
      break;
    case JC_MODEL:
      model = make_jc_model(rate, equil_freqs);
      break;
    case K2_MODEL:
      model = make_k2_model(rate, equil_freqs, transition_transversion);
      break;
    case F81_MODEL:
      model = make_f81_model(rate, equil_freqs);
      break;
    case F84_MODEL:
      model = make_f84_model(rate, equil_freqs, transition_transversion);
      break;
    case HKY_MODEL:
      model = make_hky_model(rate, 
                             equil_freqs,
                             transition_transversion);
      break;
    case TAMURA_NEI_MODEL:
      model = make_tamura_nei_model(rate, 
                                    equil_freqs, 
				    transition_transversion, 
                                    purine_pyrimidine);
      break;
  }

  if (site_specific_freqs != NULL) {
    make_halpern_bruno_model(model, site_specific_freqs);
  }
  return model;
}

/****************************************************************************
  * Allocate a data structure for Position Specific Site Scoring Model
  ****************************************************************************/
EVOMODEL_T* make_psfm_model(ARRAY_T* equil_freqs) {

  EVOMODEL_T* model = mm_malloc(sizeof(EVOMODEL_T));
  model->type = SINGLE_MODEL;
  model->equil_freqs = allocate_array(get_array_length(equil_freqs));
  copy_array(equil_freqs, model->equil_freqs);
  model->use_halpern_bruno = FALSE;

  // No parameters for the SINGLE_MODEL
  model->num_params = 0;
  model->param_names =  NULL;
  model->param_values = NULL;
  model->rate_matrix = NULL;

  return model;
}

/****************************************************************************
  * Allocate a data structure for an Average Position Specific Site 
  * Scoring Model
  ****************************************************************************/
EVOMODEL_T* make_avg_psfm_model (ARRAY_T* equil_freqs) {

  EVOMODEL_T* model = NULL;

  model = mm_malloc(sizeof(EVOMODEL_T));
  model->type = AVERAGE_MODEL;
  model->equil_freqs = allocate_array(get_array_length(equil_freqs));
  copy_array(equil_freqs, model->equil_freqs);
  model->use_halpern_bruno = FALSE;

  // No parameters for the AVERAGE_MODEL
  model->num_params = 0;
  model->param_names =  NULL;
  model->param_values = NULL;
  model->rate_matrix = NULL;

  return model;
}

/****************************************************************************
 * Allocate a data structure for the Jukes-Cantor model
 *   rate is overall substitution rate
 *   equil_freqs contains the equilibrium frequencies for the model.
 ****************************************************************************/
EVOMODEL_T* make_jc_model (
  const double rate,
  ARRAY_T* equil_freqs
) 
{

  int i = 0;
  EVOMODEL_T* model = NULL;
  ARRAY_T* jc_freqs = NULL;

  model = mm_malloc(sizeof(EVOMODEL_T));
  model->type = JC_MODEL;
  model->equil_freqs = allocate_array(get_array_length(equil_freqs));
  copy_array(equil_freqs, model->equil_freqs);
  model->use_halpern_bruno = FALSE;
  model->evaluator = jc_evaluator;

  // The only parameter for Jukes-Cantor is the subsitution substitution_rate
  model->num_params = 1;
  model->param_names = (char **) mm_malloc(model->num_params * sizeof(char *));
  *(model->param_names) = "overall substitution rate";
  model->param_values = (double *) mm_malloc(model->num_params * sizeof(double));
  (model->param_values)[0] = rate;

  // JC doesn't use the equil freqs, but assumes all transitions
  // equally likely
  jc_freqs = allocate_array(4);
  for(i = 0; i < 4; i++) {
    set_array_item(i, 0.25, jc_freqs);
  }

  // Make rate matrix from equivalent Tamura-Nei model.
  model->rate_matrix = make_rate_matrix(jc_freqs, 4 * rate / 3, 0.0, 0.0);
  free_array(jc_freqs);

  return model;
}


/****************************************************************************
 *  Allocate a data structure for the Kimura 2 parameter model
 *    rate is overall substitution rate
 *    equil_freqs contains the equilibrium frequencies for the model.
 *    transition_transversion is transition/transversion rate ratio
 ****************************************************************************/
EVOMODEL_T* make_k2_model(
  const double rate, 
  ARRAY_T* equil_freqs, 
  const double transition_transversion
) 
{
  EVOMODEL_T* model = mm_malloc(sizeof(EVOMODEL_T));
  model->type = K2_MODEL;
  model->equil_freqs = allocate_array(get_array_length(equil_freqs));
  copy_array(equil_freqs, model->equil_freqs);
  model->use_halpern_bruno = FALSE;
  model->evaluator = k2_evaluator;

  // The only parameters for the Kimura 2 parameter model
  // are the overal substitution rate transition transversion ratio.
  model->num_params = 2;
  model->param_names = (char **) mm_malloc(model->num_params * sizeof(char *));
  (model->param_names)[0] = "overall substittution rate";
  (model->param_names)[1] = "transition/transversion ratio";
  model->param_values = (double *) mm_malloc(model->num_params * sizeof(double));
  (model->param_values)[0] = rate;
  (model->param_values)[1] = transition_transversion;

  // K2 doesn't use the equil freqs, but assumes all transitions equally likely
  ARRAY_T* k2_freqs = allocate_array(4);
  int i;
  for(i = 0; i < 4; i++) {
    set_array_item(i, 0.25, k2_freqs);
  }

  // Calculate equivalent parameters for Tamura-Nei model.
  double Pa = 0.25;
  double Pc = 0.25;
  double Pg = 0.25;
  double Pt = 0.25;
  double Pr = Pa + Pg;
  double Py = Pc + Pt;
  double pyr_pyr_rate = rate * ((Pr * Py * transition_transversion) 
    - (Pa * Pg) - (Pc * Pt)) / (2.0 * (1.0 + transition_transversion)
	  * ((Py * Pa * Pg) + (Pr * Pc * Pt)));
  double pur_pur_rate = pyr_pyr_rate;
  double non_specific_rate = rate 
    / ( 2.0 * Pr * Py * (1.0 + transition_transversion));
  // Build rate matrix based on equivalent Tamura-Nei model
  model->rate_matrix = make_rate_matrix(
    k2_freqs,
    non_specific_rate,
    pur_pur_rate, 
    pyr_pyr_rate
  );
  free_array(k2_freqs);

  return model;
}

/****************************************************************************
 *  Allocate a data structure for an F81 model
 *    rate is overall substitution rate
 *    equil_freqs contains the equilibrium frequencies for the model.
 ****************************************************************************/ 
EVOMODEL_T* make_f81_model(
  const double rate,
  ARRAY_T* equil_freqs
) 
{
  EVOMODEL_T* model = mm_malloc(sizeof(EVOMODEL_T));
  model->type = F81_MODEL;
  model->equil_freqs = allocate_array(get_array_length(equil_freqs));
  copy_array(equil_freqs, model->equil_freqs);
  model->use_halpern_bruno = FALSE;
  model->evaluator = f81_evaluator;

  // The parameters for the F81 model are the equilibrium base frequencies.
  model->num_params = 5;
  model->param_names = (char **) mm_malloc(model->num_params * sizeof(char *));
  (model->param_names)[0] = "overall substitution rate";
  (model->param_names)[1] = "A frequency";
  (model->param_names)[2] = "C frequency";
  (model->param_names)[3] = "G frequency";
  (model->param_names)[4] = "T frequency";
  model->param_values = (double *) mm_malloc(model->num_params * sizeof(double));
  (model->param_values)[0] = rate;
  double Pa = get_array_item(alphabet_index('A', alphabet), model->equil_freqs);
  double Pc = get_array_item(alphabet_index('C', alphabet), model->equil_freqs);
  double Pg = get_array_item(alphabet_index('G', alphabet), model->equil_freqs);
  double Pt = get_array_item(alphabet_index('T', alphabet), model->equil_freqs);
  double Pr = Pa + Pg;
  double Py = Pc + Pt;
  (model->param_values)[1] = Pa;
  (model->param_values)[2] = Pc;
  (model->param_values)[3] = Pg;
  (model->param_values)[4] = Pt;

  // Calculate equivalent parameters for Tamura-Nei model.
  double pyr_pyr_rate = 0;
  double pur_pur_rate = 0;
  double non_specific_rate = rate / ( 3.0 * Pr * Py);

  // Build rate matrix from equivalent Tamura-Nei model
  model->rate_matrix = make_rate_matrix(
    model->equil_freqs, 
    non_specific_rate,
    pyr_pyr_rate, 
    pur_pur_rate 
  );

  return model;
}

/****************************************************************************
 *  Allocate a data structure for an F84 model
 *    rate is overall substitution rate
 *    equil_freqs contains the equilibrium frequencies for the model.
 *    transition_transversion is transition/transversion rate ratio
 ****************************************************************************/ 
EVOMODEL_T* make_f84_model(
  const double rate,
  ARRAY_T* equil_freqs,
  const double transition_transversion
) 
{
  EVOMODEL_T* model = mm_malloc(sizeof(EVOMODEL_T));
  model->type = F84_MODEL;
  model->equil_freqs = allocate_array(get_array_length(equil_freqs));
  copy_array(equil_freqs, model->equil_freqs);
  model->use_halpern_bruno = FALSE;
  model->evaluator = f84_evaluator;

  // The parameters for the F84 model are the equilibrium 
  // base frequencies and the transition/transversion ratio.
  model->num_params = 6;
  model->param_names = (char **) mm_malloc(model->num_params * sizeof(char *));
  (model->param_names)[0] = "overall substitution rate";
  (model->param_names)[1] = "A frequency";
  (model->param_names)[2] = "C frequency";
  (model->param_names)[3] = "G frequency";
  (model->param_names)[4] = "T frequency";
  (model->param_names)[5] = "ratio of transition/transversion rates ";
  model->param_values = (double *) mm_malloc(model->num_params * sizeof(double));
  (model->param_values)[0] = rate;
  double Pa = get_array_item(alphabet_index('A', alphabet), model->equil_freqs);
  double Pc = get_array_item(alphabet_index('C', alphabet), model->equil_freqs);
  double Pg = get_array_item(alphabet_index('G', alphabet), model->equil_freqs);
  double Pt = get_array_item(alphabet_index('T', alphabet), model->equil_freqs);
  (model->param_values)[1] = Pa;
  (model->param_values)[2] = Pc;
  (model->param_values)[3] = Pg;
  (model->param_values)[4] = Pt;
  (model->param_values)[5] = transition_transversion;

  // Calculate equivalent parameters for Tamura-Nei model.
  double Pr = Pa + Pg;
  double Py = Pc + Pt;
  double pyr_pyr_rate = rate * ((Pr * Py * transition_transversion) 
    - (Pa * Pg) - (Pc * Pt)) / (2.0 * (1.0 + transition_transversion)
	  * ((Py * Pa * Pg) + (Pr * Pc * Pt)));
  double pur_pur_rate = pyr_pyr_rate;
  double non_specific_rate = 
    rate / ( 2.0 * Pr * Py * (1.0 + transition_transversion));
  model->rate_matrix = make_rate_matrix(
    model->equil_freqs,
    non_specific_rate,
    pur_pur_rate, 
    pyr_pyr_rate
  );

  return model;
}

/****************************************************************************
 *  Allocate a data structure for an HKY model
 *    rate is overall substitution rate
 *    equil_freqs contains the equilibrium frequencies for the model.
 *    transition_transversion is transition/transversion rate ratio
 ****************************************************************************/ 
EVOMODEL_T* make_hky_model(
  const double rate,
  ARRAY_T* equil_freqs,
  const double transition_transversion
) 
{
  EVOMODEL_T* model = mm_malloc(sizeof(EVOMODEL_T));
  model->type = HKY_MODEL;
  model->equil_freqs = allocate_array(get_array_length(equil_freqs));
  copy_array(equil_freqs, model->equil_freqs);
  model->use_halpern_bruno = FALSE;
  model->evaluator = hky_evaluator;

  // The parameters for the HKY model are the equilibrium 
  // base frequencies and the transition/transversion ratio.
  model->num_params = 6;
  model->param_names = 
    (char **) mm_malloc(model->num_params * sizeof(char *));
  (model->param_names)[0] = "overall substitution rate";
  (model->param_names)[1] = "A frequency";
  (model->param_names)[2] = "C frequency";
  (model->param_names)[3] = "G frequency";
  (model->param_names)[4] = "T frequency";
  (model->param_names)[5] = "ratio of transition/transversion substitution"
    " rates";
  model->param_values = (double *) mm_malloc(model->num_params * sizeof(double));
  (model->param_values)[0] = rate;
  double Pa = get_array_item(alphabet_index('A', alphabet), model->equil_freqs);
  double Pc = get_array_item(alphabet_index('C', alphabet), model->equil_freqs);
  double Pg = get_array_item(alphabet_index('G', alphabet), model->equil_freqs);
  double Pt = get_array_item(alphabet_index('T', alphabet), model->equil_freqs);
  (model->param_values)[1] = Pa;
  (model->param_values)[2] = Pc;
  (model->param_values)[3] = Pg;
  (model->param_values)[4] = Pt;
  (model->param_values)[5] = transition_transversion;

  // Calculate equivalent parameters for Tamura-Nei model.
  double Pr = Pa + Pg;
  double Py = Pc + Pt;
  double pyr_pyr_rate = rate * ((Pr * Py * transition_transversion) 
    - (Pa * Pg) - (Pc * Pt)) / (2.0 * (1.0 + transition_transversion)
	  * ((Py * Pa * Pg * Pr / Py) + (Pr * Pc * Pt)));
  double pur_pur_rate = pyr_pyr_rate * Pr / Py;
  double non_specific_rate = 
    rate / ( 2.0 * Pr * Py * (1.0 + transition_transversion));
  model->rate_matrix = make_rate_matrix(
    model->equil_freqs,
    non_specific_rate,
    pur_pur_rate, 
    pyr_pyr_rate
  );

  return model;
}

/****************************************************************************
 *  Allocate a data structure for a Tamura-Nei model
 *    rate overall substitution rate
 *    transition_transversion transition/transversion rate ratio
 *    purine_pyrimidine purine/purine to pyrimidine/pyrimidine 
 *                      transition rate ratios
 ****************************************************************************/ 
EVOMODEL_T* make_tamura_nei_model(
  const double rate,
  ARRAY_T* equil_freqs,
  const double transition_transversion,
  const double purine_pyrimidine) 
{
  EVOMODEL_T* model = mm_malloc(sizeof(EVOMODEL_T));
  model->type = TAMURA_NEI_MODEL;
  model->equil_freqs = allocate_array(get_array_length(equil_freqs));
  copy_array(equil_freqs, model->equil_freqs);
  model->use_halpern_bruno = FALSE;
  model->evaluator = tamura_nei_evaluator;

  // The parameters for the Tamura-Nei model are the equilibrium 
  // base frequencies, the transition/transversion ratio, and the
  // purine-purine/pyrimidine-pyrimidine rate ratio.
  model->num_params = 7;
  model->param_names = (char **) mm_malloc(model->num_params * sizeof(char *));
  (model->param_names)[0] = "overall substitution rate";
  (model->param_names)[1] = "A frequency";
  (model->param_names)[2] = "C frequency";
  (model->param_names)[3] = "G frequency";
  (model->param_names)[4] = "T frequency";
  (model->param_names)[5] = "ratio of transition/transversion substitution"
    " rates";
  (model->param_names)[6] = "ratio of R/R to Y/Y transitions";
  model->param_values = (double *) mm_malloc(model->num_params * sizeof(double));
  (model->param_values)[0] = rate;
  double Pa = get_array_item(alphabet_index('A', alphabet), model->equil_freqs);
  double Pc = get_array_item(alphabet_index('C', alphabet), model->equil_freqs);
  double Pg = get_array_item(alphabet_index('G', alphabet), model->equil_freqs);
  double Pt = get_array_item(alphabet_index('T', alphabet), model->equil_freqs);
  (model->param_values)[1] = Pa;
  (model->param_values)[2] = Pc;
  (model->param_values)[3] = Pg;
  (model->param_values)[4] = Pt;
  (model->param_values)[5] = transition_transversion;
  (model->param_values)[6] = purine_pyrimidine;

  // Calculate equivalent parameters for Tamura-Nei model.
  double Pr = Pa + Pg;
  double Py = Pc + Pt;
  double pyr_pyr_rate = 
    rate * ((Pr * Py * transition_transversion) - (Pa * Pg) - (Pc * Pt))
            / (2.0 * (1.0 + transition_transversion)
	            * ((Py * Pa * Pg) + (Pr * Pc * Pt)));
  double pur_pur_rate = purine_pyrimidine * pyr_pyr_rate;
  double non_specific_rate = 
    rate / ( 2.0 * Pr * Py * (1.0 + transition_transversion));
  model->rate_matrix = make_rate_matrix(
    model->equil_freqs,
    non_specific_rate,
    pur_pur_rate, 
    pyr_pyr_rate
  );

  return model;
}

/****************************************************************************
 *  Turn model into Halpern-Bruno model
 ****************************************************************************/
void make_halpern_bruno_model(EVOMODEL_T* model, ARRAY_T* site_specific_freqs) {

  assert(model != NULL);
  assert(site_specific_freqs != NULL);

  const double epsilon = 1e-6;
  MATRIX_T* rate_matrix = model->rate_matrix;
  int num_rows = get_num_rows(rate_matrix);
  int num_cols = get_num_cols(rate_matrix);
  assert(num_rows == num_cols);
  MATRIX_T* hb_matrix = allocate_matrix(num_rows, num_cols);

  int i;
  for(i = 0; i < num_rows; i++) {
    int j;
    set_matrix_cell(i, i, 0, hb_matrix); 
    for(j = 0; j < num_cols; j++) {
      double value;
      if (i != j) {
        double Qij = get_matrix_cell(i, j, rate_matrix);
        double Qji = get_matrix_cell(j, i, rate_matrix);
        double fi = get_array_item(i, site_specific_freqs);
        double fj = get_array_item(j, site_specific_freqs);
        double numerator = Qij * log((fj * Qji) / (fi * Qij));
        double denominator = (1.0 - (fi * Qij) / (fj * Qji));
        if (fabs(fi * Qij - fj *Qji) < epsilon) {
          // Using L'Hopital's rule
          value = Qij;
        }
        else {
          value = numerator / denominator;
        }
        set_matrix_cell(i, j, value, hb_matrix); 
        // Set diagonal elments so rows sum to 0
        set_matrix_cell(
          i, 
          i, 
          get_matrix_cell(i, i, hb_matrix) - value,
          hb_matrix
        ); 
      }
    }
  }
  model->rate_matrix = hb_matrix;
  model->use_halpern_bruno = TRUE;
  // Now replace the equilibrium frequecies with the site-specific ones
  // so that pruning algorithm will work.
  copy_array(site_specific_freqs, model->equil_freqs);
  free_matrix(rate_matrix);
}

/****************************************************************************
  * Free a data structure for an evolutionary model
  ****************************************************************************/
void free_model(EVOMODEL_T* model) {

  if (model == NULL) {
    return;
  }

  free_array(model->equil_freqs);
  myfree(model->param_names);
  myfree(model->param_values);
  free_matrix(model->rate_matrix);
  myfree(model);
}

/****************************************************************************
  * get methods for evomodel data structure
  ****************************************************************************/
const char* get_model_name(const EVOMODEL_T* model) {

  assert(model != NULL);

  return model_names[model->type];
}

MODEL_TYPE_T get_model_type(const EVOMODEL_T* model) {

  assert(model != NULL);
  
  return model->type;
}

int get_num_model_params (const EVOMODEL_T* model) {

  assert(model != NULL);

  return model->num_params;
}

char** get_model_param_names (const EVOMODEL_T* model) {

  assert(model != NULL);

  return model->param_names;
}

double* get_model_param_values(const EVOMODEL_T* model) {

  assert(model != NULL);

  return model->param_values;
}

double get_model_equil_freq(char base, const EVOMODEL_T* model) {

  double value = NaN();

  assert(model != NULL);
  assert(strchr(alphabet, base));

  base = toupper(base);
  value = get_array_item(alphabet_index(base, alphabet), model->equil_freqs);

  return value;
};

/****************************************************************************
 *  Return the array of equilibrium frequencies used by the model.  
 *  Caller is responsible for freeing the returned array.
 ****************************************************************************/ 
ARRAY_T* get_model_equil_freqs(const EVOMODEL_T* model) {

  assert(model != NULL);

  ARRAY_T* equil_freqs = allocate_array(get_array_length(model->equil_freqs));
  copy_array(model->equil_freqs, equil_freqs);

  return equil_freqs;
}

MATRIX_T* get_model_rate_matrix(EVOMODEL_T* model) {

  assert(model != NULL);
  MATRIX_T* rate_matrix = allocate_matrix(
    get_num_rows(model->rate_matrix),
    get_num_cols(model->rate_matrix)
  );
  copy_matrix(model->rate_matrix, rate_matrix);

  return rate_matrix;
}

MATRIX_T* get_model_prob_matrix(double t, const EVOMODEL_T* model) {

  assert(model != NULL);

  MATRIX_T* prob_matrix = NULL;

  if (model->use_halpern_bruno == TRUE) {
    prob_matrix = matrix_exponential(t, model->rate_matrix);
  } else {
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
        double value = model->evaluator(old_base, new_base, t, model);
        set_matrix_cell(i, j, value, prob_matrix);
      }
    }
  }
  return prob_matrix;
}

MATRIX_T* old_get_model_prob_matrix(const EVOMODEL_T* model, double t) {

  assert(model != NULL);

  MATRIX_T* prob_matrix = NULL;

  if (model->use_halpern_bruno == TRUE) {
    prob_matrix = matrix_exponential(t, model->rate_matrix);
  } else {
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
        double value = model->evaluator(old_base, new_base, t, model);
        set_matrix_cell(i, j, value, prob_matrix);
      }
    }
  }
  return prob_matrix;
}

BOOLEAN_T uses_halpern_bruno(const EVOMODEL_T* model) {

  assert(model != NULL);

  return model->use_halpern_bruno;
}

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

/****************************************************************************
  * Calculate probability for nucleotide substitution using Jukes-Cantor
  * model. See for example Durbin et al, "Biological Sequence Analyis",
  * page 195.
  ****************************************************************************/
static double jc_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *model
)
{
  assert(model != NULL);
  assert(JC_MODEL == model->type);
  assert(model->num_params == 1);
  assert(model->param_names != NULL);
  assert(model->param_values != NULL);

  double rate = (model->param_values)[0];
  double Pa = 0.25;
  double Pc = 0.25;
  double Pg = 0.25;
  double Pt = 0.25;

  double pyr_pyr_rate = 0.0;
  double pur_pur_rate = 0.0;
  double non_specific_rate = 4.0 * rate / 3.0;

  return inner_tamura_nei_evaluator(
    old_base, 
    new_base, 
    t, 
    Pa, 
    Pc, 
    Pg, 
    Pt, 
    non_specific_rate,
    pur_pur_rate, 
    pyr_pyr_rate 
  );
}

/****************************************************************************
  * Calculate probability for nucleotide substitution using Kimura 2 parameter 
  * model. See for example Durbin et al, "Biological Sequence Analyis",
  * page 196.
  ****************************************************************************/
static double k2_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *model
)
{
  assert(model != NULL);
  assert(K2_MODEL == model->type);
  assert(model->num_params == 2);
  assert(model->param_names != NULL);
  assert(model->param_values != NULL);

  double rate = (model->param_values)[0];
  double transition_transversion = (model->param_values)[1];
  double Pa = 0.25;
  double Pc = 0.25;
  double Pg = 0.25;
  double Pt = 0.25;
  double Pr = Pa + Pg;
  double Py = Pc + Pt;

  double pyr_pyr_rate = rate * ((Pr * Py * transition_transversion) 
    - (Pa * Pg) - (Pc * Pt)) / (2.0 * (1.0 + transition_transversion)
	  * ((Py * Pa * Pg) + (Pr * Pc * Pt)));
  double pur_pur_rate = pyr_pyr_rate;
  double non_specific_rate = rate 
    / ( 2.0 * Pr * Py * (1.0 + transition_transversion));

  return inner_tamura_nei_evaluator(
    old_base,
    new_base,
    t,
    Pa,
    Pc,
    Pg,
    Pt,
    non_specific_rate,
    pur_pur_rate,
    pyr_pyr_rate
  );
}

/****************************************************************************
  * Calculate probability for nucleotide substitution using F81
  * model. See for example Ewens, "Statistical Methods in Bioinformatcs"
  * page 382
  ****************************************************************************/
static double f81_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *model
)
{
  assert(model != NULL);
  assert(F81_MODEL == model->type);
  assert(model->num_params == 5);
  assert(model->param_names != NULL);
  assert(model->param_values != NULL);

  double rate = (model->param_values)[0];
  double Pa = model->param_values[1];
  double Pc = model->param_values[2];
  double Pg = model->param_values[3];
  double Pt = model->param_values[4];
  double Pr = Pa + Pg;
  double Py = Pc + Pt;

  double pyr_pyr_rate = 0;
  double pur_pur_rate = 0;
  double non_specific_rate = rate / (3.0 * Pr * Py);

  return inner_tamura_nei_evaluator(
    old_base, 
    new_base, 
    t, 
    Pa, 
    Pc, 
    Pg, 
    Pt, 
    non_specific_rate,
    pur_pur_rate, 
    pyr_pyr_rate 
  );
}



/****************************************************************************
  * Calculate probability for nucleotide substitution using F84
  * model. This is the Tamura-Nei model with r = 1.
  * See Felsenstein, "Inferring Phylogenies" page 201.
  ****************************************************************************/
static double f84_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *model
)
{
  assert(model != NULL);
  assert(F84_MODEL == model->type);
  assert(model->num_params == 6);
  assert(model->param_names != NULL);
  assert(model->param_values != NULL);

  double rate = (model->param_values)[0];
  double Pa = (model->param_values)[1];
  double Pc = (model->param_values)[2];
  double Pg = (model->param_values)[3];
  double Pt = (model->param_values)[4];
  double Pr = Pa + Pg;
  double Py = Pc + Pt;
  double transition_transversion = (model->param_values)[5];

  assert(Py > 0.0);
  assert(Pr > 0.0);

  double pyr_pyr_rate = rate * ((Pr * Py * transition_transversion) 
    - (Pa * Pg) - (Pc * Pt)) / (2.0 * (1.0 + transition_transversion)
	  * ((Py * Pa * Pg) + (Pr * Pc * Pt)));
  double pur_pur_rate = pyr_pyr_rate;
  double non_specific_rate = 
    rate / ( 2.0 * Pr * Py * (1.0 + transition_transversion));

  return inner_tamura_nei_evaluator(
    old_base, 
    new_base, 
    t, 
    Pa, 
    Pc, 
    Pg, 
    Pt, 
    non_specific_rate,
    pur_pur_rate, 
    pyr_pyr_rate
  );
}

/****************************************************************************
  * Calculate probability for nucleotide substitution using HKY
  * model. This is the Tamura-Nei model with r = Pr/Py.
  * See Felsenstein, "Inferring Phylogenies" page 201.
  ****************************************************************************/
static double hky_evaluator(
  const char old_base, 
  const char new_base, 
  const double t,
  const EVOMODEL_T *model
)
{
  assert(model != NULL);
  assert(HKY_MODEL == model->type);
  assert(model->num_params == 6);
  assert(model->param_names != NULL);
  assert(model->param_values != NULL);

  double rate = (model->param_values)[0];
  double Pa = (model->param_values)[1];
  double Pc = (model->param_values)[2];
  double Pg = (model->param_values)[3];
  double Pt = (model->param_values)[4];
  double transition_transversion = (model->param_values)[5];
  double Pr = Pa + Pg;
  double Py = Pc + Pt;

  assert(Pr > 0.0);
  assert(Py > 0.0);

  double pyr_pyr_rate = rate * ((Pr * Py * transition_transversion) 
    - (Pa * Pg) - (Pc * Pt)) / (2.0 * (1.0 + transition_transversion)
	  * ((Py * Pa * Pg * Pr / Py) + (Pr * Pc * Pt)));
  double pur_pur_rate = pyr_pyr_rate * Pr / Py;
  double non_specific_rate = 
    rate / ( 2.0 * Pr * Py * (1.0 + transition_transversion));

  return inner_tamura_nei_evaluator(
    old_base, 
    new_base, 
    t, 
    Pa, 
    Pc, 
    Pg, 
    Pt, 
    non_specific_rate,
    pur_pur_rate, 
    pyr_pyr_rate
  );
}

/****************************************************************************
  * Calculate probability for nucleotide substitution using Tamura-Nei
  * model. See for example Felsenstein, "Inferring Phylogenies",
  * page 203.
  ****************************************************************************/
static double tamura_nei_evaluator(
  const char new_base, 
  const char old_base, 
  const double t,
  const EVOMODEL_T *model
)
{
  assert(model != NULL);
  assert(TAMURA_NEI_MODEL == model->type);
  assert(model->num_params == 7);
  assert(model->param_names != NULL);
  assert(model->param_values != NULL);

  double rate = (model->param_values)[0];
  double Pa = (model->param_values)[1];
  double Pc = (model->param_values)[2];
  double Pg = (model->param_values)[3];
  double Pt = (model->param_values)[4];
  double transition_transversion = (model->param_values)[5];
  double purine_pyrimidine = (model->param_values)[6];
  double Pr = Pa + Pg;
  double Py = Pc + Pt;

  assert(Pr > 0.0);
  assert(Py > 0.0);

  double pyr_pyr_rate = rate * ((Pr * Py * transition_transversion) 
    - (Pa * Pg) - (Pc * Pt)) / (2.0 * (1.0 + transition_transversion)
	  * ((Py * Pa * Pg * purine_pyrimidine) + (Pr * Pc * Pt)));
  double pur_pur_rate = purine_pyrimidine * pyr_pyr_rate;
  double non_specific_rate 
    = rate / ( 2.0 * Pr * Py * (1.0 + transition_transversion));

  return inner_tamura_nei_evaluator(
    new_base, 
    old_base, 
    t,
    Pa, 
    Pc, 
    Pg, 
    Pt, 
    non_specific_rate,
    pur_pur_rate, 
    pyr_pyr_rate
  );

}

/****************************************************************************
  * Calculate the instantaneous transition rate for nucleotide substitution 
  * using Tamura-Nei model. See for example Felsenstein, 
  * "Inferring Phylogenies", page 201.
  ****************************************************************************/
static  MATRIX_T* make_rate_matrix(
  ARRAY_T* equil_freqs,
  const double B,  // non-specific rate
  const double Ar, // purine-purine rate
  const double Ay  // pyrimidine-pyrimidine rate
) 
{
  // Turn parameters into intermediate values
  // used by model
  double Pa = get_array_item(alphabet_index('A', alphabet), equil_freqs);
  double Pc = get_array_item(alphabet_index('C', alphabet), equil_freqs);
  double Pg = get_array_item(alphabet_index('G', alphabet), equil_freqs);
  double Pt = get_array_item(alphabet_index('T', alphabet), equil_freqs);
  double Pr = Pa + Pg;
  double Py = Pc + Pt;

  assert(Pr > 0.0);
  assert(Py > 0.0);

  MATRIX_T* rate_matrix = allocate_matrix(4, 4);
  set_matrix_cell(0, 0, - B * Py - (Ar / Pr + B) * Pg, rate_matrix);
  set_matrix_cell(0, 1, B * Pc, rate_matrix);
  set_matrix_cell(0, 2, (Ar / Pr + B) * Pg, rate_matrix);
  set_matrix_cell(0, 3, B * Pt, rate_matrix);
  set_matrix_cell(1, 0, B * Pa, rate_matrix);
  set_matrix_cell(1, 1, -B * Pr - (Ay / Py + B) * Pt, rate_matrix);
  set_matrix_cell(1, 2, B * Pg, rate_matrix);
  set_matrix_cell(1, 3, (Ay / Py + B) * Pt, rate_matrix);
  set_matrix_cell(2, 0, (Ar / Pr + B) * Pa, rate_matrix);
  set_matrix_cell(2, 1, B * Pc, rate_matrix);
  set_matrix_cell(2, 2, -B * Py - (Ar / Pr + B) * Pa, rate_matrix);
  set_matrix_cell(2, 3, B * Pt, rate_matrix);
  set_matrix_cell(3, 0, B * Pa, rate_matrix);
  set_matrix_cell(3, 1, (Ay / Py + B) * Pc, rate_matrix);
  set_matrix_cell(3, 2, B * Pg, rate_matrix);
  set_matrix_cell(3, 3, -B * Pr - (Ay / Py + B) * Pc, rate_matrix);

  return rate_matrix;
  
}

/******************************************************************************
 * This function calculates the exponential of a square matrix
 * using the Taylor series expansion for the exponential:
    
 * exp(At) = I + t A + t^2/2 A^2 + t^3/3! A^3 + ...
******************************************************************************/
MATRIX_T* matrix_exponential(const double t, MATRIX_T *a) {

  assert(a != NULL);

  double id_raw_matrix[4][4] = { 
    { 1.0, 0.0, 0.0, 0.0},
    { 0.0, 1.0, 0.0, 0.0},
    { 0.0, 0.0, 1.0, 0.0},
    { 0.0, 0.0, 0.0, 1.0}
  };

  int num_rows = get_num_rows(a);
  int num_cols = get_num_cols(a);
  assert(num_rows == num_cols);

  MATRIX_T* temp1 = allocate_matrix(num_rows, num_cols);
  MATRIX_T* temp2 = NULL;
  MATRIX_T* result = allocate_matrix(num_rows, num_cols);
  fill_matrix((double *) &id_raw_matrix, temp1);
  fill_matrix((double *) &id_raw_matrix, result);
  int i;
  const int max_iterations = 50;
  // Taylor series for matrix exponential
  for (i = 1; i <= max_iterations; i++) {
    temp2 = matrix_multiply(a, temp1);
    free_matrix(temp1);
    temp1 = temp2;
    scalar_mult_matrix(t / ((double) i), temp1);
    sum_matrices(temp1, result);
  }
  free_matrix(temp2);
  for (i = 0; i < 0;  i++) {
    temp2 = matrix_multiply(result, result);
    free_matrix(result);
    result = temp2;
  }

  return result;

}
