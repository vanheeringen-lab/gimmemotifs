/* ---------------------------- Header ---------------------------------

  Module name   : utest.h
  Module prefix : T_

  Description: Provides funtions for unit tests.

  ---------------------------------------------------------------------

  Version:

    $Id$

  ---------------------------------------------------------------------

  Author   : S. Maetschke

  Copyright: Institute for Molecular Bioscience (IMB)

------------------------------------------------------------------------ */


#ifndef __UTEST
#define __UTEST




/* ---------------------------- Includes ------------------------------- */



/* ----------------------- Global Constants ---------------------------- */

/* ----------------------- Global Types -------------------------------- */
/*.......................................................................*/

/* ----------------------- Global Variables ---------------------------- */


/* ----------------------- Global Prototypes --------------------------- */

void T_set_verbose(int verbose_level);
/*
  Description   : Setter for the level of verbosity for the tests.
  Parameter     :
  verbose_level : 0: only failures are reported
                  1: tested file is reported
                  2: file and number of test performed are reported
                  3: the individual tests are reported.
  Example       :
    T_set_verbose(1);
*/

#define T_start() _T_start(__FILE__)
void _T_start(char *filename);
/*
  Description   : Starts the test.
  Parameter     :
    filename  : Name of the module/file to test.
  Example       :
    T_start(__FILE__);
    T_end();
*/

/*.......................................................................*/

void T_end();
/*
  Description   : Ends the test.
  Parameter     :
  Example       :
    T_start(__FILE__);
    T_end();
*/


/*.......................................................................*/


#define T_is_null(actual) _T_is_null(#actual, (actual), __LINE__, __FILE__)
void _T_is_null(char* expression, void* actual, int line, char *file);
/*
  Description   : Tests for NULL.
  Parameter     :
    acutal      : Actual value.
  Example       :
    T_is_null(p);
*/

/*.......................................................................*/

#define T_boolean(expected, actual) _T_boolean(#actual, (expected), (actual), __LINE__, __FILE__)
#define T_is_true(actual) _T_boolean(#actual, 1, (actual), __LINE__, __FILE__)
#define T_is_false(actual) _T_boolean(#actual, 0, (actual), __LINE__, __FILE__)
void _T_boolean(char* expression, int expected, int actual, int line, char *file);
/*
  Description   : Tests boolean values.
  Parameter     :
      expected  : Expected value.
      actual    : Actual value.
  Example       :
    T_boolean(1, 0);
*/

/*.......................................................................*/

#define T_string(expected, actual) _T_string(#actual, (expected), (actual), __LINE__, __FILE__)
void _T_string(char* expression, char* expected, char* actual, int line, char *file);
/*
  Description   : Tests string values.
  Parameter     :
      expected  : Expected value.
      actual    : Actual value.
  Example       :
    T_string("hello", "world");
*/

/*.......................................................................*/


#define T_int(expected, actual) _T_int(#actual, (expected), (actual), __LINE__, __FILE__)
void _T_int(char* expression, int expected, int actual, int line, char *file);
/*
  Description   : Tests integer values.
  Parameter     :
      expected  : Expected value.
      actual    : Actual value.
  Example       :
    T_int(11, 234);
*/

/*.......................................................................*/

#define T_char(expected, actual) _T_char(#actual, (expected), (actual), __LINE__, __FILE__)
void _T_char(char* expression, char expected, char actual, int line, char *file);
/*
  Description   : Tests char values.
  Parameter     :
      expected  : Expected value.
      actual    : Actual value.
  Example       :
    T_char('a', 'c');
*/

/*.......................................................................*/

#define T_double(expected, actual, eps) _T_double(#actual, (expected), (actual), eps, __LINE__, __FILE__)
void _T_double(char* expression, double expected, double actual, double eps, int line, char *file);
/*
  Description   : Tests boolean values.
  Parameter     :
      expected  : Expected value.
      actual    : Actual value.
      esp       : Permitted difference between the two values.
  Example       :
    T_double(1.00, 1.01, 0.01);
*/

/*.......................................................................*/




/*.......................................................................*/


/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/


/*.......................................................................*/


/*.......................................................................*/



#endif
