#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

/* The Perl include files perl.h redefines malloc and free. Here, we need the
 * usual malloc and free, defined in stdlib.h. So we undefine the ones in
 * perl.h.
 */

#ifdef malloc
#undef malloc
#endif
#ifdef free
#undef free
#endif

#include <stdlib.h>

#include "../src/cluster.h"


typedef struct {Node* nodes; int n;} Tree;

/* -------------------------------------------------
 * Using the warnings registry, check to see if warnings
 * are enabled for the Algorithm::Cluster module.
 */
static int
warnings_enabled(pTHX) {

    dSP;

    I32 count;
    bool isEnabled; 
    SV * mysv;

    ENTER ;
    SAVETMPS;
    PUSHMARK(SP) ;
    XPUSHs(sv_2mortal(newSVpv("Algorithm::Cluster",18)));
    PUTBACK ;

    count = perl_call_pv("warnings::enabled", G_SCALAR) ;

    if (count != 1) croak("No arguments returned from call_pv()\n") ;

    mysv = POPs; 
    isEnabled = (bool) SvTRUE(mysv); 

    PUTBACK ;
    FREETMPS ;
    LEAVE ;

    return isEnabled;
}

/* -------------------------------------------------
 * Create a row of doubles, initialized to a value
 */
static double*
malloc_row_dbl(pTHX_ int ncols, double val) {

    int j;
    double * row;

    row = malloc(ncols * sizeof(double) );
    if (!row) {
        return NULL;
    }

    for (j = 0; j < ncols; j++) { 
        row[j] = val;
    }
    return row;
}

/* -------------------------------------------------
 * Only coerce to a double if we already know it's 
 * an integer or double, or a string which is actually numeric.
 * Don't blindly run the macro SvNV, because that will coerce 
 * a non-numeric string to be a double of value 0.0, 
 * and we do not want that to happen, because if we test it again, 
 * it will then appear to be a valid double value. 
 */
static int
extract_double_from_scalar(pTHX_ SV * mysv, double * number) {

    if (SvPOKp(mysv) && SvLEN(mysv)) {  

        /* This function is not in the public perl API */
        if (Perl_looks_like_number(aTHX_ mysv)) {
            *number = SvNV( mysv );
            return 1;
        } else {
            return 0;
        } 
    } else if (SvNIOK(mysv)) {  
        *number = SvNV( mysv );
        return 1;
    } else {
        return 0;
    }
}



/* -------------------------------------------------
 * Convert a Perl 2D matrix into a 2D matrix of C doubles.
 * If no data are masked, mask can be passed as NULL.
 * NOTE: on errors this function returns a value greater than zero.
 */
static double**
parse_data(pTHX_ SV * matrix_ref, int** mask) {

    AV * matrix_av;
    SV * row_ref;
    AV * row_av;
    SV * cell;

    int type, i, j, nrows, ncols, n;

    double** matrix;

    /* NOTE -- we will just assume that matrix_ref points to an arrayref,
     * and that the first item in the array is itself an arrayref.
     * The calling perl functions must check this before we get this pointer.  
     * (It's easier to implement these checks in Perl rather than C.)
     * The value of perl_rows is now fixed. But the value of
     * rows will be decremented, if we skip any (invalid) Perl rows.
     */
    matrix_av  = (AV *) SvRV(matrix_ref);
    nrows = (int) av_len(matrix_av) + 1;

    if(nrows <= 0) {
        return NULL;
    }
    matrix   = malloc(nrows*sizeof(double*));
    if (!matrix) {
        return NULL;
    }

    row_ref  = *(av_fetch(matrix_av, (I32) 0, 0)); 
    row_av   = (AV *) SvRV(row_ref);
    ncols    = (int) av_len(row_av) + 1;


    /* ------------------------------------------------------------ 
     * Loop once for each row in the Perl matrix, and convert it to
     * C doubles. 
     */
    for (i=0; i < nrows; i++) { 

        row_ref = *(av_fetch(matrix_av, (I32) i, 0)); 

        if(! SvROK(row_ref) ) {

            if(warnings_enabled(aTHX))
                Perl_warn(aTHX_ 
                    "Row %d: Wanted array reference, but "
                    "got a scalar. No row to process?\n", i);
            break;
        }

        row_av = (AV *) SvRV(row_ref);
        type = SvTYPE(row_av); 
    
        /* Handle unexpected cases */
        if(type != SVt_PVAV ) {

             /* Reference doesn't point to an array at all. */
            if(warnings_enabled(aTHX))
                Perl_warn(aTHX_ 
                    "Row %d: Wanted array reference, but got "
                    "a reference to something else (%d)\n",
                    i, type);
            break;

        }

        n = (int) av_len(row_av) + 1;
        if (n != ncols) {
            /* All rows in the matrix should have the same
             * number of columns. */
            if(warnings_enabled(aTHX))
                Perl_warn(aTHX_ 
                    "Row %d: Contains %d columns "
                    "(expected %d)\n", i, n, ncols);
            break;
        }

        matrix[i] = malloc(ncols*sizeof(double));
        if (!matrix[i])
            break;

        /* Loop once for each cell in the row. */
        for (j=0; j < ncols; j++) { 
        
            double num;
            if (!mask || mask[i][j]) {
                cell = *(av_fetch(row_av, (I32) j, 0)); 
                if(extract_double_from_scalar(aTHX_ cell,&num) <= 0) {    
                    if(warnings_enabled(aTHX))
                        Perl_warn(aTHX_ 
                            "Row %d col %d: Value is not "
                                                    "a number.\n", i, j);
                    free(matrix[i]); /* not included below */
                    break;
                }
            }
            else {
                /* Don't read the value if it is masked.
                 * Set it to some arbitrary value. */
                num = 0.0;
            }
            matrix[i][j] = num;

        } /* End for (j=0; j < ncols; j++) */
        if (j < ncols) break;

    } /* End for (i=0; i < nrows; i++) */

    if (i < nrows) { /* encountered a break */
        nrows = i;
        for (i = 0; i < nrows; i++) free(matrix[i]);
        free(matrix);
        matrix = NULL;
    }

    return matrix;
}


/* -------------------------------------------------
 * Convert a Perl 2D matrix into a 2D matrix of C ints.
 * On errors this function returns a value greater than zero.
 */
static int**
parse_mask(pTHX_ SV * matrix_ref) {

    AV * matrix_av;
    SV * row_ref;
    AV * row_av;
    SV * cell;

    int type, i, j, nrows, ncols, n;

    int** matrix;

    /* NOTE -- we will just assume that matrix_ref points to an arrayref,
     * and that the first item in the array is itself an arrayref.
     * The calling perl functions must check this before we get this pointer.  
     * (It's easier to implement these checks in Perl rather than C.)
     * The value of perl_rows is now fixed. But the value of
     * rows will be decremented, if we skip any (invalid) Perl rows.
     */
    matrix_av = (AV *) SvRV(matrix_ref);
    nrows = (int) av_len(matrix_av) + 1;

    if(nrows <= 0) {
        return NULL;  /* Caller must handle this case!! */
    }
    matrix    = malloc(nrows * sizeof(int *) );
    if (!matrix) {
        return NULL;
    }

    row_ref   = *(av_fetch(matrix_av, (I32) 0, 0)); 
    row_av    = (AV *) SvRV(row_ref);
    ncols     = (int) av_len(row_av) + 1;



    /* ------------------------------------------------------------ 
     * Loop once for each row in the Perl matrix, and convert it to C ints. 
     */
    for (i=0; i < nrows; i++) { 

        row_ref = *(av_fetch(matrix_av, (I32) i, 0)); 

        if(! SvROK(row_ref) ) {

            if(warnings_enabled(aTHX))
                Perl_warn(aTHX_ 
                    "Row %d: Wanted array reference, but "
                    "got a scalar. No row to process?\n", i);
            break;
        }

        row_av = (AV *) SvRV(row_ref);
        type = SvTYPE(row_av); 
    
        /* Handle unexpected cases */
        if(type != SVt_PVAV ) {

             /* Reference doesn't point to an array at all. */
            if(warnings_enabled(aTHX))
                Perl_warn(aTHX_ 
                    "Row %d: Wanted array reference, but got "
                    "a reference to something else (%d)\n",
                    i, type);
            break;

        }

        n = (int) av_len(row_av) + 1;
        if (n != ncols) {
            /* All rows in the matrix should have the same
             * number of columns. */
            if(warnings_enabled(aTHX))
                Perl_warn(aTHX_ 
                    "Row %d: Contains %d columns "
                    "(expected %d)\n", i, n, ncols);
            break;
        }

        matrix[i] = malloc(ncols * sizeof(int) );
        if (!matrix[i]) {
            break;
        }

        /* Loop once for each cell in the row. */
        for (j=0; j < ncols; ++j) { 
            double num;
            cell = *(av_fetch(row_av, (I32) j, 0)); 
            if(extract_double_from_scalar(aTHX_ cell,&num) <= 0) {    
                if(warnings_enabled(aTHX))
                    Perl_warn(aTHX_
                        "Row %d col %d: Value is not "
                        "a number.\n", i, j);
                free(matrix[i]); /* not included below */
                break;
            }
            matrix[i][j] = (int) num;

        } /* End for (j=0; j < ncols; j++) */
        if (j < ncols) break;

    } /* End for (i=0; i < nrows; i++) */

    if (i < nrows) { /* break statement encountered */
        nrows = i;
        for (i = 0; i < nrows; i++) free(matrix[i]);
        free(matrix);
        matrix = NULL;
    }

    return matrix;
}


/* -------------------------------------------------
 *
 */
static void
free_matrix_int(int ** matrix, int nrows) {

    int i;
    for(i = 0; i < nrows; ++i ) {
        free(matrix[i]);
    }

    free(matrix);
}


/* -------------------------------------------------
 *
 */
static void
free_matrix_dbl(double ** matrix, int nrows) {

    int i;
    for(i = 0; i < nrows; ++i ) {
        free(matrix[i]);
    }

    free(matrix);
}


/* -------------------------------------------------
 *
 */
static void
free_ragged_matrix_dbl(double ** matrix, int nrows) {

    int i;
    for(i = 1; i < nrows; ++i ) {
        free(matrix[i]);
    }

    free(matrix);
}


/* -------------------------------------------------
 * Convert a Perl array into an array of doubles
 * On error, this function returns NULL.
 */
static double*
malloc_row_perl2c_dbl (pTHX_ SV * input, int* np) {

    int i;
    AV* array    = (AV *) SvRV(input);
    const int n  = (int) av_len(array) + 1;
    double* data = malloc(n * sizeof(double)); 
    if (!data) {
        return NULL;
    }

    /* Loop once for each item in the Perl array, and convert
         * it to a C double. 
     */
    for (i=0; i < n; i++) {
        double num;
        SV * mysv = *(av_fetch(array, (I32) i, (I32) 0));
        if(extract_double_from_scalar(aTHX_ mysv,&num) > 0) {    
            data[i] = num;
        } else {
            /* Error reading data */
            if (warnings_enabled(aTHX))
                Perl_warn(aTHX_ 
                    "Error parsing array: item %d is not a number\n", i);      
            free(data);
            return NULL;
        }
    }
    if(np) *np = n;
    return data;
}

/* -------------------------------------------------
 * Convert a Perl array into an array of ints
 * On errors this function returns NULL.
 */
static int*
malloc_row_perl2c_int (pTHX_ SV * input) {

    int i;

    AV* array = (AV *) SvRV(input);
    const int n = (int) av_len(array) + 1;
    int* data = malloc(n*sizeof(int)); 
    if (!data) {
        return NULL;
    }

    /* Loop once for each item in the Perl array,
     * and convert it to a C double. 
     */
    for (i=0; i < n; i++) {
        double num;
        SV * mysv = *(av_fetch(array, (I32) i, (I32) 0));
        if(extract_double_from_scalar(aTHX_ mysv,&num) > 0) {    
            data[i] = (int) num;
        } else {
            /* Check if the item is numeric */
            if (warnings_enabled(aTHX))
                Perl_warn(aTHX_ "Error when parsing array: item %d is"
                    " not a number, skipping\n", i);      
            free(data);
            return NULL;
        }
    }

    return data;
}

/* -------------------------------------------------
 * Copy a Perl array into an array of ints.
 * If an error occurs, return 0; otherwise return 1.
 */
static int
copy_row_perl2c_int (pTHX_ SV * input, int* output) {

    int i;

    AV* array = (AV *) SvRV(input);
    const int n = (int) av_len(array) + 1;

    /* Loop once for each item in the Perl array,
     * and convert it to a C double. 
     */
    for (i=0; i < n; i++) {
        double num;
        SV * mysv = *(av_fetch(array, (I32) i, (I32) 0));
        if(extract_double_from_scalar(aTHX_ mysv,&num) > 0) {    
            output[i] = (int) num;
        } else {
            /* Skip any items which are not numeric */
            if (warnings_enabled(aTHX))
                Perl_warn(aTHX_ 
                    "Error when parsing array: item %d is"
                    " not a number\n", i);      
            return 0;
        }
    }
    return 1;
}
/* -------------------------------------------------
 *
 */
static SV *
row_c2perl_dbl(pTHX_ double * row, int ncols) {

    int j;
    AV * row_av = newAV();
    for(j=0; j<ncols; ++j) {
        av_push(row_av, newSVnv(row[j]));
    }
    return newRV_noinc((SV*) row_av);
}

/* -------------------------------------------------
 *
 */
static SV*
row_c2perl_int(pTHX_ int * row, int ncols) {

    int j;
    AV * row_av = newAV();
    for(j=0; j<ncols; ++j) {
        av_push(row_av, newSVnv(row[j]));
    }
    return ( newRV_noinc( (SV*) row_av ) );
}

/* -------------------------------------------------
 *
 */
static SV*
matrix_c2perl_dbl(pTHX_ double ** matrix, int nrows, int ncols) {

    int i;
    AV * matrix_av = newAV();
    SV * row_ref;
    for(i=0; i<nrows; ++i) {
        row_ref = row_c2perl_dbl(aTHX_ matrix[i], ncols);
        av_push(matrix_av, row_ref);
    }
    return ( newRV_noinc( (SV*) matrix_av ) );
}

/* -------------------------------------------------
 *
 */
static SV*
matrix_c2perl_int(pTHX_ int ** matrix, int nrows, int ncols) {

    int i;
    AV * matrix_av = newAV();
    SV * row_ref;
    for(i=0; i<nrows; ++i) {
        row_ref = row_c2perl_int(aTHX_ matrix[i], ncols);
        av_push(matrix_av, row_ref);
    }
    return ( newRV_noinc( (SV*) matrix_av ) );
}

/* -------------------------------------------------
 *
 */
static SV*
ragged_matrix_c2perl_dbl(pTHX_ double ** matrix, int nobjects) {

    int i;
    AV * matrix_av = newAV();
    SV * row_ref;
    for(i=0; i<nobjects; ++i) {
        row_ref = row_c2perl_dbl(aTHX_ matrix[i], i);
        av_push(matrix_av, row_ref);
    }
    return ( newRV_noinc( (SV*) matrix_av ) );
}

/* -------------------------------------------------
 * Check if the data matrix is a distance matrix, or
 * a raw distance matrix.
 */
static int
is_distance_matrix(pTHX_ SV * data_ref)
{
    /* We don't check data_ref because we expect the caller to check it 
     */
    AV * matrix_av  = (AV *) SvRV(data_ref);
    SV * row_ref    = *(av_fetch(matrix_av, (I32) 0, 0)); 
    AV * row_av     = (AV *) SvRV(row_ref);
    const int ncols = (int) av_len(row_av) + 1;
    if (ncols==0) return 1;

    return 0;
}


/* -------------------------------------------------
 * Convert the 'data' and 'mask' matrices and the 'weight' array
 * from C to Perl.  Also check for errors, and bail out if there are any.
 */
static int
malloc_matrices(pTHX_
    SV *  weight_ref, double  ** weight, int nweights, 
    SV *  data_ref,   double *** matrix,
    SV *  mask_ref,   int    *** mask,
    int   nrows,      int        ncols
) {

    if(SvROK(mask_ref) && SvTYPE(SvRV(mask_ref)) == SVt_PVAV) { 
        *mask = parse_mask(aTHX_ mask_ref);
        if(*mask==NULL) return 0;
    } else {
        int i,j;
        int** p = malloc(nrows*sizeof(int*));
        if(!p) return 0;
        for (i = 0; i < nrows; ++i) { 
            p[i] = malloc(ncols*sizeof(int));
            if(!p[i]) {
                while(--i >= 0) free(p[i]);
                free(p);
                return 0;
            }
            for (j = 0; j < ncols; j++) p[i][j] = 1;
        }
        *mask = p;
    }

    /* We don't check data_ref because we expect the caller to check it 
     */
    *matrix = parse_data(aTHX_ data_ref, *mask);
    if(*matrix==NULL) {
        free_matrix_int(*mask,     nrows);
        return 0;
    }

    if(weight_ref==NULL) return 1; /* Weights not needed */
    if(SvROK(weight_ref) && SvTYPE(SvRV(weight_ref)) == SVt_PVAV) { 
        *weight = malloc_row_perl2c_dbl(aTHX_ weight_ref, NULL);
    } else {
        *weight = malloc_row_dbl(aTHX_ nweights,1.0);
    }

    if(!(*weight)) {
        free_matrix_int(*mask,     nrows);
        free_matrix_dbl(*matrix,   nrows);
        return 0;
    }

    return 1;
}

static double**
parse_distance(pTHX_ SV* matrix_ref, int nobjects)
{
    int i,j;

    AV* matrix_av  = (AV *) SvRV(matrix_ref);
    double** matrix = malloc(nobjects*sizeof(double*));
    if (!matrix) {
        return NULL;
    }

    matrix[0] = NULL;
    for (i=1; i < nobjects; i++) { 
        SV* row_ref = *(av_fetch(matrix_av, (I32) i, 0)); 
        AV* row_av  = (AV *) SvRV(row_ref);
        matrix[i] = malloc(i * sizeof(double));
        if (!matrix[i]) {
            break;
        }
        /* Loop once for each cell in the row. */
        for (j=0; j < i; j++) { 
            double num;
            SV* cell = *(av_fetch(row_av, (I32) j, 0)); 
            if(extract_double_from_scalar(aTHX_ cell,&num) > 0) {    
                matrix[i][j] = num;
            } else {
                if(warnings_enabled(aTHX))
                    Perl_warn(aTHX_ 
                        "Row %d col %d: Value is not "
                                                "a number.\n", i, j);
                break;
            }
        }
    }

    if (i < nobjects) {
        nobjects = i+1;
        for (i = 1; i < nobjects; i++) free(matrix[i]);
        free(matrix);
        matrix = NULL;
    }

    return matrix;
}

/******************************************************************************/
/**                                                                          **/
/** XS code begins here                                                      **/
/**                                                                          **/
/******************************************************************************/
/******************************************************************************/

MODULE = Algorithm::Cluster PACKAGE = Algorithm::Cluster::Node
PROTOTYPES: ENABLE

SV*
new (class, left, right, distance)
    char* class
    int left
    int right
    double distance
    PREINIT:
    Node* node;
    SV* obj;
    CODE:
    node = malloc(sizeof(Node));
    RETVAL = newSViv(0);
    obj = newSVrv(RETVAL, class);
    node->left = left;
    node->right = right;
    node->distance = distance;

    sv_setiv(obj, PTR2IV(node));
    SvREADONLY_on(obj);
    OUTPUT:
    RETVAL


int
left (obj)
    SV* obj
    CODE:
    RETVAL = (INT2PTR(Node*,SvIV(SvRV(obj))))->left;
    OUTPUT:
    RETVAL

int
right (obj)
    SV* obj
    CODE:
    RETVAL = (INT2PTR(Node*,SvIV(SvRV(obj))))->right;
    OUTPUT:
    RETVAL

double
distance (obj)
    SV* obj
    CODE:
    RETVAL = (INT2PTR(Node*,SvIV(SvRV(obj))))->distance;
    OUTPUT:
    RETVAL

void
set_left (obj, left)
    SV* obj
    int left
    PREINIT:
    Node* node;
    CODE:
    if (!sv_isa(obj, "Algorithm::Cluster::Node")) {
        croak("set_left should be applied to an Algorithm::Cluster::Node object");
    }
    node = INT2PTR(Node*,SvIV(SvRV(obj)));
    node->left = left;

void
set_right (obj, right)
    SV* obj
    int right
    PREINIT:
    Node* node;
    CODE:
    if (!sv_isa(obj, "Algorithm::Cluster::Node")) {
        croak("set_right should be applied to an Algorithm::Cluster::Node object");
    }
    node = INT2PTR(Node*,SvIV(SvRV(obj)));
    node->right = right;

void
set_distance (obj, distance)
    SV* obj
    double distance
    PREINIT:
    Node* node;
    CODE:
    if (!sv_isa(obj, "Algorithm::Cluster::Node")) {
        croak("set_distance should be applied to an Algorithm::Cluster::Node object");
    }
    node = INT2PTR(Node*,SvIV(SvRV(obj)));
    node->distance = distance;

void DESTROY (obj)
    SV* obj
    PREINIT:
    I32* temp;
    Node* node;
    PPCODE:
    temp = PL_markstack_ptr++;
    node = INT2PTR(Node*, SvIV(SvRV(obj)));
    free(node);
    if (PL_markstack_ptr != temp) {
        /* truly void, because dXSARGS not invoked */
        PL_markstack_ptr = temp;
        XSRETURN_EMPTY;
        /* return empty stack */
    }  /* must have used dXSARGS; list context implied */
    return;  /* assume stack size is correct */


MODULE = Algorithm::Cluster PACKAGE = Algorithm::Cluster::Tree
PROTOTYPES: ENABLE

SV*
new (class, nodes)
    char* class
    SV* nodes

    PREINIT:
    Tree* tree;
    SV* obj;
        int i;
    int n;
    AV* array;
    int* flag;

    CODE:
    if(!SvROK(nodes) || SvTYPE(SvRV(nodes)) != SVt_PVAV) { 
        croak("Algorithm::Cluster::Tree::new expects an array of nodes\n");
    }
    array = (AV *) SvRV(nodes);
    n = (int) av_len(array) + 1;
    tree = malloc(sizeof(Tree));
    if (tree) {
        tree->n = n;
        tree->nodes = malloc(n*sizeof(Node));
    }
    if (! tree || !tree->nodes) {
        if (tree) free(tree);
        croak("Algorithm::Cluster::Tree::new memory error\n");
    }

        for (i = 0; i < n; i++) {
        Node* node;
        SV* node_ref = *(av_fetch(array, (I32) i, 0)); 
        if (!sv_isa(node_ref, "Algorithm::Cluster::Node")) break;
        node = INT2PTR(Node*,SvIV(SvRV(node_ref)));
        tree->nodes[i].left = node->left;
        tree->nodes[i].right = node->right;
        tree->nodes[i].distance = node->distance;
    }

    if (i < n) {
        /* break encountered */
        free(tree->nodes);
        free(tree);
        croak("Algorithm::Cluster::Tree::new expects an array of nodes\n");
    }

    flag = malloc((2*n+1)*sizeof(int));
    if(flag) {
         int j;
        for (i = 0; i < 2*n+1; i++) flag[i] = 0;
        for (i = 0; i < n; i++) {
            j = tree->nodes[i].left;
            if (j < 0) {
                j = -j-1;
                if (j>=i) break;
            }
            else j+=n;
            if (flag[j]) break;
            flag[j] = 1;
            j = tree->nodes[i].right;
            if (j < 0) {
                j = -j-1;
                if (j>=i) break;
            }
            else j+=n;
            if (flag[j]) break;
            flag[j] = 1;
        }
        free(flag);
    }

    if (!flag || i < n) {
        /* break encountered */
        free(tree->nodes);
        free(tree);
        croak("the array of nodes passed to Algorithm::Cluster::Tree::new does not represent a valid tree\n");
    }

    RETVAL = newSViv(0);
    obj = newSVrv(RETVAL, class);
    sv_setiv(obj, PTR2IV(tree));
    SvREADONLY_on(obj);

    OUTPUT:
    RETVAL

int
length (obj)
    SV* obj
    CODE:
    RETVAL = (INT2PTR(Tree*,SvIV(SvRV(obj))))->n;
    OUTPUT:
    RETVAL


SV *
get (obj, index)
    SV* obj
    int index
    PREINIT:
    Tree* tree;
    Node* node;
    SV* scalar;
    CODE:
    tree = INT2PTR(Tree*,SvIV(SvRV(obj)));
    if (index < 0 || index >= tree->n) {
        croak("Index out of bounds in Algorithm::Cluster::Tree::get\n");
    }
    RETVAL = newSViv(0);
    scalar = newSVrv(RETVAL, "Algorithm::Cluster::Node");
    node = malloc(sizeof(Node));
    node->left = tree->nodes[index].left;
    node->right = tree->nodes[index].right;
    node->distance = tree->nodes[index].distance;
    sv_setiv(scalar, PTR2IV(node));
    SvREADONLY_on(scalar);
    OUTPUT:
    RETVAL

void
scale(obj)
    SV* obj
    PREINIT:
    int i;
    int n;
    Tree* tree;
    Node* nodes;
    double maximum;
    CODE:
    if (!sv_isa(obj, "Algorithm::Cluster::Tree")) {
        croak("scale can only be applied to an Algorithm::Cluster::Tree object");
    }
    tree = INT2PTR(Tree*,SvIV(SvRV(obj)));
    n = tree->n;
    nodes = tree->nodes;
    maximum = DBL_MIN;
    for (i = 0; i < n; i++) {
        double distance = nodes[i].distance;
        if (distance > maximum) maximum = distance;
    }
    if (maximum!=0.0) {
        for (i = 0; i < n; i++) nodes[i].distance /= maximum;
    }

AV *
cut(obj, nclusters)
    SV* obj
    int nclusters
    PREINIT:
    int i;
    int n;
    Tree* tree;
    int* clusterid;
    CODE:
    if (!sv_isa(obj, "Algorithm::Cluster::Tree")) {
        croak("cut can only be applied to an Algorithm::Cluster::Tree object");
    }
    tree = INT2PTR(Tree*,SvIV(SvRV(obj)));
    n = tree->n + 1;
    if (nclusters < 1) {
        croak("cut: Requested number of clusters should be positive");
    }
    if (nclusters > n) {
        croak("cut: More clusters requested than items available");
    }
    clusterid = malloc(n*sizeof(int));
    if (!clusterid) {
        croak("cut: Insufficient memory");
    }
        /* --------------------------------------------------------------- */
    cuttree(n, tree->nodes, nclusters, clusterid);
    /* -- Check for errors flagged by the C routine ------------------ */
    if (clusterid[0]==-1) {
        free(clusterid);
        croak("cut: Error in the cuttree routine");
    }
    RETVAL = newAV();
    for(i=0; i<n; i++) {
        av_push(RETVAL, newSVnv(clusterid[i]));
    }
    free(clusterid);
    sv_2mortal((SV*)RETVAL);
    OUTPUT:
    RETVAL


void DESTROY (obj)
    SV* obj
    PREINIT:
    I32* temp;
    Tree* tree;
    PPCODE:
    temp = PL_markstack_ptr++;
    tree = INT2PTR(Tree*, SvIV(SvRV(obj)));
    free(tree->nodes);
    free(tree);
    if (PL_markstack_ptr != temp) {
        /* truly void, because dXSARGS not invoked */
        PL_markstack_ptr = temp;
        XSRETURN_EMPTY;
        /* return empty stack */
    }  /* must have used dXSARGS; list context implied */
    return;  /* assume stack size is correct */


MODULE = Algorithm::Cluster    PACKAGE = Algorithm::Cluster
PROTOTYPES: ENABLE


SV *
_version()
    CODE:
    RETVAL = newSVpv( CLUSTERVERSION , 0);

    OUTPUT:
    RETVAL



SV *
_mean(input)
    SV * input;

    PREINIT:
    int array_length;
    double * data;  /* one-dimensional array of doubles */

    CODE:
    if(SvTYPE(SvRV(input)) != SVt_PVAV) { 
        XSRETURN_UNDEF;
    }

    data = malloc_row_perl2c_dbl (aTHX_ input, &array_length);
    if (data) {
        RETVAL = newSVnv( mean(array_length, data) );
        free(data);
    } else {
        croak("memory allocation failure in _mean\n");
    }

    OUTPUT:
    RETVAL


SV *
_median(input)
    SV * input;

    PREINIT:
    int array_length;
    double * data;  /* one-dimensional array of doubles */

    CODE:
    if(SvTYPE(SvRV(input)) != SVt_PVAV) { 
        XSRETURN_UNDEF;
    }

    data = malloc_row_perl2c_dbl (aTHX_ input, &array_length);
    if (data) {
        RETVAL = newSVnv( median(array_length, data) );
        free(data);
    } else {
        croak("memory allocation failure in _median\n");
    }

    OUTPUT:
    RETVAL


SV *
_treecluster(nrows,ncols,data_ref,mask_ref,weight_ref,transpose,dist,method)
    int      nrows;
    int      ncols;
    SV *     data_ref;
    SV *     mask_ref;
    SV *     weight_ref;
    int      transpose;
    char *   dist;
    char *   method;

    PREINIT:
    Node*    nodes;

    double  * weight = NULL;
    double ** matrix = NULL;
    int    ** mask   = NULL;
    double ** distancematrix = NULL;
    const int ndata = transpose ? nrows : ncols;
    const int nelements = transpose ? ncols : nrows;

    CODE:
    /* ------------------------
     * Don't check the parameters, because we rely on the Perl
     * caller to check most paramters.
     */
    /* ------------------------
     * Convert data and mask matrices and the weight array
     * from C to Perl.  Also check for errors, and ignore the
     * mask or the weight array if there are any errors. 
     */
    if (is_distance_matrix(aTHX_ data_ref)) {
        distancematrix = parse_distance(aTHX_ data_ref, nelements);
        if (!distancematrix) {
                croak("memory allocation failure in _treecluster\n");
        }
    } else {
        int ok;
        ok = malloc_matrices(aTHX_ weight_ref, &weight, ndata, 
                    data_ref,   &matrix,
                    mask_ref,   &mask,  
                    nrows,      ncols);
        if (!ok) {
            croak("failed to read input data for _treecluster\n");
        }
        }

    /* ------------------------
     * Run the library function
     */
    nodes = treecluster(nrows, ncols, matrix, mask, weight, transpose,
                dist[0], method[0], distancematrix);

    /* ------------------------
     * Check result to make sure we didn't run into memory problems
     */
    if(!nodes) {
        /* treecluster failed due to insufficient memory */
        if (matrix) {
            free_matrix_int(mask,     nrows);
            free_matrix_dbl(matrix,   nrows);
            free(weight);
        } else {
            free_ragged_matrix_dbl(distancematrix, nelements);
        }
        croak("memory allocation failure in treecluster\n");
    }
    else {

        /* ------------------------
          * Convert generated C matrices to Perl matrices
          */
        const int n = nelements-1;
        int i;
        SV* obj;
        Tree* tree;
        RETVAL = newSViv(0);
        obj = newSVrv(RETVAL, "Algorithm::Cluster::Tree");
        tree = malloc(sizeof(Tree));
        tree->n = n;
        tree->nodes = malloc(n*sizeof(Node));
        sv_setiv(obj, PTR2IV(tree));
        SvREADONLY_on(obj);
        for(i=0; i<n; i++) {
            tree->nodes[i].left = nodes[i].left;
            tree->nodes[i].right = nodes[i].right;
            tree->nodes[i].distance = nodes[i].distance;
        }
        free(nodes);
    }

    /* ------------------------
     * Free what we've malloc'ed 
     */
    if (matrix) {
        free_matrix_int(mask,     nrows);
        free_matrix_dbl(matrix,   nrows);
        free(weight);
    } else {
        free_ragged_matrix_dbl(distancematrix, nelements);
    }

    /* Finished _treecluster() */
    OUTPUT:
    RETVAL


void
_kcluster(nclusters,nrows,ncols,data_ref,mask_ref,weight_ref,transpose,npass,method,dist,initialid_ref)
    int      nclusters;
    int      nrows;
    int      ncols;
    SV *     data_ref;
    SV *     mask_ref;
    SV *     weight_ref;
    int      transpose;
    int      npass;
    char *   method;
    char *   dist;
    SV *     initialid_ref;

    PREINIT:
    SV  *    clusterid_ref;
    int *    clusterid;
    int      nobjects;
    int      ndata;
    double   error;
    int      ifound;
    int      ok;

    double  * weight;
    double ** matrix;
    int    ** mask;


    PPCODE:
    /* ------------------------
     * Don't check the parameters, because we rely on the Perl
     * caller to check most parameters.
     */

    /* ------------------------
     * Malloc space for the return values from the library function
     */
    if (transpose==0) {
        nobjects = nrows;
        ndata = ncols;
    } else {
        nobjects = ncols;
        ndata = nrows;
    }
        clusterid = malloc(nobjects * sizeof(int) );
        if (!clusterid) {
        croak("memory allocation failure in _kcluster\n");
    }

    /* ------------------------
     * Convert data and mask matrices and the weight array
     * from C to Perl.  Also check for errors, and ignore the
     * mask or the weight array if there are any errors. 
     */
    ok = malloc_matrices( aTHX_ weight_ref, &weight, ndata, 
                data_ref,   &matrix,
                mask_ref,   &mask,  
                nrows,      ncols);
        if (!ok) {
        free(clusterid);
            croak("failed to read input data for _kcluster\n");
    }

    /* ------------------------
     * Copy initialid to clusterid, if needed
     */

    if (npass==0) {
        copy_row_perl2c_int(aTHX_ initialid_ref, clusterid);
    }

    /* ------------------------
     * Run the library function
     */
    kcluster( 
        nclusters, nrows, ncols, 
        matrix, mask, weight, transpose,
        npass, method[0], dist[0], clusterid,  &error, &ifound
        
    );

    /* ------------------------
     * Convert generated C matrices to Perl matrices
     */
    clusterid_ref =    row_c2perl_int(aTHX_ clusterid, nobjects);

    /* ------------------------
     * Push the new Perl matrices onto the return stack
     */
    XPUSHs(sv_2mortal( clusterid_ref   ));
    XPUSHs(sv_2mortal( newSVnv(error) ));
    XPUSHs(sv_2mortal( newSViv(ifound) ));

    /* ------------------------
     * Free what we've malloc'ed 
     */
    free(clusterid);
    free_matrix_int(mask,     nrows);
    free_matrix_dbl(matrix,   nrows);
    free(weight);

    /* Finished _kcluster() */



void
_kmedoids(nclusters,nobjects,distancematrix_ref,npass,initialid_ref)
    int      nclusters;
    int      nobjects;
    SV *     distancematrix_ref;
    int      npass;
    SV *     initialid_ref;


    PREINIT:
    double** distancematrix;
    SV  *    clusterid_ref;
    int *    clusterid;
    double   error;
    int      ifound;



    PPCODE:
    /* ------------------------
     * Don't check the parameters, because we rely on the Perl
     * caller to check most parameters.
     */

    /* ------------------------
     * Malloc space for the return values from the library function
     */
    clusterid = malloc(nobjects * sizeof(int) );
        if (!clusterid) {
            croak("memory allocation failure in _kmedoids\n");
    }

    /* ------------------------
     * Convert data and mask matrices and the weight array
     * from C to Perl.  Also check for errors, and ignore the
     * mask or the weight array if there are any errors. 
     */
    distancematrix = parse_distance(aTHX_ distancematrix_ref, nobjects);
    if (!distancematrix) {
        free(clusterid);
            croak("failed to allocate memory for distance matrix in _kmedoids\n");
    }

    /* ------------------------
     * Copy initialid to clusterid, if needed
     */

    if (npass==0) {
        copy_row_perl2c_int(aTHX_ initialid_ref, clusterid);
    }

    /* ------------------------
     * Run the library function
     */
    kmedoids( 
        nclusters, nobjects, 
        distancematrix, npass, clusterid, 
        &error, &ifound
    );

    if(ifound==-1) {
        free(clusterid);
        free_ragged_matrix_dbl(distancematrix, nobjects);
            croak("memory allocation failure in _kmedoids\n");
    }
    else if(ifound==0) {
        free(clusterid);
        free_ragged_matrix_dbl(distancematrix, nobjects);
            croak("error in input arguments in kmedoids\n");
    }
    else {

        /* ------------------------
         * Convert generated C matrices to Perl matrices
         */
        clusterid_ref =    row_c2perl_int(aTHX_ clusterid, nobjects);

        /* ------------------------
         * Push the new Perl matrices onto the return stack
         */
        XPUSHs(sv_2mortal( clusterid_ref   ));
        XPUSHs(sv_2mortal( newSVnv(error) ));
        XPUSHs(sv_2mortal( newSViv(ifound) ));

    }
    /* ------------------------
     * Free what we've malloc'ed 
     */
    free(clusterid);
    free_ragged_matrix_dbl(distancematrix, nobjects);

    /* Finished _kmedoids() */



double
_clusterdistance(nrows,ncols,data_ref,mask_ref,weight_ref,cluster1_len,cluster2_len,cluster1_ref,cluster2_ref,dist,method,transpose)
    int      nrows;
    int      ncols;
    SV *     data_ref;
    SV *     mask_ref;
    SV *     weight_ref;
    int      cluster1_len;
    int      cluster2_len;
    SV *     cluster1_ref;
    SV *     cluster2_ref;
    char *   dist;
    char *   method;
    int      transpose;

    PREINIT:
    int   nweights;

    int     * cluster1;
    int     * cluster2;

    double  * weight;
    double ** matrix;
    int    ** mask;

    double distance;

    int ok;

    CODE:

    /* ------------------------
     * Don't check the parameters, because we rely on the Perl
     * caller to check most paramters.
     */

    /* ------------------------
     * Convert cluster index Perl arrays to C arrays
     */
    cluster1 = malloc_row_perl2c_int(aTHX_ cluster1_ref);
    cluster2 = malloc_row_perl2c_int(aTHX_ cluster2_ref);
    if (!cluster1 || !cluster2) {
        if (cluster1) free(cluster1);
        if (cluster2) free(cluster2);
        croak("memory allocation failure in _clusterdistance\n");
    }

    /* ------------------------
     * Convert data and mask matrices and the weight array
     * from C to Perl.  Also check for errors, and ignore the
     * mask or the weight array if there are any errors. 
     * Set nweights to the correct number of weights.
     */
    nweights = (transpose==0) ? ncols : nrows;
    ok = malloc_matrices( aTHX_ weight_ref, &weight, nweights, 
                data_ref,   &matrix,
                mask_ref,   &mask,  
                nrows,      ncols);
    if (!ok) {
        free(cluster1);
        free(cluster2);
            croak("failed to read input data for _clusterdistance\n");
    }

    /* ------------------------
     * Run the library function
     */
    distance = clusterdistance( 
        nrows, ncols, 
        matrix, mask, weight,
        cluster1_len, cluster2_len, cluster1, cluster2,
        dist[0], method[0], transpose
    );

    RETVAL = distance;

    /* ------------------------
     * Free what we've malloc'ed 
     */
    free_matrix_int(mask,     nrows);
    free_matrix_dbl(matrix,   nrows);
    free(weight);
    free(cluster1);
    free(cluster2);

    /* Finished _clusterdistance() */

    OUTPUT:
    RETVAL



void
_clustercentroids(nclusters,nrows,ncols,data_ref,mask_ref,clusterid_ref,transpose,method)
    int      nclusters;
    int      nrows;
    int      ncols;
    SV *     data_ref;
    SV *     mask_ref;
    SV *     clusterid_ref;
    int      transpose;
    char *   method;

    PREINIT:
    SV  *    cdata_ref;
    SV  *    cmask_ref;
    int     * clusterid;

    double ** matrix;
    int    ** mask;

    double ** cdata;
    int    ** cmask;
    int       cnrows = 0; /* Initialize to make the compiler shut up */
    int       cncols = 0; /* Initialize to make the compiler shut up */

    int i;
    int ok;

    PPCODE:

    /* ------------------------
     * Don't check the parameters, because we rely on the Perl
     * caller to check most paramters.
     */
        if (transpose==0)
        {    cnrows = nclusters;
             cncols = ncols;
        }
        else if (transpose==1)
        {    cnrows = nrows;
             cncols = nclusters;
        }

    /* ------------------------
     * Convert cluster index Perl arrays to C arrays
     */
    clusterid = malloc_row_perl2c_int(aTHX_ clusterid_ref);
    if (!clusterid) {
        croak("memory allocation failure in _clustercentroids\n");
    }

    /* ------------------------
     * Convert data and mask matrices and the weight array
     * from C to Perl.  Also check for errors, and ignore the
     * mask or the weight array if there are any errors. 
     * Set nweights to the correct number of weights.
     */
    ok = malloc_matrices( aTHX_ NULL, NULL, 0, 
                data_ref,   &matrix,
                mask_ref,   &mask,  
                nrows,      ncols);
    if (!ok) {
        free(clusterid);
            croak("failed to read input data for _clustercentroids\n");
    }


    /* ------------------------
     * Create the output variables cdata and cmask.
     */
    i = 0;
    cdata = malloc(cnrows * sizeof(double*));
    cmask = malloc(cnrows * sizeof(int*));
    if (cdata && cmask) {
        for ( ; i < cnrows; i++) {
            cdata[i] = malloc(cncols*sizeof(double));
            cmask[i] = malloc(cncols*sizeof(int));
            if (!cdata[i] || !cmask[i]) break;
        }
    }
    if (i < cnrows)
    {
        if (cdata[i]) free(cdata[i]);
        if (cmask[i]) free(cmask[i]);
        while (--i >= 0) {
            free(cdata[i]);
            free(cmask[i]);
        }
        if (cdata) free(cdata);
        if (cmask) free(cmask);
        free(clusterid);
        free_matrix_int(mask,     nrows);
        free_matrix_dbl(matrix,   nrows);
        croak("memory allocation failure in _clustercentroids\n");
    }

    /* ------------------------
     * Run the library function
     */
    ok = getclustercentroids(
               nclusters, nrows, ncols,
               matrix, mask, clusterid,
               cdata, cmask, transpose, method[0]);

    if (ok) {
            /* ------------------------
             * Convert generated C matrices to Perl matrices
             */
            cdata_ref = matrix_c2perl_dbl(aTHX_ cdata, cnrows, cncols);
            cmask_ref = matrix_c2perl_int(aTHX_ cmask, cnrows, cncols);

            /* ------------------------
             * Push the new Perl matrices onto the return stack
             */
            XPUSHs(sv_2mortal( cdata_ref   ));
            XPUSHs(sv_2mortal( cmask_ref   ));
    }

    /* ------------------------
     * Free what we've malloc'ed 
     */
    free_matrix_int(mask,     nrows);
    free_matrix_dbl(matrix,   nrows);
    free_matrix_int(cmask,    cnrows);
    free_matrix_dbl(cdata,    cnrows);
    free(clusterid);

    if (!ok) {
        croak("memory allocation failure in _clustercentroids\n");
    }
    /* Finished _clustercentroids() */

void
_distancematrix(nrows,ncols,data_ref,mask_ref,weight_ref,transpose,dist)
    int      nrows;
    int      ncols;
    SV *     data_ref;
    SV *     mask_ref;
    SV *     weight_ref;
    int      transpose;
    char *   dist;

    PREINIT:
    SV  *    matrix_ref;
    int      nobjects;
    int      ndata;

    double ** data;
    int    ** mask;
    double  * weight;
    double ** matrix;

    int       ok;


    PPCODE:
    /* ------------------------
     * Don't check the parameters, because we rely on the Perl
     * caller to check most parameters.
     */

    /* ------------------------
     * Malloc space for the return values from the library function
     */
        if (transpose==0) {
        nobjects = nrows;
        ndata = ncols;
    } else {
        nobjects = ncols;
        ndata = nrows;
    }

    /* ------------------------
     * Convert data and mask matrices and the weight array
     * from C to Perl.  Also check for errors, and ignore the
     * mask or the weight array if there are any errors. 
     */
    ok = malloc_matrices( aTHX_
        weight_ref, &weight, ndata, 
        data_ref,   &data,
        mask_ref,   &mask,  
        nrows,      ncols
    );
    if (!ok) {
            croak("failed to read input data for _distancematrix");
    }

    /* ------------------------
     * Run the library function
     */
        matrix = distancematrix (nrows,
                                 ncols,
                                 data,
                                 mask,
                                 weight,
                                 dist[0],
                                 transpose);

    /* ------------------------
     * Convert generated C matrices to Perl matrices
     */
    matrix_ref  = ragged_matrix_c2perl_dbl(aTHX_ matrix,  nobjects);

    /* ------------------------
     * Push the new Perl matrices onto the return stack
     */
    XPUSHs(sv_2mortal(matrix_ref));

    /* ------------------------
     * Free what we've malloc'ed 
     */
    free_ragged_matrix_dbl(matrix, nobjects);
    free_matrix_int(mask, nrows);
    free_matrix_dbl(data, nrows);
    free(weight);

    /* Finished _distancematrix() */


void
_somcluster(nrows,ncols,data_ref,mask_ref,weight_ref,transpose,nxgrid,nygrid,inittau,niter,dist)
    int      nrows;
    int      ncols;
    SV *     data_ref;
    SV *     mask_ref;
    SV *     weight_ref;
    int      transpose;
    int      nxgrid;
    int      nygrid;
    double   inittau;
    int      niter;
    char *   dist;

    PREINIT:
    int      (*clusterid)[2];
    SV *  clusterid_ref;
    double*** celldata;

    double  * weight;
    double ** matrix;
    int    ** mask;

    int ok;

    int i;
    AV * matrix_av;
    const int ndata = transpose ? nrows : ncols;
    const int nelements = transpose ? ncols : nrows;

    PPCODE:
    /* ------------------------
     * Don't check the parameters, because we rely on the Perl
     * caller to check most paramters.
     */

    /* ------------------------
     * Allocate space for clusterid[][2]. 
     */
    clusterid = malloc(nelements*sizeof(int[2]));
    if (!clusterid) {
        croak("memory allocation failure in _somcluster\n");
    }
    celldata  =  0;
    /* Don't return celldata, for now at least */


    /* ------------------------
     * Convert data and mask matrices and the weight array
     * from C to Perl.  Also check for errors, and ignore the
     * mask or the weight array if there are any errors. 
     * Set nweights to the correct number of weights.
     */
    ok = malloc_matrices( aTHX_ weight_ref, &weight, ndata, 
                data_ref,   &matrix,
                mask_ref,   &mask,  
                nrows,      ncols);
    if (!ok) {
            croak("failed to read input data for _somcluster\n");
    }

    /* ------------------------
     * Run the library function
     */
    somcluster( 
        nrows, ncols, 
        matrix, mask, weight,
        transpose, nxgrid, nygrid, inittau, niter,
        dist[0], celldata, clusterid
    );

    /* ------------------------
     * Convert generated C matrices to Perl matrices
     */
    matrix_av = newAV();
    for(i=0; i<nelements; ++i) {
        SV* row_ref;
        AV* row_av = newAV();
        av_push(row_av, newSViv(clusterid[i][0]));
        av_push(row_av, newSViv(clusterid[i][1]));
        row_ref = newRV((SV*)row_av);
        av_push(matrix_av, row_ref);
    }
    clusterid_ref = newRV_noinc((SV*)matrix_av);

    /* ------------------------
     * Push the new Perl matrices onto the return stack
     */
    XPUSHs(sv_2mortal( clusterid_ref   ));

    /* ------------------------
     * Free what we've malloc'ed 
     */
    free_matrix_int(mask,     nrows);
    free_matrix_dbl(matrix,   nrows);
    free(weight);
    free(clusterid);

    /* Finished _somcluster() */


void
_pca(nrows, ncols, data_ref)
    int      nrows;
    int      ncols;
    SV * data_ref;

    PREINIT:
    double** u;
    double** v;
    double* w;
    double* m;
    int i;
    int j;
    int nmin;
    int error;
    SV * mean_ref;
    SV * coordinates_ref;
    SV * pc_ref;
    SV * eigenvalues_ref;

    PPCODE:
    if(SvTYPE(SvRV(data_ref)) != SVt_PVAV) { 
        croak("argument to _pca is not an array reference\n");
    }
    nmin = nrows < ncols ? nrows : ncols;
    /* -- Create the output variables -------------------------------------- */
    u = parse_data(aTHX_ data_ref, NULL);
    w = malloc(nmin*sizeof(double));
    v = malloc(nmin*sizeof(double*));
    m = malloc(ncols*sizeof(double));
    if (v) {
        for (i = 0; i < nmin; i++) {
            v[i] = malloc(nmin*sizeof(double));
            if (!v[i]) break;
        }
        if (i < nmin) { /* then we encountered the break */
            while (i-- > 0) free(v[i]);
            free(v);
            v = NULL;
        }
    }
    if (!u || !v || !w || !m) {
        if (u) free(u);
        if (v) free(v);
        if (w) free(w);
        if (m) free(m);
        croak("memory allocation failure in _pca\n");
    }
    /* -- Calculate the mean of each column ------------------------------ */
    for (j = 0; j < ncols; j++) {
        m[j] = 0.0;
        for (i = 0; i < nrows; i++) m[j] += u[i][j];
        m[j] /= nrows;
    }
    /* -- Subtract the mean of each column ------------------------------- */
    for (i = 0; i < nrows; i++)
        for (j = 0; j < ncols; j++)
            u[i][j] -= m[j];
    error = pca(nrows, ncols, u, v, w);
    if (error==0) {
        /* Convert the C variables to Perl variables */
        mean_ref = row_c2perl_dbl(aTHX_ m, ncols);
        if (nrows >= ncols) {
            coordinates_ref = matrix_c2perl_dbl(aTHX_ u, nrows, ncols);
            pc_ref = matrix_c2perl_dbl(aTHX_ v, nmin, nmin);
        }
        else /* nrows < ncols */ {
            pc_ref = matrix_c2perl_dbl(aTHX_ u, nrows, ncols);
            coordinates_ref = matrix_c2perl_dbl(aTHX_ v, nmin, nmin);
        }
        eigenvalues_ref = row_c2perl_dbl(aTHX_ w, nmin);
    }
    for (i = 0; i < nrows; i++) free(u[i]);
    for (i = 0; i < nmin; i++) free(v[i]);
    free(u);
    free(v);
    free(w);
    free(m);
    if (error==-1)
        croak("Insufficient memory for principal components analysis");
    if (error > 0)
        croak("Singular value decomposition failed to converge");
    /* ------------------------
     * Push the new Perl matrices onto the return stack
     */
    XPUSHs(sv_2mortal(mean_ref));
    XPUSHs(sv_2mortal(coordinates_ref));
    XPUSHs(sv_2mortal(pc_ref));
    XPUSHs(sv_2mortal(eigenvalues_ref));
