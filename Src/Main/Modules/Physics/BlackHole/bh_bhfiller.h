#include "bh_header.h"
#include "utilities_lib.h"
#include "maths_spectral_methods_lib.h"
#include "maths_special_functions_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_linear_algebra_lib.h"

/* constants */
#define MAX_STR  (100)
#define MAX_STR2 (200)
#define MAX_STR_LARGE (1000)
#define MAX_COEFFS (20)

/* 2 indices */
#define IJ(i,j,n)  ((j)+(i)*(n))
/* 2 indices with symmetry for n = 3 */
#define IJsymm3(i,j)  \
 ((j)>=(i) ? \
  (5*(((j)+(i)*2)/6)+((j)+(i)*2)%6) : \
  (5*(((i)+(j)*2)/6)+((i)+(j)*2)%6))

/* math */
#define Power(a,b) pow(a,b)
#define Sqrt(a) sqrt(a)

/* all needed items for bhfiller function */
struct BHFiller_S
{
  Grid_T *grid;/* the grid */
  Patch_T **patches_outBH;/* patches outside the BH */
  Patch_T **patches_inBH;/* patches inside the BH */
  unsigned npi;/* number of patches inside the BH */
  unsigned npo;/* number of patches outside the BH */
  unsigned lmax;/* max l in Ylm expansion */
  unsigned Ntheta;/* number of points in theta direction */
  unsigned Nphi;/* number of points in phi direction */
  unsigned NCoeffs;/* number of coeffs for in extrapolant function */
  unsigned nf;/* number of fields */
  char method[MAX_STR];/* methods: 
                       // O. TnYlm_C2 => f(r,th,ph) = R(r)*T(th,ph) 
                       // and it is asked for C2 continitui along r.
                       // R expanded in Cheb_Tn and T in Ylm.
                       //
                       // O. WTGR => it's C1 in which using Tanh function.
                       //
                       // O. EllEq => it's C^n in which elliptic equation
                       // of n-th order is solved. */
  struct
  {
    char f[MAX_STR];/* f */
    char df[3][MAX_STR];/* df/dx */
    char ddf[6][MAX_STR];/* d^2f/dx^2 */
    double *radial_coeffs[MAX_COEFFS];/* radial coeffs to ensure continuity */
    double *realYlm_coeffs[MAX_COEFFS];/* Ylm coeffs of radial real part */
    double *imagYlm_coeffs[MAX_COEFFS];/* Ylm coeffs of radial imag part */
    double f_r0;/* value of the field at r = 0 */
    double (*func_r0)(void *const params);/* function as we close to r = 0 */
    unsigned did_add_df:1;/* 1 if automatically adds required df fields  , otherwise 0. */
    unsigned did_add_ddf:1;/* 1 if automatically adds required ddf fields, otherwise 0. */
  }**fld;/* field info */
  int (*bhfiller)(struct BHFiller_S *const bhf);/* the method to fill the BH */
};

/* parameters for func_r0 */
struct Param_S
{
  double r;/* radius */
  double rfill;/* BH radius */
  double M;/* M BH */
  double eps;/* cut off */
};

/* initialize the bhfiller struct */
static struct BHFiller_S* 
bhf_init
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  );

/* free bhfiller struct */
static void bhf_free(struct BHFiller_S *const bhf);

int 
bh_black_hole_filler
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  );
  
int 
bh_bhfiller
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  );
  
  
int bh_fill_inside_black_hole(Physics_T *const phys);  
static int bhf_ChebTn_Ylm_pefect_S2_CS(struct BHFiller_S *const bhf);
double bh_bhf_poly_smoother(const double r,const double rmax,const double rmin);
double bh_bhf_smoother(const double r, const double rmax,const double rmin);
static void collect_names(struct BHFiller_S *const bhf,const char **const fields_name,const unsigned nf);
static double polynomial5(const double r, const double rmax,const double rmin);
static double polynomial7(const double r, const double rmax,const double rmin);
static double punc_psi(void *const params);
static double punc_eta(void *const params);
static double punc_K(void *const params);
static double punc_Beta_U0(void *const params);
static double punc_Beta_U1(void *const params);
static double punc_Beta_U2(void *const params);
static double punc_gConf_D0D0(void *const params);
static double punc_gConf_D0D1(void *const params);
static double punc_gConf_D0D2(void *const params);
static double punc_gConf_D1D1(void *const params);
static double punc_gConf_D1D2(void *const params);
static double punc_gConf_D2D2(void *const params);






