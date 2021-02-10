#include "bh_header.h"
#include "utilities_lib.h"
#include "maths_spectral_methods_lib.h"
#include "maths_special_functions_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_linear_algebra_lib.h"

/* constants */
#define MAX_STR0  (50)
#define MAX_STR   (99)
#define MAX_STR2  (199)
#define MAX_STR_LARGE (999)
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

/* prefix the given string str with the given prefix pre and 
// return string name. */
#define PrefixIt(name,pre,str) (sprintf(name,"%s_%s",pre,str) ? name : 0)


/* global variable for this file */
static int Verbose = 0;

/* struct for general purposes */
struct Demand_S
{
 double r0,r1,r;/* for extrapolation. */
 double fr0,fr1,dfr0,ddfr0;/* f(r0),f(r1), df(r0)/dr, d^2f(r0)/dr^2 */
};

/* all needed items for bhfiller function */
struct BHFiller_S
{
  Physics_T *phys;
  Grid_T *grid;/* the grid */
  Patch_T **patches_outBH;/* patches outside the BH */
  Patch_T **patches_inBH;/* patches inside the BH */
  Uint npi;/* number of patches inside the BH */
  Uint npo;/* number of patches outside the BH */
  Uint lmax;/* max l in Ylm expansion */
  Uint Ntheta;/* number of points in theta direction */
  Uint Nphi;/* number of points in phi direction */
  Uint NCoeffs;/* number of coeffs for in extrapolant function */
  Uint nf;/* number of fields */
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
    Uint did_add_df:1;/* 1 if automatically adds required df fields  , otherwise 0. */
    Uint did_add_ddf:1;/* 1 if automatically adds required ddf fields, otherwise 0. */
  }**fld;/* field info */
  int (*bhfiller)(struct BHFiller_S *const bhf);/* the method to fill the BH */
  double (*extrap)(struct Demand_S *const demand);/* the function used for extrapolation (approximation) */
  Uint C1: 1;/* if C1 == 1 means, first order derivatives needed */
  Uint C2: 1;/* if C2 == 1 means, second order derivatives needed */
};


/* parameters for func_r0 */
struct Param_S
{
  double r;/* radius */
  double rfill;/* BH radius */
  double M;/* M BH */
  double eps;/* cut off */
};


int bh_add_patch_inside_black_hole(Physics_T *const phys,const char *const region);

/* initialize the bhfiller struct */
static struct BHFiller_S* 
bhf_init
  (
  Physics_T *const phys/* physics of interest */,
  char **const fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  );

/* free bhfiller struct */
static void bhf_free(struct BHFiller_S *const bhf);

int 
bh_bhfiller
  (
  Physics_T *const phys/* physics of interest */,
  char **const fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  );
  

void 
bh_interpolating_fields_on_a_line
  (
  Physics_T *const phys/* physics of interest */,
  const char *const sfields_name/* comma separated fields */,
  const char *const dir/* output directory */,
  const char *const stem_g/* if stem of a metric given => test det(g) > 0 */
  );
  
int bh_fill_inside_black_hole(Physics_T *const phys);  
static int bhf_ChebTn_Ylm_pefect_S2_CS(struct BHFiller_S *const bhf);
static int bhf_ChebTn_general_S2_CS(struct BHFiller_S *const bhf);
double bh_bhf_poly_smoother(const double r,const double rmax,const double rmin);
double bh_bhf_smoother(const double r, const double rmax,const double rmin);
static void collect_names(struct BHFiller_S *const bhf,char **const fields_name,const Uint nf);
static double polynomial5(const double r, const double rmax,const double rmin);
static double polynomial7(const double r, const double rmax,const double rmin);
static int bhf_f_df_ddf_perfect_s2_CS(struct BHFiller_S *const bhf);
static double approx_expmr_C0(struct Demand_S *const demand);
static double approx_r_expmr_C1(struct Demand_S *const demand);









