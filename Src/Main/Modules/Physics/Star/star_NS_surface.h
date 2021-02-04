#include "star_header.h"
#include "maths_spectral_methods_lib.h"
#include "maths_special_functions_lib.h"
#include "maths_linear_algebra_lib.h"

/* constants */
#define MAX_STR  (100)
#define MAX_STR2 (200)
#define MAX_STR_LARGE (1000)

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

/* global variable for this file */
static int Verbose = 0;

/* root finder struct for NS surface eq */
struct NS_surface_RootFinder_S
{
  Patch_T **hpatches;/* all possible patches to find h */
  Uint Nph;/* number of hpatches */
  void *root_finder;
  double x0[3];/* (x,y,z) at the surface */
  double *N;/* the direction of increasing or decreasing of x = x0+N*d */
  //double Euler_C;/* Euler equation const. */
  //double scale;/* to avoid long step in root finder */
  //double maxR;/* max R allowed for NS around */
};

/* struct for general purposes */
struct Demand_S
{
 double r0,r1,r;/* for extrapolation. */
 double fr0,fr1,dfr0,ddfr0;/* f(r0),f(r1), df(r0)/dr, d^2f(r0)/dr^2 */
};

/* all needed items for extrapolation function */
struct Extrap_S
{
  Physics_T *phys;/* the physics of interest */
  Grid_T *grid;/* the grid */
  Patch_T **patches_out;/* patches outside where extrapolation takes place */
  Patch_T **patches_in;/* patches inside where fields already exist */
  Uint npi;/* number of patches inside  */
  Uint npo;/* number of patches outside */
  Uint lmax;/* max l in Ylm expansion */
  Uint Ntheta;/* number of points in theta direction */
  Uint Nphi;/* number of points in phi direction */
  Uint NCoeffs;/* number of coeffs for in extrapolant function */
  Uint nf;/* number of fields */
  char method[MAX_STR];/* the specified methods */
  Uint C2: 1;/* if C2 == 1 means, second order derivatives needed */
  struct
  {
    char f[MAX_STR];/* f */
    char df[3][MAX_STR];/* df/dx */
    char ddf[6][MAX_STR];/* d^2f/dx^2 */
    double (*func_r0)(void *const params);/* function as we close to r0 */
    Uint did_add_df:1;/* 1 if automatically adds required df fields  , otherwise 0. */
    Uint did_add_ddf:1;/* 1 if automatically adds required ddf fields, otherwise 0. */
  }**fld;/* field info */
  double (*extrap)(struct Demand_S *const demand);/* function used for extrapolation (approximation) */
  int (*fmain)(struct Extrap_S *const extap);/* call this to extrapolate  */
};


int 
star_NS_extrapolate
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method to be used for extrapolating */
  );
  

static struct Extrap_S* 
extrap_init
  (
  Physics_T *const phys/* physics of interest */,
  const char **fields_name/* ends determined by 0 */,
  const char *const method/* the method or instruction 
                          // to be used for extrapolating */
  );

static void collect_names(struct Extrap_S *const extrap,
                          const char **const fields_name,
                          const Uint nf);
static void extrap_free(struct Extrap_S *const extrap);
static int fmain_f_df_ddf_CS(struct Extrap_S *const extrap);
static double approx_exp2(struct Demand_S *const demand);
static double approx_poly2(struct Demand_S *const demand);
static void find_NS_surface_Ylm_bisect_CS(Physics_T *const phys);
static double NS_surface_enthalpy_root_finder_eq(void *params,const double *const x);
static double NS_surface_denthalpy_dr_root_finder(void *params,const double *const x,const Uint dir)  __attribute__((unused));
static void find_NS_surface_perfect_s2(Physics_T *const phys);
double star_NS_mass_shedding_indicator(Physics_T *const phys);
int star_NS_find_star_surface(Physics_T *const phys);
void star_start_off_TOV (Physics_T *const phys);
static double approx_inverse_r2(struct Demand_S *const demand);
static double approx_inverse_r2_expmr(struct Demand_S *const demand);
static double approx_inverse_r2_expmAr(struct Demand_S *const demand);
static int extrapolate_expmr_C0_CS(struct Extrap_S *const extrap);
static double approx_inverse_r_expmAr(struct Demand_S *const demand);
static double approx_inverse_r_expmr(struct Demand_S *const demand);





