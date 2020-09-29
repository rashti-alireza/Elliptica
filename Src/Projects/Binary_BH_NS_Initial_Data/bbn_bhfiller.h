#include "bbn_headers.h"
#include "utilities_lib.h"
#include "maths_approximation_lib.h"
#include "maths_analytic_lib.h"
#include "maths_equation_solvings_lib.h"

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

/* _gamma inverse */
#define COMPUTE_gammaI(a00,a01,a02,a10,a11,a12,a20,a21,a22) \
  { \
  _gammaI_U0U0[ijk] = (a11*a22 - a12*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  _gammaI_U0U1[ijk] = (-a01*a22 + a02*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  _gammaI_U0U2[ijk] = (a01*a12 - a02*a11)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  _gammaI_U1U1[ijk] = a00*(a00*a22 - a02*a20)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  _gammaI_U1U2[ijk] =-a00*(a00*a12 - a02*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  _gammaI_U2U2[ijk] = a00*(a00*a11 - a01*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  }

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
struct BHFiller_S* 
bhf_init
  (
  Grid_T *const grid/* the whole grid */,
  const char *const method/* the method to be used for extrapolating */
  );

/* free bhfiller struct */
static void bhf_free(struct BHFiller_S *const bhf);

int 
bbn_bhfiller
  (
  Grid_T *const grid/* the whole grid */,
  const char *const method/* the method to be used for extrapolating */
  );

static void 
find_XYZ_and_patch_of_theta_phi_BH_CS
 (
 double *const X,Patch_T **const ppatch,
 const double theta,const double phi,Grid_T *const grid
 );

static int bhf_ChebTn_Ylm(struct BHFiller_S *const bhf);
static int bhf_WTGR(struct BHFiller_S *const bhf);
static double interpolate_from_patch_prim(const char *const field,const double *const X,Patch_T *const patch);
void bbn_bam_error(const char *const msg,const char *const file,const int line);
static void collect_names(struct BHFiller_S *const bhf,const unsigned nf);
static int bhf_4th_Poly_Ylm(struct BHFiller_S *const bhf);
static int bhf_ell_Brown(struct BHFiller_S *const bhf);
void *bbn_bhf_bc_Brown(void *vp1,void *vp2);
void *bbn_bhf_eq_Brown(void *vp1,void *vp2);
void *bbn_bhf_jacobian_bc_Brown(void *vp1,void *vp2);
void *bbn_bhf_jacobian_eq_Brown(void *vp1,void *vp2);
static int find_X_and_patch_outside_BH(const double *const x,const char *const hint,Grid_T *const grid,double *const X,Patch_T **const ppatch);

static void 
bbn_bhf_ell_Brown_field_update
  (
  Patch_T *const patch,
  const char *const name
  );
  
static void 
bbn_bhf_ell_Brown_source_update
  (
  Grid_T *const grid,
  const char *const name
  );


static double punc_psi(void *const params);
static double punc_eta(void *const params);
static double punc_K(void *const params);
static double punc_Beta_U0(void *const params);
static double punc_Beta_U1(void *const params);
static double punc_Beta_U2(void *const params);
static double punc_gamma_D0D0(void *const params);
static double punc_gamma_D0D1(void *const params);
static double punc_gamma_D0D2(void *const params);
static double punc_gamma_D1D1(void *const params);
static double punc_gamma_D1D2(void *const params);
static double punc_gamma_D2D2(void *const params);
double bbn_bhf_smoother(const double r, const double rmax,const double rmin);
static double polynomial5(const double r, const double rmax,const double rmin);
static double polynomial7(const double r, const double rmax,const double rmin);
double bbn_bhf_poly_smoother(const double r,const double rmax,const double rmin);



