#include "bbn_headers.h"
#include "utilities_lib.h"
#include "maths_approximation_lib.h"
#include "maths_analytic_lib.h"

/* constants */
#define MAX_STR  (100)
#define MAX_STR2 (200)
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
  unsigned nf;/* number of fields */
  char method[MAX_STR];/* methods: 
                       // O. TnYlm_C2 => f(r,th,ph) = R(r)*T(th,ph) 
                       // and it is asked for C2 continitui along r.
                       // R expanded in Cheb_Tn and T in Ylm. */
  struct
  {
    char f[MAX_STR];/* f */
    char df[3][MAX_STR];/* df/dx */
    char ddf[6][MAX_STR];/* d^2f/dx^2 */
    double *ChebTn_coeffs[4];/* ChebTn coeffs to ensure C2 continuity */
    double *realYlm_coeffs[4];/* Ylm coeffs of ChebTn real part */
    double *imagYlm_coeffs[4];/* Ylm coeffs of ChebTn imag part */
    double f_r0;/* value of the field at r = 0 */
  }**fld;/* field info */
  int (*bhfiller)(struct BHFiller_S *const bhf);/* the method to fill the BH */
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

static int bhf_ChebTnYlm_C2(struct BHFiller_S *const bhf);


