#include "bbn_headers.h"
#include "utilities_lib.h"
#include "maths_approximation_lib.h"
#include "maths_analytic_lib.h"

#define MAX_STR  (100)
#define MAX_STR2 (200)

/* 2 indices */
#define IJ(i,j,n)  ((j)+(i)*(n))

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
  }*fld;/* field info */
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

