#include "spectral_filter.h"

/* a quick tutorial by example:
  // filter rho0 field in the given patch using "erfclog" method
  spectral_filter_T args;
  args.patch     = patch;
  args.field     = "rho0";
  args.filter    = "erfclog";
  args.erfclog_p = 12;
  spectral_filter(&args);
*/

// filter the field of interest from high wave-numbers
void spectral_filter(const spectral_filter_T *const args)
{
  // sanity checks
  assert(args->patch);
  if (strcmp_i(args->filter,"erfclog"))
  {
    filter_erfclog(args);
  }
  else
  {
    Error0(NO_OPTION);
  }
}

/*  u(x) = sum_{j=0}^{N} filter_erfclog(j/N , p)  c_j  base_j(x) */
static void filter_erfclog(const spectral_filter_T *const args)
{
  Patch_T *const patch = args->patch;
  Field_T *const field = patch->fields[Ind(args->field)];
  double *const v      = field->v;
  const Uint *const n  = patch->n;
  const int p = args->erfclog_p; 
  double *C = 0;/* coeffs of expansion of field in Tn */ 
  Uint i,j,k;
  
  make_coeffs_3d(field);
  C = field->v2;

  // adjust coeffs
  for (i = 0; i < n[0]; ++i)
  {
    double th1 = (double)i/(n[0]-1.);
    for (j = 0; j < n[1]; ++j)
    {
      double th2 = (double)j/(n[1]-1.);
      for (k = 0; k < n[2]; ++k)
      {
        double th3 = (double)k/(n[2]-1.);
        C[i_j_k_to_ijk(n,i,j,k)] *= erfclog(th1,p)*erfclog(th2,p)*erfclog(th3,p);
      }
    }
  }

  // now rewrite the field with new coeffs
  for (Uint ijk_ = 0; ijk_ < patch->nn; ++ijk_)
  {
    double X = General2ChebyshevExtrema(patch->node[ijk_]->X[0],0,patch);
    double Y = General2ChebyshevExtrema(patch->node[ijk_]->X[1],1,patch);
    double Z = General2ChebyshevExtrema(patch->node[ijk_]->X[2],2,patch);
    v[ijk_]  = 0.;
    
    for (i = 0; i < n[0]; ++i)
    {
      double tx = Tx(i,X);
      for (j = 0; j < n[1]; ++j)
      {
        double ty = Ty(j,Y);
        for (k = 0; k < n[2]; ++k)
        {
          Uint ijk = i_j_k_to_ijk(n,i,j,k);
          
          v[ijk_] += C[ijk]*tx*ty*Tz(k,Z);
        }
      }
    }
  }// end of for (Uint ijk_ = 0; ijk_ < patch->nn; ++ijk_)
}

/* logerfc filter (adapted from whisky THC)
 *
 *  It is an approximate Vandeven filter, see:
 *
 *  J. B. Boyd, The Erfc-Log Filter and the Asymptotics of the Euler
 *  and Vandeven Sequence Accelerations,In Proceedings
 *  of the Third International Conference on Spectral and
 *  High Order Methods (1996) 267-276
 *
 *  equation (15).
 *
 *  u(x) = sum_{j=0}^{N} filter_erfclog(j/N , p)  c_j  base_j(x)
 * 
 *  note: the smaller the p the wider the Gaussian part of the erfclog.
 */
static double erfclog(const double th, const int p/* e.g., p = 8 */)
{
  const double epsilon = ROUND_OFF_ERR;
  const double th_bar = fabs(th) - 0.5;
  const double den = 4.0*Pow2(th_bar);
  double tmp;

  if (th < epsilon)
  {
    return 1.;
  }
  if (th > 1. - epsilon)
  {
    return 0.;
  }
  if(den < epsilon)
  {
    return 0.5;
  }

  tmp = -log(1. - den);
  tmp /= den;
  tmp = sqrt(tmp);
  tmp *= 2.0*sqrt(p)*th_bar;

  return 0.5*erfc(tmp);
}

