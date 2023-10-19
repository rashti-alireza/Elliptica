#include "spectral_filter.h"

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
    double th1 = (double)i/n[0];
    for (j = 0; j < n[1]; ++j)
    {
      double th2 = (double)j/n[1];
      for (k = 0; k < n[2]; ++k)
      {
        double th3 = (double)k/n[2];
        C[i_j_k_to_ijk(n,i,j,k)] *= erfclog(th1,p)*erfclog(th2,p)*erfclog(th3,p);
      }
    }
  }
  
  // now rewrite the field
  for (i = 0; i < n[0]; ++i)
  {
    for (j = 0; j < n[1]; ++j)
    {
      for (k = 0; k < n[2]; ++k)
      {
        Uint ijk = i_j_k_to_ijk(n,i,j,k);
        double X = patch->node[ijk]->X[0];
        double Y = patch->node[ijk]->X[1];
        double Z = patch->node[ijk]->X[2];
        v[ijk] = C[ijk]*Tx(i,X)*Ty(j,Y)*Tz(k,Z);
      }
    }
  }
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

