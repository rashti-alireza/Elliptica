#include "core_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "fields_lib.h"
#include "maths_spectral_methods_lib.h"
#include "maths_special_functions_lib.h"

#define Tx(i,x) T_cheb_combined(n[0],i,x)
#define Ty(j,y) T_cheb_combined(n[1],j,y)
#define Tz(k,z) T_cheb_combined(n[2],k,z)

void spectral_filter(const spectral_filter_T *const args);
double T_cheb_combined(const Uint n,const Uint i,const double x);
static void filter_erfclog(const spectral_filter_T *const args);
static double erfclog(const double th, const int p/* e.g., p = 8 */);
