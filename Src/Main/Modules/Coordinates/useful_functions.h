#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "maths_general_lib.h"

typedef enum MODE_T
{
  GUESS,
  FORCE_IN
}Mode_T;

enum Limit
{
  MIN0 = 0,
  MAX0,
  MIN1,
  MAX1,
  MIN2,
  MAX2,
  TOT_Limit
};

static void IsConsistent(Needle_T *needle);
static void find(Needle_T *needle,Mode_T mode);
static int IsInside(double *x,double *lim);
static void fill_limits(double *lim, Patch_T *patch);
int X_of_x(double *X,double *x,Patch_T *patch);
int X_of_x_Cartesian_coord(double *X,double *x,Patch_T *patch);
void needle_in(Needle_T *needle,Patch_T *patch);
