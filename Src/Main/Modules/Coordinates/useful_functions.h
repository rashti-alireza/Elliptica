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

void point_finder(Needle_T *const needle);
void needle_ex(Needle_T *const needle,const Patch_T *const patch);
void needle_in(Needle_T *const needle,const Patch_T *const patch);
void needle_guess(Needle_T *const needle,const Patch_T *const patch);
void needle_ans(Needle_T *const needle,const Patch_T *const patch);
static void find(Needle_T *const needle,Mode_T mode);
static int IsInside(const double *const x,const double *const lim);
static void fill_limits(double *const lim, const Patch_T *const patch);
int X_of_x(double *const X,const double *const x,const Patch_T *const patch);
int X_of_x_Cartesian_coord(double *const X,const double *const x);
unsigned find_node(const double *const x, const Patch_T *const patch,Flag_T *const flg);
