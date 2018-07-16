#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "memory_managing_lib.h"

struct Collocation_s
{
  double min;
  double max;
  unsigned n;
  double stp;
  double a;
  double b;
  Collocation_T c;
};


int fill_nodes(Grid_T *const grid);
static void make_coords_Cartesian_coord(Patch_T *const patch);
static void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const colloc,const unsigned dir);
static double point(const unsigned i, const struct Collocation_s *const coll_s);
double *make_normalized_collocation_1d(const Patch_T *const patch,const unsigned dir);
