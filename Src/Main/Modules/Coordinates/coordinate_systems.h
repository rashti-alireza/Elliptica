#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"

struct Collocation_s
{
  double min;
  double max;
  unsigned n;
  double stp;
  double a;
  double b;
  Collocation_T f;
};


int fill_nodes(Grid_T *const grid);
static void make_coords_Cartesian_coord(Patch_T *const patch);
static void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const colloc);
static double point(const unsigned i, const struct Collocation_s *const coll_s);
