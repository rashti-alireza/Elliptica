#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"

struct Collocation_s
{
  double min;
  double max;
  int n;
  double stp;
  double phi_in;
  double phi_fi;
  double a;
  double b;
  Collocation_T f;
};

static void make_coords_Cartesian_coord(Patch_T *patch);
static void initialize_collocation_struct(Patch_T *patch,struct Collocation_s *colloc);
static double point(int i, struct Collocation_s *coll_s);