#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"

enum Flow
{
  EQUISPACED,// Equispaced collocation
  CHEBYSHEV_ZERO // Chebyshev_zero collocation
};

struct Collocation
{
  double min;
  double max;
  int n;
  double stp;
  double phi_in;
  double phi_fi;
  double a;
  double b;
  enum Flow f;
};

static void make_coords_Cartesian(Patch_T *patch,const int U, enum Flow coll_e);
static void initialize_collocation_struct(Patch_T *patch,struct Collocation *colloc);
static double point(int i, struct Collocation *coll_s);