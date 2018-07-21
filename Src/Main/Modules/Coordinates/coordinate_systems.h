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


int make_nodes(Grid_T *const grid);
int make_JacobianT(Grid_T *const grid);
static void make_nodes_Cartesian_coord(Patch_T *const patch);
static void make_JacobianT_Cartesian_coord(Patch_T *const patch);
static void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const colloc,const unsigned dir);
static double point(const unsigned i, const struct Collocation_s *const coll_s);
double dq2_dq1(const Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
static double dN_dq(const Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
static double dN_dX(const Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_Cartesian_patch(const Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);