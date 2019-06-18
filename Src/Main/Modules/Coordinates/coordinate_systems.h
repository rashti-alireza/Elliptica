#include "coordinate_shared_lib.h"

struct R_max
{
    double R_max_l;
    double R_max_r;
};

int make_nodes(Grid_T *const grid);
int make_JacobianT(Grid_T *const grid);
void make_nodes_Cartesian_coord(Patch_T *const patch);
void make_nodes_Spherical_coord(Patch_T *const patch);
void make_nodes_CubedSpherical_coord(Patch_T *const patch);
void make_JacobianT_Cartesian_coord(Patch_T *const patch);
void make_JacobianT_CubedSpherical_coord(Patch_T *const patch);
void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const colloc,const unsigned dir);
double point_value(const unsigned i, const struct Collocation_s *const coll_s);
double dq2_dq1(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
static double dN_dq(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
static double dN_dX(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double General2ChebyshevExtrema(const double X,const unsigned dir,const Patch_T *const patch);
void grid_characteristics_example(Grid_T *const grid);
static void characteristics_Cartesian_grid(Grid_T *const grid);
static void characteristics_BNS_Spherical_grid(Grid_T *const grid);
static void NS_radii_BNS_Spherical_grid(Grid_T *const grid,void *vp);
static void characteristics_BNS_CubedSpherical_grid(Grid_T *const grid);
static void NS_surface_BNS_CubedSpherical_grid(Grid_T *const grid);
static void characteristics_BBN_CubedSpherical_grid(Grid_T *const grid);
static void NS_BH_surface_CubedSpherical_grid(Grid_T *const grid,const double R_NS_l,const double R_BH_r,const double a_BH);
