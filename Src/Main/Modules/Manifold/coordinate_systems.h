#include "manifold_header.h"

#define ADJ_00  (a11*a22 - a12*a21)
#define ADJ_01 (-a01*a22 + a02*a21)
#define ADJ_02  (a01*a12 - a02*a11)
#define ADJ_10 (-a10*a22 + a12*a20)
#define ADJ_11  (a00*a22 - a02*a20)
#define ADJ_12 (-a00*a12 + a02*a10)
#define ADJ_20  (a10*a21 - a11*a20)
#define ADJ_21 (-a00*a21 + a01*a20)
#define ADJ_22  (a00*a11 - a01*a10)

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
void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const colloc,const Uint dir);
double point_value(const Uint i, const struct Collocation_s *const coll_s);
double dq2_dq1(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p);
static double dN_dq(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p);
static double dN_dX(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p);
double General2ChebyshevExtrema(const double X,const Uint dir,const Patch_T *const patch);
void grid_characteristics_example(Grid_T *const grid);
static void characteristics_Cartesian_grid_eg(Grid_T *const grid);
static void characteristics_BNS_Spherical_grid_eg(Grid_T *const grid);
static void NS_radii_BNS_Spherical_grid_eg(Grid_T *const grid,void *vp);
static void characteristics_BNS_CS_grid_eg(Grid_T *const grid);
static void NS_surface_BNS_CS_grid_eg(Grid_T *const grid);
static void characteristics_BBN_CS_grid_eg(Grid_T *const grid);
static void NS_BH_surface_CS_grid_eg(Grid_T *const grid,const double R_NS_l,const double R_BH_r,const double a_BH);
static double dq_dN(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p);
void test_dq_dN(Grid_T *const grid);
static void characteristics_SCS_eg(Grid_T *const grid);

