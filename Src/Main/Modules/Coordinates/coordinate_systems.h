#include "coordinate_shared_lib.h"

int make_nodes(Grid_T *const grid);
int make_JacobianT(Grid_T *const grid);
void make_nodes_Cartesian_coord(Patch_T *const patch);
void make_nodes_Spherical_coord(Patch_T *const patch);
void make_nodes_CubedSpherical_coord(Patch_T *const patch);
void make_nodes_ProjectiveHemisphereUp_coord(Patch_T *const patch);
void make_nodes_ProjectiveHemisphereDown_coord(Patch_T *const patch);
void make_nodes_StereographicSphereLeft_coord(Patch_T *const patch);
void make_nodes_StereographicSphereRight_coord(Patch_T *const patch);
void make_JacobianT_Cartesian_coord(Patch_T *const patch);
void make_JacobianT_ProjectiveHemisphere_coord(Patch_T *const patch);
void make_JacobianT_StereographicSphereRight_coord(Patch_T *const patch);
void make_JacobianT_StereographicSphereLeft_coord(Patch_T *const patch);
void make_JacobianT_CubedSpherical_coord(Patch_T *const patch);
void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const colloc,const unsigned dir);
double point_value(const unsigned i, const struct Collocation_s *const coll_s);
double dq2_dq1(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
static double dN_dq(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
static double dN_dX(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double General2ChebyshevExtrema(const double X,const unsigned dir,const Patch_T *const patch);
