#include "coordinate_shared_lib.h"

void make_nodes_Cartesian_coord(Patch_T *const patch);
void make_JacobianT_Cartesian_coord(Patch_T *const patch);
double JT_Cartesian_patch(const Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double dN0_dx_Cartesian_patch(const Patch_T *const patch,const double *const X);
double dN0_dy_Cartesian_patch(const Patch_T *const patch,const double *const X);
double dN0_dz_Cartesian_patch(const Patch_T *const patch,const double *const X);
double dN1_dx_Cartesian_patch(const Patch_T *const patch,const double *const X);
double dN1_dy_Cartesian_patch(const Patch_T *const patch,const double *const X);
double dN1_dz_Cartesian_patch(const Patch_T *const patch,const double *const X);
double dN2_dx_Cartesian_patch(const Patch_T *const patch,const double *const X);
double dN2_dy_Cartesian_patch(const Patch_T *const patch,const double *const X);
double dN2_dz_Cartesian_patch(const Patch_T *const patch,const double *const X);
