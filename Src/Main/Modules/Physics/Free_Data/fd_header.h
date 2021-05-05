#ifndef fd_LIB_H
#define fd_LIB_H

#include "core_lib.h"
#include "fields_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "maths_diff_geom_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "physics_lib.h"

/* parameter prefix 
// PLEASE keep it capitalized. */
#define P_ "FREE_DATA_"


/* struct for various functions args. */
struct Analytic_Func_Arg_S
{
  double x,y,z;
  double X,Y,Z,R;
  double dX_D0,dX_D1,dX_D2;
  double dY_D0,dY_D1,dY_D2;
  double dZ_D0,dZ_D1,dZ_D2;
};

/* struct for transition functions */
struct Transition_S
{
 double r;/* independent variable */
 double rmin;/* constant r, like AH radius */
 double rmax;/* constant r, like roll-off rmax */
 double p;/* for example: rolloff power */
 double (*lambda)(struct Transition_S *const ts);/* if need more function */
};

int fd_main(Physics_T *const phys);
void fd_add_fields_gConf_igConf_dgConf(Grid_T *const grid);
void fd_add_fields_ChrisConf_dChrisConf(Grid_T *const grid);
void fd_add_fields_trK_dtrK(Grid_T *const grid);
void fd_add_fields_RicciConf(Grid_T *const grid);
void fd_KerrSchild_set_params(Physics_T *const phys);
void 
fd_populate_gConf_igConf_dgConf_KerrSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 );
void fd_kerr_schild_g_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);
        
void fd_kerr_schild_dg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);
	
void fd_kerr_schild_ddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);

void fd_kerr_schild_dddg_analytic(
        Patch_T *const patch,
        const double BH_center_x,
        const double BH_center_y,
        const double BH_center_z,
        const char *const stem);

void fd_compatible_Christoffel_symbol(Physics_T *const phys,const char *const region,const char *const ig,const char *const dg, const char *const Chris);
void fd_1st_derivative_Christoffel_symbol(Physics_T *const phys,const char *const region,const char *const dChris);
void fd_conformal_Ricci(Physics_T *const phys,
                          const char *const region,
                          const char *const ig,
                          const char *const Chris,
                          const char *const dChris,
                          const char *const RicciConf,
                          const char *const trRicciConf);

void fd_Kij_trK_KerrSchild(Patch_T *const patch,
 const double BH_center_x,const double BH_center_y,
 const double BH_center_z,const char *const ig,
 const char *const Chris,const char *const Kij,
 const char *const trK);


void fd_extrinsic_curvature_KerrSchild(Physics_T *const phys,
                                         const char *const region,
                                         const char *const ig,
                                         const char *const Chris,
                                         const char *const Kij,
                                         const char *const trK,
                                         const char *const dtrK);

void fd_add_fields_MConfIJ(Grid_T *const grid);


void fd_psi_alphaPsi_beta_KerrSchild_patch(Patch_T *const patch,
 const double BH_center_x,const double BH_center_y,
 const double BH_center_z,const char *const ig,
 const char *const Psi,const char *const AlphaPsi,
 const char *const Beta);
 
void 
fd_populate_psi_alphaPsi_beta_KerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta,
 const char *const ig
 );

void 
fd_populate_gConf_igConf_dgConf_IsoSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 );
 
void fd_extrinsic_curvature_Minkowski(Physics_T *const phys,
                                      const char *const region,
                                      const char *const Kij,
                                      const char *const trK,
                                      const char *const dtrK);


void fd_extrinsic_curvature_IsoSchild(Physics_T *const phys,
                                      const char *const region,
                                      const char *const ig,
                                      const char *const Chris,
                                      const char *const Kij,
                                      const char *const trK,
                                      const char *const dtrK);

void 
fd_populate_psi_alphaPsi_beta_IsoSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta,
 const char *const ig/*(inverse metric) if ig is null, it makes them */
 );

void 
fd_populate_gConf_igConf_dgConf_flat
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 );
 
void 
fd_populate_gConf_igConf_dgConf_PGSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 );


void fd_extrinsic_curvature_PGSchild(Physics_T *const phys,
                                      const char *const region,
                                      const char *const ig,
                                      const char *const Chris,
                                      const char *const Kij,
                                      const char *const trK,
                                      const char *const dtrK);

void 
fd_populate_psi_alphaPsi_beta_PGSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta,
 const char *const ig/*(inverse metric) if ig is null, it makes them */
 );

void fd_Kij_trK_PGSchild(Patch_T *const patch,
 const char *const dbeta,const char *const Kij,
 const char *const trK);

void fd_beta_and_dbeta_PGSchild(Physics_T *const phys,
                                const char *const region,
                                const char *const beta,
                                const char *const dbeta);
  
double *fd_find_KerrSchild_BH_surface(Physics_T *const phys);

void 
fd_populate_gConf_igConf_dgConf_ConfKerrSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *g/* given metric stem (if is null it makes it)*/,
 const char *const gConf/* stem of conformal metric */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 );
 
void 
fd_populate_psi_alphaPsi_beta_ConfKerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta
 );


void fd_alpha_KerrSchild_patch(Patch_T *const patch,
 const double BH_center_x,const double BH_center_y,
 const double BH_center_z,const char *const Alpha);

void fd_beta_KerrSchild_patch(Patch_T *const patch,
 const double BH_center_x,const double BH_center_y,
 const double BH_center_z,const char *const ig,
 const char *const Beta);

void 
fd_populate_beta_KerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Beta
 );

void 
fd_populate_beta_ConfKerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Beta
 );

void 
fd_modify_gConf_igConf_dgConf_to_w1flat_w2KS
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 );

void 
fd_modify_trK_to_wtrK_compute_dtrK
 (
 Physics_T *const phys,
 const char *const region,
 const char *const trK,
 const char *const dtrK
 );

void fd_trace_extrinsic_curvature_zero(Physics_T *const phys,
                                       const char *const region,
                                       const char *const trK,
                                       const char *const dtrK);

#endif

