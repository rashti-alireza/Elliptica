#include "fd_header.h"
#include "maths_linear_algebra_lib.h"

void 
fd_populate_gConf_dgConf_igConf_KerrSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 );
void fd_compatible_Christoffel_symbol(Physics_T *const phys,const char *const region,const char *const ig,const char *const dg, const char *const Chris);
void fd_1st_derivative_Christoffel_symbol(Physics_T *const phys,const char *const region,const char *const dChris);
void fd_conformal_Ricci(Physics_T *const phys,
                          const char *const region,
                          const char *const ig,
                          const char *const Chris,
                          const char *const dChris,
                          const char *const RicciConf,
                          const char *const trRicciConf);

void fd_extrinsic_curvature_KerrSchild(Physics_T *const phys,
                                         const char *const region,
                                         const char *const ig,
                                         const char *const Chris,
                                         const char *const Kij,
                                         const char *const trK,
                                         const char *const dtrK);



                                         
void 
fd_populate_psi_alphaPsi_beta_KerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta,
 const char *const ig/*(inverse metric) if ig is null, it makes them */
 );
 
 