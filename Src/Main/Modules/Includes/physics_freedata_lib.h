#ifndef physics_freedata_LIB_H
#define physics_freedata_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct PHYSICS_T;

int fd_main(struct PHYSICS_T *const phys);

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

void 
fd_populate_gConf_igConf_dgConf_KerrSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 );

void fd_1st_derivative_Christoffel_symbol(Physics_T *const phys,
                                            const char *const region,
                                            const char *const dChris);

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
fd_populate_psi_alphaPsi_beta_IsoSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta,
 const char *const ig
 );
 
void 
fd_populate_psi_alphaPsi_beta_PGSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta,
 const char *const ig
 );

double *fd_find_KerrSchild_BH_surface(Physics_T *const phys);

void 
fd_populate_psi_alphaPsi_beta_ConfKerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta
 );

void 
fd_populate_alpha_ConfKerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Alpha
 );

void 
fd_populate_alpha_KerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Alpha
 );
 

void 
fd_populate_beta_ConfKerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Beta
 );


void 
fd_populate_beta_KerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Beta
 );

void 
fd_populate_alpha_wKerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Alpha
 );
 
#endif


