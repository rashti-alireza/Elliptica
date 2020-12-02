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


#endif


