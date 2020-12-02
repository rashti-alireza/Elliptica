#ifndef physics_freedata_LIB_H
#define physics_freedata_LIB_H

/* forward declaration */
struct PHYSICS_T;

int fd_main(struct PHYSICS_T *const phys);

void 
fd_populate_psi_alpha_beta_KerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const Alpha,
 const char *const Beta,
 const char *const ig/*(inverse metric) if ig is null, it makes them */
 );


#endif


