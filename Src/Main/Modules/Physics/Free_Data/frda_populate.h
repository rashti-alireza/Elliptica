#include "frda_header.h"
#include "maths_linear_algebra_lib.h"

void 
frda_populate_gConf_dgConf_igConf_KerrSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 );
void frda_compatible_Christoffel_symbol(Physics_T *const phys,const char *const region,const char *const ig,const char *const dg, const char *const Chris);
void frda_1st_derivative_Christoffel_symbol(Physics_T *const phys,const char *const region,const char *const dChris);
void frda_conformal_Ricci(Physics_T *const phys,
                          const char *const region,
                          const char *const ig,
                          const char *const Chris,
                          const char *const dChris,
                          const char *const RicciConf,
                          const char *const trRicciConf);

void frda_extrinsic_curvature_KerrSchild(Physics_T *const phys,
                                         const char *const region,
                                         const char *const ig,
                                         const char *const Chris,
                                         const char *const Kij,
                                         const char *const trK,
                                         const char *const dtrK);



                                         
