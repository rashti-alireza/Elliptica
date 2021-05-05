#include "fd_header.h"
#include "maths_linear_algebra_lib.h"
#include "physics_transformation_lib.h"

#define STR_LEN (99)

/* handy macro for transition function */
#define SET_TRANSITION_FUNC_BH_TYPE0  \
  const double r_min = Getd("min_radius");\
  const double r_max = Getd("RollOff_rmax");\
  const double p_att = Getd("RollOff_power");\
  double (*transit)(struct Transition_S *const ts)  = 0;\
  double (*lambda)(struct Transition_S *const ts)   = 0;\
  \
  assert(r_min > 0); assert(r_max > 0);\
  \
  IF_sval("RollOff_function","exp(-lambda*(r/rmax)^p):r<rmax")\
   {transit = f_exp_type1;}\
  else IF_sval("RollOff_function","exp(-lambda*(r/rmax)^p)")\
   {transit = f_exp_type2;}\
  else\
   {Error0("No such option for RollOff_function.");}\
  \
  IF_sval("RollOff_lambda","|(r-rmin)/(rmax-r)|")\
   {lambda = f_ratio_type1;}\
  else IF_sval("RollOff_lambda","constant_1")\
   {lambda = f_constant_1;}\
  else\
   {Error0("No such option for RollOff_lambda.");}
  

void 
fd_populate_gConf_igConf_dgConf_KerrSchild
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
 const char *const ig
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
 const char *const ig
 );


void fd_beta_and_dbeta_PGSchild(Physics_T *const phys,
                                const char *const region,
                                const char *const beta,
                                const char *const dbeta);

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

void 
fd_populate_alpha_wKerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Alpha
 );

static double f_constant_1(struct Transition_S *const ts);
static double f_ratio_type1(struct Transition_S *const ts);
static double f_exp_type1(struct Transition_S *const ts);
static double f_exp_type2(struct Transition_S *const ts);
void fd_trace_extrinsic_curvature_zero(Physics_T *const phys,
                                       const char *const region,
                                       const char *const trK,
                                       const char *const dtrK);

