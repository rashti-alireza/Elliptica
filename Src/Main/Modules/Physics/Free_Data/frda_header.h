#ifndef frda_LIB_H
#define frda_LIB_H

#include "core_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "physics_lib.h"
#include "fields_lib.h"

/* parameter prefix */
#define P_ "free_data_"


/* struct for various functions args. */
struct Analytic_Func_Arg_S
{
  double x,y,z;
  double X,Y,Z,R;
  double dX_D0,dX_D1,dX_D2;
  double dY_D0,dY_D1,dY_D2;
  double dZ_D0,dZ_D1,dZ_D2;
};

int frda_main(Physics_T *const phys);
void frda_add_fields_gConf_dgConf_igConf(Grid_T *const grid);
void frda_add_fields_ChrisConf_dChrisConf(Grid_T *const grid);
void frda_add_fields_trK_dtrK(Grid_T *const grid);
void frda_add_fields_RicciConf(Grid_T *const grid);
void frda_ks_free_data_set_params(Physics_T *const phys);


#endif
