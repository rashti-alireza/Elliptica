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


int frda_main(Physics_T *const phys);
void frda_add_fields_gConf_dgConf_igConf(Grid_T *const grid);
void frda_add_fields_ChrisConf_dChrisConf(Grid_T *const grid);
void frda_add_fields_trK_dtrK(Grid_T *const grid);
void frda_add_fields_RicciConf(Grid_T *const grid);

#endif

