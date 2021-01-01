#ifndef eq_LIB_H
#define eq_LIB_H

#include "core_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "physics_EoS_lib.h"
#include "physics_lib.h"
#include "fields_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_spectral_methods_lib.h"
#include "physics_equation_lib.h"

/* string leng */
#define EQ__STR__LEN0 (33)
#define EQ__STR__LEN1 (99)

/* prefix of internal parameters of this project */
#define P_  "EQ_"

/* backup field name prefix */
#define P_Backup_  P_"backup_"

/* residual field name suffix */
#define EQ_Residual_Suffix  "_residual"

/* define paramter and prefix used in EQ_Set_Prefix and EQ_PrefixIt*/
#define EQ_Def_Param_Prefix_Char \
char EQ__param__prefix[EQ__STR__LEN0] = {'\0'}; \
char EQ__param[EQ__STR__LEN1] = {'\0'};\
char EQ__Temp1[EQ__STR__LEN0] = {'\0'};\
char EQ__Temp2[EQ__STR__LEN0] = {'\0'};

/* define paramter prefix for given object X,
// this will be used by EQ_PrefixIt macro.
// the purpose of this macro is to find the correct
// index when multiple objects involved and 
// there are some paramters for each object for instance BHBH
// system which has BH1 and BH2, in this case, for example,
// if X is "BH" and patch->name is "grid1_BH2_around_IB_X0Y0Z0"
// EQ__param__prefix will be "BH2" according to the patch->name.
// if given X has already index then EQ__param__prefix = X
// this is working because while setting up the region on each equation 
// only pertinent patches have been chosen. 
// note: if X = system, it prefixes with system like "BHNS_" 
// for BHNS system. */
#define EQ_Set_Prefix(X) \
sprintf(EQ__Temp1,"%s1",X);\
sprintf(EQ__Temp2,"%s2",X);\
if (strchr(X,'1') || strchr(X,'2'))\
  sprintf(EQ__param__prefix,"%s",X);\
else if (strstr(patch->name,EQ__Temp1))\
  sprintf(EQ__param__prefix,"%s1",X);\
else if (strstr(patch->name,EQ__Temp2))\
  sprintf(EQ__param__prefix,"%s2",X);\
else if (strcmp_i("system",X))\
  sprintf(EQ__param__prefix,"%s",Pgets(P_"system_prefix"));\
else \
  sprintf(EQ__param__prefix,"%s",X);

/* prefix given parameter X by EQ_prefix_par defined in 
// EQ_Set_Prefix. */
#define EQ_PrefixIt(X) \
  (sprintf(EQ__param,"%s_%s",EQ__param__prefix,X) ? EQ__param : NULL)
  
void eq_solve_elliptic_equation(Physics_T *const phys);


#endif

