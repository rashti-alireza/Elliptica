/* msimplify(...) */
/* msimplify(...) */
/* msimplify(...) */
#include "fd_KerrSchild_header.h"
KS_func_def_macro(K0)KS_func_args_macro;
KS_func_def_macro(dK0_D2)KS_func_args_macro;
KS_func_def_macro(dK0_D0)KS_func_args_macro;
KS_func_def_macro(dK0_D1)KS_func_args_macro;
KS_func_def_macro(ddK0_D1D1)KS_func_args_macro;
KS_func_def_macro(ddK0_D1D2)KS_func_args_macro;
KS_func_def_macro(ddK0_D0D1)KS_func_args_macro;
KS_func_def_macro(ddK0_D0D0)KS_func_args_macro;
KS_func_def_macro(ddK0_D0D2)KS_func_args_macro;
KS_func_def_macro(ddK0_D2D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D2D2D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D2D2D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D1D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D1D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D1D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D0D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D2D1)KS_func_args_macro;
KS_func_def_macro(K1)KS_func_args_macro;
KS_func_def_macro(dK1_D1)KS_func_args_macro;
KS_func_def_macro(dK1_D0)KS_func_args_macro;
KS_func_def_macro(dK1_D2)KS_func_args_macro;
KS_func_def_macro(ddK1_D2D2)KS_func_args_macro;
KS_func_def_macro(ddK1_D0D0)KS_func_args_macro;
KS_func_def_macro(ddK1_D0D2)KS_func_args_macro;
KS_func_def_macro(ddK1_D1D2)KS_func_args_macro;
KS_func_def_macro(ddK1_D1D1)KS_func_args_macro;
KS_func_def_macro(ddK1_D0D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D1D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D1D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D1D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D0D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D2D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D2D2D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D2D2D2)KS_func_args_macro;
KS_func_def_macro(K2)KS_func_args_macro;
KS_func_def_macro(dK2_D0)KS_func_args_macro;
KS_func_def_macro(dK2_D1)KS_func_args_macro;
KS_func_def_macro(dK2_D2)KS_func_args_macro;
KS_func_def_macro(ddK2_D0D0)KS_func_args_macro;
KS_func_def_macro(ddK2_D2D2)KS_func_args_macro;
KS_func_def_macro(ddK2_D1D2)KS_func_args_macro;
KS_func_def_macro(ddK2_D0D2)KS_func_args_macro;
KS_func_def_macro(ddK2_D0D1)KS_func_args_macro;
KS_func_def_macro(ddK2_D1D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D1D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D0D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D2D2D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D2D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D2D2D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D1D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D1D1)KS_func_args_macro;
KS_func_def_macro(H)KS_func_args_macro;
KS_func_def_macro(dH_D2)KS_func_args_macro;
KS_func_def_macro(dH_D0)KS_func_args_macro;
KS_func_def_macro(dH_D1)KS_func_args_macro;
KS_func_def_macro(ddH_D2D2)KS_func_args_macro;
KS_func_def_macro(ddH_D0D0)KS_func_args_macro;
KS_func_def_macro(ddH_D1D2)KS_func_args_macro;
KS_func_def_macro(ddH_D0D2)KS_func_args_macro;
KS_func_def_macro(ddH_D0D1)KS_func_args_macro;
KS_func_def_macro(ddH_D1D1)KS_func_args_macro;
KS_func_def_macro(dddH_D2D2D2)KS_func_args_macro;
KS_func_def_macro(dddH_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddH_D0D1D1)KS_func_args_macro;
KS_func_def_macro(dddH_D0D0D1)KS_func_args_macro;
KS_func_def_macro(dddH_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddH_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddH_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddH_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddH_D1D2D1)KS_func_args_macro;
KS_func_def_macro(dddH_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddH_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddH_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddH_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddH_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddH_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddH_D1D1D2)KS_func_args_macro;
KS_func_def_macro(dddH_D2D2D1)KS_func_args_macro;
KS_func_def_macro(dddH_D0D1D0)KS_func_args_macro;
KS_func_def_macro(K0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// R
// X
// Y
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))/(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2));
}
KS_func_def_macro(dK0_D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// R
// X
// Y
((pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) - 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2);
}
KS_func_def_macro(dK0_D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// R
// X
// Y
((pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) - 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2);
}
KS_func_def_macro(dK0_D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// R
// X
// Y
((pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) - 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2);
}
KS_func_def_macro(ddK0_D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(a_BH*(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) )) - 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 2*(a_BH*
((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 
8*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(ddK0_D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(a_BH*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) )) - 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 
8*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(ddK0_D0D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(a_BH*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) )) - 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 
8*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(ddK0_D0D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(a_BH*(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) )) - 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + 2*(a_BH*
((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 
8*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(ddK0_D0D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(a_BH*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) )) - 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 
8*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(ddK0_D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(a_BH*(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) )) - 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) + 
8*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(dddK0_D0D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D2D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ) + 2*((fd_ks_dX_D2_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D1D1D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D1D1 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ) + 2*((fd_ks_dX_D1_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 2*(a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + (a_BH*
(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 2*(a_BH*
((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D0D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D0D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D0D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dX_D0_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + 
(a_BH*(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D0D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dX_D0_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dX_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dX_D2_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + (a_BH*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(a_BH*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(a_BH*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 3*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D2D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D1D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D1D2D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + 2*((fd_ks_dX_D2_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D1D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D1D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D1D1D1 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D1D1D1 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 3*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) )) - 2*pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dddR_D1D1D1 KS_func_pass_args_macro ) + 3*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
3*(a_BH*((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
3*(a_BH*((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 
3*(a_BH*(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*
(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )) + 24*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (a_BH*
((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) - 48*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 3))/
pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D2D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D2D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D2D2D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D2D2D2 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 3*(fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) )) - 2*pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dddR_D2D2D2 KS_func_pass_args_macro ) + 3*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
3*(a_BH*((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
3*(a_BH*((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 
3*(a_BH*(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*
(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dX_D2_) ))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro )) + 24*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 48*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 3))/
pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D2D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D2D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ) + 2*((fd_ks_dX_D2_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D0D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D1D1 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ) + 2*((fd_ks_dX_D1_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 2*(a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + (a_BH*
(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 2*(a_BH*
((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D0D1D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D0D1 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dX_D0_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + 
(a_BH*(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D1D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D1D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D1D1D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 2*((fd_ks_dX_D1_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 
(a_BH*(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D0D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dX_D0_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dX_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dX_D2_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + (a_BH*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(a_BH*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(a_BH*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 3*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D1D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D1D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D1D2D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + 2*((fd_ks_dX_D2_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D0D0D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D0D0 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D0D0 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D0D0 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 3*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) )) - 2*pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dddR_D0D0D0 KS_func_pass_args_macro ) + 3*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 
3*(a_BH*((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 
3*(a_BH*((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + 
3*(a_BH*(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*
(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro )) + 24*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (a_BH*
((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) - 48*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 3))/
pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D0D0D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D0D1 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dX_D0_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + 
(a_BH*(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D0D0D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D0D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D0D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dX_D0_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + 
(a_BH*(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D1D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dX_D0_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dX_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dX_D2_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + (a_BH*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(a_BH*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(a_BH*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 3*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D0_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dY_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK0_D1D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(a_BH*(fd_ks_dddY_D1D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddX_D1D1D2 KS_func_pass_args_macro ) + fd_ks_X(x, y, z)*
(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 2*((fd_ks_dX_D1_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_Y(x, y, z) + 
fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
(a_BH*((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 
(a_BH*(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*
(a_BH*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + 
fd_ks_X(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dX_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dX_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(2*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*
fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 3*(a_BH*
fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(a_BH*((fd_ks_dY_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dX_D1_) ) + fd_ks_X(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (a_BH*
((fd_ks_dY_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dX_D2_) ) + 
fd_ks_X(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) - 48*(a_BH*fd_ks_Y(x, y, z) + fd_ks_R(x, y, z)*fd_ks_X(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(K1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// R
// X
// Y
(-a_BH*fd_ks_X(x, y, z) + fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))/(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2));
}
KS_func_def_macro(dK1_D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// R
// X
// Y
((pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2);
}
KS_func_def_macro(dK1_D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// R
// X
// Y
((pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2);
}
KS_func_def_macro(dK1_D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// R
// X
// Y
((pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) + 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2);
}
KS_func_def_macro(ddK1_D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(-a_BH*(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) )) + 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) - 2*(-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) + 
8*(-a_BH*fd_ks_X(x, y, z) + fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(ddK1_D0D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(-a_BH*(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) )) + 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) - 2*(-a_BH*
((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 
8*(-a_BH*fd_ks_X(x, y, z) + fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(ddK1_D0D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(-a_BH*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) )) + 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-
a_BH*((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 
(-a_BH*((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) - 
8*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(ddK1_D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(-a_BH*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) )) + 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-
a_BH*((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 
(-a_BH*((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) - 
8*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(ddK1_D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(-a_BH*(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) )) + 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) - 2*(-a_BH*
((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 
8*(-a_BH*fd_ks_X(x, y, z) + fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(ddK1_D0D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*(-a_BH*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) )) + 2*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*
((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) - (-
a_BH*((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) - 
(-a_BH*((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) - 
8*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 3);
}
KS_func_def_macro(dddK1_D0D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dY_D0_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dY_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dY_D2_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - (-
a_BH*((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
(-a_BH*((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
(-a_BH*((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - (-a_BH*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-
a_BH*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) - (-
a_BH*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - (a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 3*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-
a_BH*((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D0D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D0D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D0D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dY_D0_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*(-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D2_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) - (-a_BH*(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-
a_BH*((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D1D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D1D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D1D1D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D1D1D1 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + 3*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) )) + 2*pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*
((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dddR_D1D1D1 KS_func_pass_args_macro ) + 3*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) - 
3*(-a_BH*((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
3*(a_BH*((fd_ks_dX_D1_) ) - fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) - 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) - 
3*(-a_BH*(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*
(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )) + 24*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) - (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (-a_BH*
((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 48*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 3))/
pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D0D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D2D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ) + 2*((fd_ks_dY_D2_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D0_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) - 2*(-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) - 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D1D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D1D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D1D1D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + 2*((fd_ks_dY_D1_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*(-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D2_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) - (-a_BH*(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-
a_BH*((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D1D1D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D1D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ) + 2*((fd_ks_dY_D1_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D0_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) - 2*(-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - (-a_BH*(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) - 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 2*(-a_BH*
((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D0D1D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D0D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dY_D0_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - 2*(-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D1_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) - (-a_BH*(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*(-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + (-
a_BH*((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D0D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D1D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ) + 2*((fd_ks_dY_D1_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D0_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) - 2*(-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - (-a_BH*(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) - 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 2*(-a_BH*
((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D0D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dY_D0_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dY_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dY_D2_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - (-
a_BH*((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
(-a_BH*((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
(-a_BH*((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - (-a_BH*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-
a_BH*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) - (-
a_BH*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - (a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 3*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-
a_BH*((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D0D0D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D0D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D0D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dY_D0_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*(-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D2_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) - (-a_BH*(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-
a_BH*((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D0D0D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D0D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dY_D0_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - 2*(-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D1_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) - (-a_BH*(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*(-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + (-
a_BH*((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D0D0D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D0D0 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D0D0 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D0D0 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + 3*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) )) + 2*pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*
((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dddR_D0D0D0 KS_func_pass_args_macro ) + 3*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) - 
3*(-a_BH*((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 
3*(a_BH*((fd_ks_dX_D0_) ) - fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) - 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) - 
3*(-a_BH*(fd_ks_ddX_D0D0 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*
(fd_ks_ddY_D0D0 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro )) + 24*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) - (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (-a_BH*
((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 48*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 3))/
pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D2D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D2D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ) + 2*((fd_ks_dY_D2_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D0_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) - 2*(-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) - 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D1D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dY_D0_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dY_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dY_D2_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - (-
a_BH*((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
(-a_BH*((fd_ks_dX_D0_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
(-a_BH*((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - (-a_BH*(fd_ks_ddX_D0D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D0D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-
a_BH*(fd_ks_ddX_D0D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) - (-
a_BH*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - (a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 3*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-a_BH*((fd_ks_dX_D0_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D0_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-
a_BH*((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D1D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D1D1D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D1D1D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + 2*((fd_ks_dY_D1_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 2*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*(-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D1_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D2_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) - (-a_BH*(fd_ks_ddX_D1D1 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D1D1 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + (-
a_BH*((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D1D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D1D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D1D2D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + 2*((fd_ks_dY_D2_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D1_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) - 2*(-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) - 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D2D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D1D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D1D2D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ) + 2*((fd_ks_dY_D2_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro )) + 2*pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 2)*((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*
fd_ks_R(x, y, z)*(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (a_BH*((fd_ks_dX_D1_) ) - 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) - fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) - 2*(-a_BH*((fd_ks_dX_D2_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 2*(-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - (-a_BH*(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) - 2*(-
a_BH*(fd_ks_ddX_D1D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*(fd_ks_ddY_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Y(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dY_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )) + 8*
(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) - 
2*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 3*(a_BH*
fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (-a_BH*((fd_ks_dX_D1_) ) + 
fd_ks_R(x, y, z)*((fd_ks_dY_D1_) ) + fd_ks_Y(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 2*(-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z) + 48*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 3)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2) + 
pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(dddK1_D2D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// X
// Y
(pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 3)*(-a_BH*(fd_ks_dddX_D2D2D2 KS_func_pass_args_macro ) + 
fd_ks_R(x, y, z)*(fd_ks_dddY_D2D2D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*
(fd_ks_dddR_D2D2D2 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + 3*(fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dY_D2_) )) + 2*pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 2)*
((a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_dddR_D2D2D2 KS_func_pass_args_macro ) + 3*(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) - 
3*(-a_BH*((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
3*(a_BH*((fd_ks_dX_D2_) ) - fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) - 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) - 
3*(-a_BH*(fd_ks_ddX_D2D2 KS_func_pass_args_macro ) + fd_ks_R(x, y, z)*
(fd_ks_ddY_D2D2 KS_func_pass_args_macro ) + fd_ks_Y(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dY_D2_) ))*fd_ks_R(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro )) + 24*(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2))*(-
(a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) - (a_BH*fd_ks_X(x, y, z) - fd_ks_R(x, y, z)*
fd_ks_Y(x, y, z))*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (-a_BH*
((fd_ks_dX_D2_) ) + fd_ks_R(x, y, z)*((fd_ks_dY_D2_) ) + 
fd_ks_Y(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 48*(a_BH*fd_ks_X(x, y, z) - 
fd_ks_R(x, y, z)*fd_ks_Y(x, y, z))*pow(fd_ks_R(x, y, z), 3)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 3))/
pow(pow(a_BH, 2) + pow(fd_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(K2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// R
// Z
fd_ks_Z(x, y, z)/fd_ks_R(x, y, z);
}
KS_func_def_macro(dK2_D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// R
// Z
(fd_ks_R(x, y, z)*((fd_ks_dZ_D0_) ) - fd_ks_Z(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 2);
}
KS_func_def_macro(dK2_D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// R
// Z
(fd_ks_R(x, y, z)*((fd_ks_dZ_D1_) ) - fd_ks_Z(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 2);
}
KS_func_def_macro(dK2_D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// R
// Z
(fd_ks_R(x, y, z)*((fd_ks_dZ_D2_) ) - fd_ks_Z(x, y, z)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 2);
}
KS_func_def_macro(ddK2_D0D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(-(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) + pow(fd_ks_R(x, y, z), 2)*
(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))/
pow(fd_ks_R(x, y, z), 3);
}
KS_func_def_macro(ddK2_D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(-(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z) + pow(fd_ks_R(x, y, z), 2)*
(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/
pow(fd_ks_R(x, y, z), 3);
}
KS_func_def_macro(ddK2_D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(-(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z) + pow(fd_ks_R(x, y, z), 2)*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 3);
}
KS_func_def_macro(ddK2_D0D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(-(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) + pow(fd_ks_R(x, y, z), 2)*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 3);
}
KS_func_def_macro(ddK2_D0D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(-(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) + pow(fd_ks_R(x, y, z), 2)*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 3);
}
KS_func_def_macro(ddK2_D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(-(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z) + pow(fd_ks_R(x, y, z), 2)*
(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))/
pow(fd_ks_R(x, y, z), 3);
}
KS_func_def_macro(dddK2_D0D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(2*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D0_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D0D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D0D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D2_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))/
pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D0D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ) + 2*((fd_ks_dZ_D2_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D2D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D0D1D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(2*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D1_) ) + 2*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D0_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D0D1 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D0D0D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(2*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D0_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D0D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D0D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D2_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))/
pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D0D0D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(6*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) - 
(fd_ks_Z(x, y, z)*(fd_ks_dddR_D0D0D0 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 3*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D0D0 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 3))/
pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D0D0D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(2*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D1_) ) + 2*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D0_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D0D1 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D1D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + 2*((fd_ks_dZ_D2_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D1D2D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D2D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(6*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 
(fd_ks_Z(x, y, z)*(fd_ks_dddR_D2D2D2 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 3*(fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D2D2D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 3))/
pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D1D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(2*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D1_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D1D1D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D2D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ) + 2*((fd_ks_dZ_D2_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D2D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D2D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + 2*((fd_ks_dZ_D2_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D1D2D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D1D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D2_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))/
pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D1D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(2*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + (fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D1_) )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D1D1D2 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D1D1D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ) + 2*((fd_ks_dZ_D1_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D1D1 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D1D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(6*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) - 
(fd_ks_Z(x, y, z)*(fd_ks_dddR_D1D1D1 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 3*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D1D1D1 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 3))/
pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(dddK2_D0D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
(2*(fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z) - (fd_ks_Z(x, y, z)*
(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ) + 2*((fd_ks_dZ_D1_) )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 2) + pow(fd_ks_R(x, y, z), 3)*
(fd_ks_dddZ_D0D1D1 KS_func_pass_args_macro ) - 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))/pow(fd_ks_R(x, y, z), 4);
}
KS_func_def_macro(H)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// R
// Z
Lambda*M_BH*pow(fd_ks_R(x, y, z), 3)/(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4));
}
KS_func_def_macro(dH_D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// R
// Z
-Lambda*M_BH*(2*pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 2)/(pow(a_BH, 4)*pow(fd_ks_Z(x, y, z), 4) + 2*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 4)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 8));
}
KS_func_def_macro(dH_D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// R
// Z
-Lambda*M_BH*(2*pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 2)/(pow(a_BH, 4)*pow(fd_ks_Z(x, y, z), 4) + 2*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 4)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 8));
}
KS_func_def_macro(dH_D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// R
// Z
-Lambda*M_BH*(2*pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 2)/(pow(a_BH, 4)*pow(fd_ks_Z(x, y, z), 4) + 2*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 4)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 8));
}
KS_func_def_macro(ddH_D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*
fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) - 4*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D2_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) + 8*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z)/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3);
}
KS_func_def_macro(ddH_D0D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) - 4*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro ) - 2*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) + 8*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z)/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3);
}
KS_func_def_macro(ddH_D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(-2*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 6*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 2)*((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ) + 
3*pow(a_BH, 4)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
6*pow(a_BH, 4)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) - 6*pow(a_BH, 4)*fd_ks_R(x, y, z)*
pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 
6*pow(a_BH, 4)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 6)*
fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 6)*((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ) + 
2*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 5)*
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 5)*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) - 24*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 4)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 
pow(fd_ks_R(x, y, z), 9)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 2*pow(fd_ks_R(x, y, z), 8)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)/
(pow(a_BH, 6)*pow(fd_ks_Z(x, y, z), 6) + 3*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 4)*
pow(fd_ks_Z(x, y, z), 4) + 3*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 8)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 12));
}
KS_func_def_macro(ddH_D0D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(-2*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 6*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 2)*((fd_ks_dZ_D0_) )*((fd_ks_dZ_D2_) ) + 
3*pow(a_BH, 4)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 
6*pow(a_BH, 4)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) - 6*pow(a_BH, 4)*fd_ks_R(x, y, z)*
pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 
6*pow(a_BH, 4)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 6)*
fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 6)*((fd_ks_dZ_D0_) )*((fd_ks_dZ_D2_) ) + 
2*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 5)*
fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 5)*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ) - 24*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 4)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 
pow(fd_ks_R(x, y, z), 9)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 2*pow(fd_ks_R(x, y, z), 8)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)/
(pow(a_BH, 6)*pow(fd_ks_Z(x, y, z), 6) + 3*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 4)*
pow(fd_ks_Z(x, y, z), 4) + 3*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 8)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 12));
}
KS_func_def_macro(ddH_D0D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(-2*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 6*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 2)*((fd_ks_dZ_D0_) )*((fd_ks_dZ_D1_) ) + 
3*pow(a_BH, 4)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
6*pow(a_BH, 4)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) - 6*pow(a_BH, 4)*fd_ks_R(x, y, z)*
pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 
6*pow(a_BH, 4)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 6)*
fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 6)*((fd_ks_dZ_D0_) )*((fd_ks_dZ_D1_) ) + 
2*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 5)*
fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 5)*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ) - 24*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 4)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) - 
pow(fd_ks_R(x, y, z), 9)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 2*pow(fd_ks_R(x, y, z), 8)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*fd_ks_R(x, y, z)/
(pow(a_BH, 6)*pow(fd_ks_Z(x, y, z), 6) + 3*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 4)*
pow(fd_ks_Z(x, y, z), 4) + 3*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 8)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 12));
}
KS_func_def_macro(ddH_D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) - 4*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - 2*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) + 8*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
fd_ks_R(x, y, z)/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3);
}
KS_func_def_macro(dddH_D2D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D2D2D2 KS_func_pass_args_macro ) - 6*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) - 6*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 4*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + pow(((fd_ks_dZ_D2_) ), 2))*
pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*pow(fd_ks_R(x, y, z), 4))*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(9*pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_dddZ_D2D2D2 KS_func_pass_args_macro ) + 3*
((fd_ks_dZ_D2_) )*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 3) - pow(a_BH, 2)*(pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dddR_D2D2D2 KS_func_pass_args_macro ) + 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 6*fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 6*(fd_ks_dR_D2 KS_func_pass_args_macro )*
pow(((fd_ks_dZ_D2_) ), 2))*pow(fd_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 3) + 
(fd_ks_R(x, y, z)*(fd_ks_dddR_D2D2D2 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 5)) + 24*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 24*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D2_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) - 48*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ), 3))/
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D2D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D2_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D0 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D0D2D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D2_) )*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 4*
fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*pow(((fd_ks_dZ_D2_) ), 2) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 
(fd_ks_R(x, y, z)*(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D2_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D0_) )*((fd_ks_dZ_D2_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro )) - 
48*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*pow(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D0D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D0 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D0D1D1 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D1_) )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + 2*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 4*
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*pow(((fd_ks_dZ_D1_) ), 2) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ))*pow(fd_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 
(fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D0_) )*((fd_ks_dZ_D1_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 2*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro )) - 
48*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*pow(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D1 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D0D0D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D1 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D0D0D1 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddZ_D0D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
2*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D1_) ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 4*fd_ks_Z(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D0_) )*((fd_ks_dZ_D1_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 2*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )) - 
48*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D1 KS_func_pass_args_macro )))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D0D0D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D0D0D2 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddZ_D0D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
2*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D2_) ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 4*fd_ks_Z(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D0_) )*((fd_ks_dZ_D2_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )) - 
48*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D2 KS_func_pass_args_macro )))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D0D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(-2*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 5)*
(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) + 6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*
pow(fd_ks_Z(x, y, z), 4)*((fd_ks_dZ_D0_) )*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 
6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 4)*
((fd_ks_dZ_D1_) )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 4)*
((fd_ks_dZ_D2_) )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) - 24*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D0_) )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ) + 3*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) - 6*pow(a_BH, 6)*
pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) - 6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) - 
6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) - 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
((fd_ks_dZ_D0_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
((fd_ks_dZ_D1_) )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
((fd_ks_dZ_D2_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 18*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 4)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ) + 18*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ) + 18*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ) + 6*pow(a_BH, 6)*fd_ks_R(x, y, z)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
6*pow(a_BH, 6)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 6*pow(a_BH, 6)*fd_ks_R(x, y, z)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
12*pow(a_BH, 6)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) - 12*pow(a_BH, 6)*
fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) - 12*pow(a_BH, 6)*
fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 6*pow(a_BH, 6)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - 4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*
pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) + 4*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 7)*pow(fd_ks_Z(x, y, z), 2)*((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*
pow(fd_ks_Z(x, y, z), 2)*((fd_ks_dZ_D1_) )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 
4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*pow(fd_ks_Z(x, y, z), 2)*
((fd_ks_dZ_D2_) )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 24*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ) + 5*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + 
4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D0_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D1_) )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D2_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 68*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ) - 68*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ) - 68*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ) - 18*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 5)*
pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
18*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 4)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 18*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 4)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) - 186*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 4)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 11)*fd_ks_Z(x, y, z)*(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) - 2*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 11)*((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 11)*
((fd_ks_dZ_D1_) )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) - 2*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 11)*((fd_ks_dZ_D2_) )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 10)*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ) + 10*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 10)*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ) - 22*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
22*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 22*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 60*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) - 60*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 9)*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) - 60*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 9)*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 186*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 8)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 
pow(fd_ks_R(x, y, z), 14)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + 2*
pow(fd_ks_R(x, y, z), 13)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
2*pow(fd_ks_R(x, y, z), 13)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
2*pow(fd_ks_R(x, y, z), 13)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
6*pow(fd_ks_R(x, y, z), 12)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/(pow(a_BH, 8)*pow(fd_ks_Z(x, y, z), 8) + 4*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 4)*pow(fd_ks_Z(x, y, z), 6) + 6*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 8)*pow(fd_ks_Z(x, y, z), 4) + 4*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 12)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 16));
}
KS_func_def_macro(dddH_D1D1D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D0 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D0D1D1 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D1_) )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + 2*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 4*
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*pow(((fd_ks_dZ_D1_) ), 2) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ))*pow(fd_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + 
(fd_ks_R(x, y, z)*(fd_ks_dddR_D0D1D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D0_) )*((fd_ks_dZ_D1_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 2*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro )) - 
48*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*pow(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D1 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D1D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D2_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D1 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D1D2D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D1_) )*
(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D2_) )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 4*
fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*pow(((fd_ks_dZ_D2_) ), 2) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 
(fd_ks_R(x, y, z)*(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D2_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro )) - 
48*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*pow(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D1D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D1D1D2 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D1_) )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddZ_D1D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
2*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D2_) ) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 4*fd_ks_Z(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro )) - 
48*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D2 KS_func_pass_args_macro )))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D1D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(-2*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 5)*
(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) + 6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*
pow(fd_ks_Z(x, y, z), 4)*((fd_ks_dZ_D0_) )*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 
6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 4)*
((fd_ks_dZ_D1_) )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 4)*
((fd_ks_dZ_D2_) )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) - 24*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D0_) )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ) + 3*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) - 6*pow(a_BH, 6)*
pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) - 6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) - 
6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) - 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
((fd_ks_dZ_D0_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
((fd_ks_dZ_D1_) )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
((fd_ks_dZ_D2_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 18*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 4)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ) + 18*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ) + 18*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ) + 6*pow(a_BH, 6)*fd_ks_R(x, y, z)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
6*pow(a_BH, 6)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 6*pow(a_BH, 6)*fd_ks_R(x, y, z)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
12*pow(a_BH, 6)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) - 12*pow(a_BH, 6)*
fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) - 12*pow(a_BH, 6)*
fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 6*pow(a_BH, 6)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - 4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*
pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) + 4*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 7)*pow(fd_ks_Z(x, y, z), 2)*((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*
pow(fd_ks_Z(x, y, z), 2)*((fd_ks_dZ_D1_) )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 
4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*pow(fd_ks_Z(x, y, z), 2)*
((fd_ks_dZ_D2_) )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 24*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ) + 5*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + 
4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D0_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D1_) )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D2_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 68*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ) - 68*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ) - 68*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ) - 18*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 5)*
pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
18*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 4)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 18*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 4)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) - 186*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 4)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 11)*fd_ks_Z(x, y, z)*(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) - 2*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 11)*((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 11)*
((fd_ks_dZ_D1_) )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) - 2*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 11)*((fd_ks_dZ_D2_) )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 10)*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ) + 10*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 10)*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ) - 22*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
22*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 22*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 60*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) - 60*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 9)*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) - 60*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 9)*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 186*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 8)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 
pow(fd_ks_R(x, y, z), 14)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + 2*
pow(fd_ks_R(x, y, z), 13)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
2*pow(fd_ks_R(x, y, z), 13)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
2*pow(fd_ks_R(x, y, z), 13)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
6*pow(fd_ks_R(x, y, z), 12)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/(pow(a_BH, 8)*pow(fd_ks_Z(x, y, z), 8) + 4*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 4)*pow(fd_ks_Z(x, y, z), 6) + 6*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 8)*pow(fd_ks_Z(x, y, z), 4) + 4*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 12)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 16));
}
KS_func_def_macro(dddH_D0D2D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D2_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D0 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D0D2D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D2_) )*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 4*
fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*pow(((fd_ks_dZ_D2_) ), 2) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 
(fd_ks_R(x, y, z)*(fd_ks_dddR_D0D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D2_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D0_) )*((fd_ks_dZ_D2_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro )) - 
48*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*pow(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D1D1D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D1D1D1 KS_func_pass_args_macro ) - 6*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) - 6*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 4*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + pow(((fd_ks_dZ_D1_) ), 2))*
pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 
pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*pow(fd_ks_R(x, y, z), 4))*
(fd_ks_dR_D1 KS_func_pass_args_macro ) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(9*pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_dddZ_D1D1D1 KS_func_pass_args_macro ) + 3*
((fd_ks_dZ_D1_) )*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 3) - pow(a_BH, 2)*(pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dddR_D1D1D1 KS_func_pass_args_macro ) + 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 6*fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + 6*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 3) + 
(fd_ks_R(x, y, z)*(fd_ks_dddR_D1D1D1 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 5)) + 24*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + 24*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) - 48*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ), 3))/
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D0D2D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D0D0D2 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + (fd_ks_ddZ_D0D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
2*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D2_) ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 4*fd_ks_Z(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D0_) )*((fd_ks_dZ_D2_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )) - 
48*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D2 KS_func_pass_args_macro )))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D0D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(-2*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 5)*
(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) + 6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*
pow(fd_ks_Z(x, y, z), 4)*((fd_ks_dZ_D0_) )*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 
6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 4)*
((fd_ks_dZ_D1_) )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 4)*
((fd_ks_dZ_D2_) )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) - 24*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 3)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D0_) )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ) + 3*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) - 6*pow(a_BH, 6)*
pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) - 6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) - 
6*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) - 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
((fd_ks_dZ_D0_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
((fd_ks_dZ_D1_) )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 6*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 5)*
((fd_ks_dZ_D2_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 18*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*pow(fd_ks_Z(x, y, z), 4)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ) + 18*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ) + 18*pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 2)*
pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ) + 6*pow(a_BH, 6)*fd_ks_R(x, y, z)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
6*pow(a_BH, 6)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 6*pow(a_BH, 6)*fd_ks_R(x, y, z)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
12*pow(a_BH, 6)*fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) - 12*pow(a_BH, 6)*
fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) - 12*pow(a_BH, 6)*
fd_ks_R(x, y, z)*pow(fd_ks_Z(x, y, z), 5)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 6*pow(a_BH, 6)*
pow(fd_ks_Z(x, y, z), 6)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) - 4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*
pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) + 4*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 7)*pow(fd_ks_Z(x, y, z), 2)*((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*
pow(fd_ks_Z(x, y, z), 2)*((fd_ks_dZ_D1_) )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 
4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*pow(fd_ks_Z(x, y, z), 2)*
((fd_ks_dZ_D2_) )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 24*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 7)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ) + 5*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + 
4*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D0_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D1_) )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 4*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 3)*
((fd_ks_dZ_D2_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 68*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ) - 68*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ) - 68*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 6)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ) - 18*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 5)*
pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
18*pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 4)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 18*
pow(a_BH, 4)*pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 4)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 5)*pow(fd_ks_Z(x, y, z), 3)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) - 186*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 4)*pow(fd_ks_Z(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 11)*fd_ks_Z(x, y, z)*(fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ) - 2*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 11)*((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) - 2*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 11)*
((fd_ks_dZ_D1_) )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) - 2*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 11)*((fd_ks_dZ_D2_) )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 10)*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D2 KS_func_pass_args_macro ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) )*
(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 10*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ) + 10*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 10)*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D2_) ) + 10*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 10)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ) - 22*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) - 
22*pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) - 22*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 60*
pow(a_BH, 2)*pow(fd_ks_R(x, y, z), 9)*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) - 60*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 9)*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) - 60*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 9)*fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ) + 186*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 8)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) - 
pow(fd_ks_R(x, y, z), 14)*(fd_ks_dddR_D0D1D2 KS_func_pass_args_macro ) + 2*
pow(fd_ks_R(x, y, z), 13)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
2*pow(fd_ks_R(x, y, z), 13)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D0D2 KS_func_pass_args_macro ) + 
2*pow(fd_ks_R(x, y, z), 13)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) - 
6*pow(fd_ks_R(x, y, z), 12)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))/(pow(a_BH, 8)*pow(fd_ks_Z(x, y, z), 8) + 4*
pow(a_BH, 6)*pow(fd_ks_R(x, y, z), 4)*pow(fd_ks_Z(x, y, z), 6) + 6*pow(a_BH, 4)*
pow(fd_ks_R(x, y, z), 8)*pow(fd_ks_Z(x, y, z), 4) + 4*pow(a_BH, 2)*
pow(fd_ks_R(x, y, z), 12)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 16));
}
KS_func_def_macro(dddH_D0D0D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D0D0D0 KS_func_pass_args_macro ) - 6*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) - 6*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 4*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + pow(((fd_ks_dZ_D0_) ), 2))*
pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 
pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*pow(fd_ks_R(x, y, z), 4))*
(fd_ks_dR_D0 KS_func_pass_args_macro ) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(9*pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_dddZ_D0D0D0 KS_func_pass_args_macro ) + 3*
((fd_ks_dZ_D0_) )*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 3) - pow(a_BH, 2)*(pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dddR_D0D0D0 KS_func_pass_args_macro ) + 6*fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 6*fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ) + 6*(fd_ks_dR_D0 KS_func_pass_args_macro )*
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 3) + 
(fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D0 KS_func_pass_args_macro ) + 3*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 5)) + 24*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + 24*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) - 48*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ), 3))/
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D1D1D2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D2 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D1D1D2 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D1_) )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + (fd_ks_ddZ_D1D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
2*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D2_) ) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 4*fd_ks_Z(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_dddR_D1D1D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(fd_ks_ddR_D1D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D1 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D1_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D1D1 KS_func_pass_args_macro ) + pow((fd_ks_dR_D1 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro )) - 
48*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ), 2)*(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D2 KS_func_pass_args_macro )))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D2D2D1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D2 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D2_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D1 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D1D2D2 KS_func_pass_args_macro ) + ((fd_ks_dZ_D1_) )*
(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D2_) )*
(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 
2*fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro )*
((fd_ks_dZ_D2_) ) + 2*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2)*
((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 4*
fd_ks_Z(x, y, z)*(fd_ks_dR_D2 KS_func_pass_args_macro )*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) )*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*pow(((fd_ks_dZ_D2_) ), 2) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) )*
((fd_ks_dZ_D2_) ))*pow(fd_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + 
(fd_ks_R(x, y, z)*(fd_ks_dddR_D1D2D2 KS_func_pass_args_macro ) + (fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D2 KS_func_pass_args_macro )*
(fd_ks_ddR_D1D2 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D2D2 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D2_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D2D2 KS_func_pass_args_macro ) + pow((fd_ks_dR_D2 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D1D2 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D1_) )*((fd_ks_dZ_D2_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D2_) ) + 2*
(fd_ks_dR_D2 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro )*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D1D2 KS_func_pass_args_macro ) + 
(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D2 KS_func_pass_args_macro )) - 
48*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*pow(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D2_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D2 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D2 KS_func_pass_args_macro ), 2))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}
KS_func_def_macro(dddH_D0D1D0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// R
// Z
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 3)*
pow(fd_ks_R(x, y, z), 2)*(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) - 2*
pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ))*
fd_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + ((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ))*pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 2*(fd_ks_dR_D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
(fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + (fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ))*pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))*(fd_ks_dR_D1 KS_func_pass_args_macro )) - 2*pow(pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2) + pow(fd_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*
(fd_ks_dddZ_D0D0D1 KS_func_pass_args_macro ) + 2*((fd_ks_dZ_D0_) )*
(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + (fd_ks_ddZ_D0D0 KS_func_pass_args_macro )*
((fd_ks_dZ_D1_) ))*pow(fd_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ) + 
2*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*((fd_ks_dZ_D1_) ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro )*
((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 4*fd_ks_Z(x, y, z)*
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 2*fd_ks_Z(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 2*
fd_ks_Z(x, y, z)*(fd_ks_dR_D1 KS_func_pass_args_macro )*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
4*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
4*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) )*
((fd_ks_dZ_D1_) ) + 2*(fd_ks_dR_D1 KS_func_pass_args_macro )*
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_dddR_D0D0D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(fd_ks_ddR_D0D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D1 KS_func_pass_args_macro ) + 
((fd_ks_dZ_D0_) )*((fd_ks_dZ_D1_) ))*
pow(fd_ks_R(x, y, z), 2) - pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
2*(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D1_) ) + 2*
(fd_ks_dR_D1 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro )*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + (fd_ks_R(x, y, z)*(fd_ks_ddR_D0D1 KS_func_pass_args_macro ) + 
(fd_ks_dR_D0 KS_func_pass_args_macro )*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
pow(fd_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(-pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + 4*
(fd_ks_dR_D0 KS_func_pass_args_macro )*((fd_ks_dZ_D0_) ))*fd_ks_R(x, y, z)*
fd_ks_Z(x, y, z) + pow(a_BH, 2)*(fd_ks_Z(x, y, z)*(fd_ks_ddZ_D0D0 KS_func_pass_args_macro ) + 
pow(((fd_ks_dZ_D0_) ), 2))*pow(fd_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2) + (fd_ks_R(x, y, z)*
(fd_ks_ddR_D0D0 KS_func_pass_args_macro ) + pow((fd_ks_dR_D0 KS_func_pass_args_macro ), 2))*
pow(fd_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4))*((pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ))*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + 2*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D1 KS_func_pass_args_macro ))*
(fd_ks_dR_D0 KS_func_pass_args_macro ))*(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*
((fd_ks_dZ_D0_) ) - pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*
(fd_ks_dR_D0 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro )) - 
48*pow(pow(a_BH, 2)*fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D0_) ) - 
pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D0 KS_func_pass_args_macro ) + 
pow(fd_ks_R(x, y, z), 4)*(fd_ks_dR_D0 KS_func_pass_args_macro ), 2)*(pow(a_BH, 2)*
fd_ks_R(x, y, z)*fd_ks_Z(x, y, z)*((fd_ks_dZ_D1_) ) - pow(a_BH, 2)*
pow(fd_ks_Z(x, y, z), 2)*(fd_ks_dR_D1 KS_func_pass_args_macro ) + pow(fd_ks_R(x, y, z), 4)*
(fd_ks_dR_D1 KS_func_pass_args_macro )))/pow(pow(a_BH, 2)*pow(fd_ks_Z(x, y, z), 2) + 
pow(fd_ks_R(x, y, z), 4), 4);
}

