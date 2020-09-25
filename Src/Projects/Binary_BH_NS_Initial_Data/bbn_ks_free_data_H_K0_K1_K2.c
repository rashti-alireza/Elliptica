/* msimplify(...) */
/* msimplify(...) */
/* msimplify(...) */
#include "bbn_ks_free_data_analytic.h"
KS_func_def_macro(K0)KS_func_args_macro;
KS_func_def_macro(dK0_D1)KS_func_args_macro;
KS_func_def_macro(dK0_D0)KS_func_args_macro;
KS_func_def_macro(dK0_D2)KS_func_args_macro;
KS_func_def_macro(ddK0_D1D2)KS_func_args_macro;
KS_func_def_macro(ddK0_D0D0)KS_func_args_macro;
KS_func_def_macro(ddK0_D2D2)KS_func_args_macro;
KS_func_def_macro(ddK0_D1D1)KS_func_args_macro;
KS_func_def_macro(ddK0_D0D2)KS_func_args_macro;
KS_func_def_macro(ddK0_D0D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D2D2D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D1D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D2D2D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D2D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D0D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D1D0)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddK0_D1D1D2)KS_func_args_macro;
KS_func_def_macro(K1)KS_func_args_macro;
KS_func_def_macro(dK1_D2)KS_func_args_macro;
KS_func_def_macro(dK1_D0)KS_func_args_macro;
KS_func_def_macro(dK1_D1)KS_func_args_macro;
KS_func_def_macro(ddK1_D0D1)KS_func_args_macro;
KS_func_def_macro(ddK1_D1D1)KS_func_args_macro;
KS_func_def_macro(ddK1_D2D2)KS_func_args_macro;
KS_func_def_macro(ddK1_D0D2)KS_func_args_macro;
KS_func_def_macro(ddK1_D1D2)KS_func_args_macro;
KS_func_def_macro(ddK1_D0D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D2D2D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D2D2D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D2D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D1D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D1D0)KS_func_args_macro;
KS_func_def_macro(dddK1_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D1D1)KS_func_args_macro;
KS_func_def_macro(dddK1_D0D0D1)KS_func_args_macro;
KS_func_def_macro(K2)KS_func_args_macro;
KS_func_def_macro(dK2_D0)KS_func_args_macro;
KS_func_def_macro(dK2_D1)KS_func_args_macro;
KS_func_def_macro(dK2_D2)KS_func_args_macro;
KS_func_def_macro(ddK2_D1D2)KS_func_args_macro;
KS_func_def_macro(ddK2_D0D1)KS_func_args_macro;
KS_func_def_macro(ddK2_D2D2)KS_func_args_macro;
KS_func_def_macro(ddK2_D1D1)KS_func_args_macro;
KS_func_def_macro(ddK2_D0D0)KS_func_args_macro;
KS_func_def_macro(ddK2_D0D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D2D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D1D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D0D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D1D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D1D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddK2_D2D2D1)KS_func_args_macro;
KS_func_def_macro(dddK2_D2D2D2)KS_func_args_macro;
KS_func_def_macro(H)KS_func_args_macro;
KS_func_def_macro(dH_D1)KS_func_args_macro;
KS_func_def_macro(dH_D0)KS_func_args_macro;
KS_func_def_macro(dH_D2)KS_func_args_macro;
KS_func_def_macro(ddH_D1D1)KS_func_args_macro;
KS_func_def_macro(ddH_D0D2)KS_func_args_macro;
KS_func_def_macro(ddH_D2D2)KS_func_args_macro;
KS_func_def_macro(ddH_D0D1)KS_func_args_macro;
KS_func_def_macro(ddH_D0D0)KS_func_args_macro;
KS_func_def_macro(ddH_D1D2)KS_func_args_macro;
KS_func_def_macro(dddH_D0D2D2)KS_func_args_macro;
KS_func_def_macro(dddH_D1D2D0)KS_func_args_macro;
KS_func_def_macro(dddH_D0D1D0)KS_func_args_macro;
KS_func_def_macro(dddH_D1D1D2)KS_func_args_macro;
KS_func_def_macro(dddH_D0D2D0)KS_func_args_macro;
KS_func_def_macro(dddH_D1D2D1)KS_func_args_macro;
KS_func_def_macro(dddH_D0D0D0)KS_func_args_macro;
KS_func_def_macro(dddH_D2D2D2)KS_func_args_macro;
KS_func_def_macro(dddH_D1D2D2)KS_func_args_macro;
KS_func_def_macro(dddH_D0D1D2)KS_func_args_macro;
KS_func_def_macro(dddH_D0D0D1)KS_func_args_macro;
KS_func_def_macro(dddH_D1D1D1)KS_func_args_macro;
KS_func_def_macro(dddH_D0D2D1)KS_func_args_macro;
KS_func_def_macro(dddH_D1D1D0)KS_func_args_macro;
KS_func_def_macro(dddH_D0D0D2)KS_func_args_macro;
KS_func_def_macro(dddH_D2D2D0)KS_func_args_macro;
KS_func_def_macro(dddH_D0D1D1)KS_func_args_macro;
KS_func_def_macro(dddH_D2D2D1)KS_func_args_macro;
KS_func_def_macro(K0)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// R
// X
// Y
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))/(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2));
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
((pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) - 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2);
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
((pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) - 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2);
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
((pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) - 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(a_BH*Derivative(bbn_ks_Y(x, y, z), y, z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y, z) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), y)) - 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y, z) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 
8*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (x, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (x, 2)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), x)) - 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + 2*(a_BH*
Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 
8*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
pow(Derivative(bbn_ks_R(x, y, z), x), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (z, 2)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), z)) - 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) + 
8*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (y, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (y, 2)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), y)) - 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 2*(a_BH*
Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 
8*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(a_BH*Derivative(bbn_ks_Y(x, y, z), x, z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, z) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x)) - 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 
8*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(a_BH*Derivative(bbn_ks_Y(x, y, z), x, y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, y) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x)) - 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 
8*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), y, (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y, (z, 2)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_X(x, y, z), y) + 2*Derivative(bbn_ks_X(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + (a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), (z, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (z, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), y, z) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), x, (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, (z, 2)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_X(x, y, z), x) + 2*Derivative(bbn_ks_X(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + (a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), (z, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (z, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, z) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), x, y, z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, y, z) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x, y) + Derivative(bbn_ks_X(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_X(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_X(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y, z) + (a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x, y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, y) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y) + 3*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), x, (y, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, (y, 2)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_X(x, y, z), x) + 2*Derivative(bbn_ks_X(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + (a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), (y, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (y, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), x, y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, y) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 2*(a_BH*
Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (z, 3)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (z, 3)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 3)) + 3*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), (z, 2)) + 3*Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_X(x, y, z), z)) - 2*pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 3)) + 3*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
3*(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
3*(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 
3*(a_BH*Derivative(bbn_ks_Y(x, y, z), (z, 2)) + bbn_ks_R(x, y, z)*
Derivative(bbn_ks_X(x, y, z), (z, 2)) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_X(x, y, z), z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z)) + 24*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) - 48*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), z), 3))/
pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), y, (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y, (z, 2)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_X(x, y, z), y) + 2*Derivative(bbn_ks_X(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + (a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), (z, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (z, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), y, z) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), x, (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, (z, 2)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_X(x, y, z), x) + 2*Derivative(bbn_ks_X(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + (a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), (z, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (z, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, z) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (y, 2), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (y, 2), z) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_X(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), y, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 2*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), (y, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (y, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_R(x, y, z), z) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (x, 3)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (x, 3)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 3)) + 3*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), (x, 2)) + 3*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_X(x, y, z), x)) - 2*pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 3)) + 3*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 
3*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 
3*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + 
3*(a_BH*Derivative(bbn_ks_Y(x, y, z), (x, 2)) + bbn_ks_R(x, y, z)*
Derivative(bbn_ks_X(x, y, z), (x, 2)) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x)) + 24*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) - 48*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), x), 3))/
pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (x, 2), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (x, 2), y) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_X(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_X(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 2*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_R(x, y, z), y) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), (x, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (x, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), x, y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), y) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), x, (y, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, (y, 2)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_X(x, y, z), x) + 2*Derivative(bbn_ks_X(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + (a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), (y, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (y, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), x, y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, y) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 2*(a_BH*
Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), x, y, z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, y, z) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x, y) + Derivative(bbn_ks_X(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_X(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_X(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y, z) + (a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x, y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, y) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y) + 3*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (x, 2), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (x, 2), z) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_X(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 2*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), (x, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (x, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), z) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (x, 2), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (x, 2), y) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_X(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_X(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 2*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_R(x, y, z), y) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), (x, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (x, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), x, y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), y) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), x, y, z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, y, z) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x, y) + Derivative(bbn_ks_X(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_X(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_X(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y, z) + (a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x, y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, y) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y) + 3*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (x, 2), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (x, 2), z) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_X(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 2*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), (x, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (x, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), z) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (y, 3)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (y, 3)) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 3)) + 3*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), (y, 2)) + 3*Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_X(x, y, z), y)) - 2*pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 3)) + 3*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
3*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
3*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 
3*(a_BH*Derivative(bbn_ks_Y(x, y, z), (y, 2)) + bbn_ks_R(x, y, z)*
Derivative(bbn_ks_X(x, y, z), (y, 2)) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)) + 24*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + (a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) - 48*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), y), 3))/
pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(a_BH*Derivative(bbn_ks_Y(x, y, z), (y, 2), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (y, 2), z) + bbn_ks_X(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_X(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), y, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 2*(a_BH*bbn_ks_Y(x, y, z) + 
bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 
(a_BH*Derivative(bbn_ks_Y(x, y, z), (y, 2)) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), (y, 2)) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) + 2*
(a_BH*Derivative(bbn_ks_Y(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y, z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_X(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_X(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(2*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*
bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_R(x, y, z), z) + 3*(a_BH*
bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_X(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
bbn_ks_R(x, y, z) - 48*(a_BH*bbn_ks_Y(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_X(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(K1)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// R
// X
// Y
(-a_BH*bbn_ks_X(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))/(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2));
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
((pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) + 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2);
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
((pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2);
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
((pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(-a_BH*Derivative(bbn_ks_X(x, y, z), x, y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, y) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x)) + 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y) + (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) - 
(-a_BH*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) - 
8*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (y, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (y, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), y)) + 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), y), 2) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 
8*(-a_BH*bbn_ks_X(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (z, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), z)) + 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), z), 2) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) + 
8*(-a_BH*bbn_ks_X(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(-a_BH*Derivative(bbn_ks_X(x, y, z), x, z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, z) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x)) + 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) - 
(-a_BH*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) - 
8*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(-a_BH*Derivative(bbn_ks_X(x, y, z), y, z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y, z) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), y)) + 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y, z) + (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) - 
(-a_BH*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) - 
8*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (x, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (x, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), x)) + 2*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*
((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 
8*(-a_BH*bbn_ks_X(x, y, z) + bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 2)*
pow(Derivative(bbn_ks_R(x, y, z), x), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 3);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), x, y, z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, y, z) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x, y) + Derivative(bbn_ks_Y(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_Y(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_Y(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y, z) + (a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) - 
(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) - 
(-a_BH*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_X(x, y, z), x, y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, y) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - (a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y) - 3*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (-
a_BH*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), x, (y, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, (y, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_Y(x, y, z), x) + 2*Derivative(bbn_ks_Y(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + (a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y) - (-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), x) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
pow(Derivative(bbn_ks_R(x, y, z), y), 2) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (y, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (y, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), x, y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) - 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, y) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (z, 3)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (z, 3)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 3)) + 3*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), (z, 2)) + 3*Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_Y(x, y, z), z)) + 2*pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*
((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 3)) + 3*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) - 
3*(-a_BH*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
3*(a_BH*Derivative(bbn_ks_X(x, y, z), z) - bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) - 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*pow(Derivative(bbn_ks_R(x, y, z), z), 2) - 
3*(-a_BH*Derivative(bbn_ks_X(x, y, z), (z, 2)) + bbn_ks_R(x, y, z)*
Derivative(bbn_ks_Y(x, y, z), (z, 2)) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Y(x, y, z), z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z)) + 24*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) - (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) + 48*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), z), 3))/
pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), x, y, z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, y, z) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x, y) + Derivative(bbn_ks_Y(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_Y(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_Y(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y, z) + (a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) - 
(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) - 
(-a_BH*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_X(x, y, z), x, y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, y) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - (a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y) - 3*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (-
a_BH*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (y, 3)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (y, 3)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 3)) + 3*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), (y, 2)) + 3*Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_Y(x, y, z), y)) + 2*pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*
((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 3)) + 3*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (y, 2)) - 
3*(-a_BH*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
3*(a_BH*Derivative(bbn_ks_X(x, y, z), y) - bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) - 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*pow(Derivative(bbn_ks_R(x, y, z), y), 2) - 
3*(-a_BH*Derivative(bbn_ks_X(x, y, z), (y, 2)) + bbn_ks_R(x, y, z)*
Derivative(bbn_ks_Y(x, y, z), (y, 2)) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)) + 24*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) - (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (-a_BH*
Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) + 48*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), y), 3))/
pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), x, y, z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, y, z) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x, y) + Derivative(bbn_ks_Y(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_Y(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_Y(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y, z) + (a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) - 
(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) - 
(-a_BH*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_X(x, y, z), x, y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, y) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) - (-
a_BH*Derivative(bbn_ks_X(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - (a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y) - 3*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (-
a_BH*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), x, (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, (z, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_Y(x, y, z), x) + 2*Derivative(bbn_ks_Y(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + (a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), x) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (z, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) - 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, z) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (x, 3)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (x, 3)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 3)) + 3*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), (x, 2)) + 3*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_Y(x, y, z), x)) + 2*pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 2)*
((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 3)) + 3*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (x, 2)) - 
3*(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 
3*(a_BH*Derivative(bbn_ks_X(x, y, z), x) - bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) - 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*pow(Derivative(bbn_ks_R(x, y, z), x), 2) - 
3*(-a_BH*Derivative(bbn_ks_X(x, y, z), (x, 2)) + bbn_ks_R(x, y, z)*
Derivative(bbn_ks_Y(x, y, z), (x, 2)) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x)) + 24*(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) - (a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (-a_BH*
Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) + 48*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), x), 3))/
pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), x, (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, (z, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_Y(x, y, z), x) + 2*Derivative(bbn_ks_Y(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + (a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), x) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (z, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) - 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, z) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), y, (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y, (z, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_Y(x, y, z), y) + 2*Derivative(bbn_ks_Y(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + (a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), y) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (z, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) - 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), y, z) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (x, 2), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (x, 2), z) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_Y(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 2*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_R(x, y, z), z) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), z) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
pow(Derivative(bbn_ks_R(x, y, z), x), 2) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (x, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (x, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) - 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), z) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z) + 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + (-
a_BH*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (x, 2), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (x, 2), z) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_Y(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 2*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_R(x, y, z), z) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), z) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
pow(Derivative(bbn_ks_R(x, y, z), x), 2) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (x, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (x, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), x, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) - 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), z) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z) + 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + (-
a_BH*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (y, 2), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (y, 2), z) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_Y(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), y, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 2*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_R(x, y, z), z) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), z) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
pow(Derivative(bbn_ks_R(x, y, z), y), 2) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (y, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (y, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) - 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_R(x, y, z), z) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z) + 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (-
a_BH*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (y, 2), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (y, 2), z) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_Y(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), y, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 2*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_R(x, y, z), z) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), z) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
pow(Derivative(bbn_ks_R(x, y, z), y), 2) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (y, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (y, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) - 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_R(x, y, z), z) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z) + 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + (-
a_BH*Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (x, 2), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (x, 2), y) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_Y(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_Y(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 2*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_R(x, y, z), y) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), y) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
pow(Derivative(bbn_ks_R(x, y, z), x), 2) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (x, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (x, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), x, y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) - 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), y) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y) + 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + (-
a_BH*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), y, (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y, (z, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_Y(x, y, z), y) + 2*Derivative(bbn_ks_Y(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + (a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), y) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), z) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (z, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (z, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), y, z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y, z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) - 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), y, z) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), z) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), z) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), x, (y, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, (y, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_Y(x, y, z), x) + 2*Derivative(bbn_ks_Y(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + (a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y) - (-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), x) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
pow(Derivative(bbn_ks_R(x, y, z), y), 2) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (y, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (y, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), x, y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) - 
2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, y) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(bbn_ks_R(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 4);
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
(pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 3)*(-a_BH*Derivative(bbn_ks_X(x, y, z), (x, 2), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (x, 2), y) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_Y(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_Y(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_R(x, y, z), 2), 2)*((a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 2*(a_BH*bbn_ks_X(x, y, z) - 
bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_R(x, y, z), y) - 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) - 2*(-a_BH*
Derivative(bbn_ks_X(x, y, z), x) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_X(x, y, z), y) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + (a_BH*Derivative(bbn_ks_X(x, y, z), y) - 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) - bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*
pow(Derivative(bbn_ks_R(x, y, z), x), 2) - (-a_BH*Derivative(bbn_ks_X(x, y, z), (x, 2)) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), (x, 2)) + bbn_ks_Y(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) - 2*(-
a_BH*Derivative(bbn_ks_X(x, y, z), x, y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x, y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Y(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Y(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*
bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) - 
(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), y) - 3*(a_BH*
bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y) + 2*(-a_BH*Derivative(bbn_ks_X(x, y, z), x) + 
bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), x) + bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + (-
a_BH*Derivative(bbn_ks_X(x, y, z), y) + bbn_ks_R(x, y, z)*Derivative(bbn_ks_Y(x, y, z), y) + 
bbn_ks_Y(x, y, z)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
bbn_ks_R(x, y, z) + 48*(a_BH*bbn_ks_X(x, y, z) - bbn_ks_R(x, y, z)*bbn_ks_Y(x, y, z))*
pow(bbn_ks_R(x, y, z), 3)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y))/pow(pow(a_BH, 2) + pow(bbn_ks_R(x, y, z), 2), 4);
}
KS_func_def_macro(K2)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// R
// Z
bbn_ks_Z(x, y, z)/bbn_ks_R(x, y, z);
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
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x))/pow(bbn_ks_R(x, y, z), 2);
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
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y))/pow(bbn_ks_R(x, y, z), 2);
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
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), z))/pow(bbn_ks_R(x, y, z), 2);
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
(-(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z) + pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_Z(x, y, z), y, z) + 2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))/pow(bbn_ks_R(x, y, z), 3);
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
(-(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) + pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_Z(x, y, z), x, y) + 2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y))/pow(bbn_ks_R(x, y, z), 3);
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
(-(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z) + pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 2*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), z), 2))/
pow(bbn_ks_R(x, y, z), 3);
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
(-(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z) + pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 2*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 2))/
pow(bbn_ks_R(x, y, z), 3);
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
(-(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) + pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 2*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2))/
pow(bbn_ks_R(x, y, z), 3);
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
(-(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) + pow(bbn_ks_R(x, y, z), 2)*
Derivative(bbn_ks_Z(x, y, z), x, z) + 2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z))/pow(bbn_ks_R(x, y, z), 3);
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
(2*(2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), z) + 
pow(Derivative(bbn_ks_R(x, y, z), x), 2)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_Z(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, z))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), (x, 2), z) - 6*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(bbn_ks_R(x, y, z), 4);
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
(2*(2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_R(x, y, z), z) + 
pow(Derivative(bbn_ks_R(x, y, z), y), 2)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_Z(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), y, z))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), (y, 2), z) - 6*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(bbn_ks_R(x, y, z), 4);
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
(2*(2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), y) + 
pow(Derivative(bbn_ks_R(x, y, z), x), 2)*Derivative(bbn_ks_Z(x, y, z), y) + 2*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_Z(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, y))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), (x, 2), y) - 6*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y))/pow(bbn_ks_R(x, y, z), 4);
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
(2*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z) + pow(Derivative(bbn_ks_R(x, y, z), z), 2)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_Z(x, y, z), x) + 2*Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x, (z, 2)) - 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(bbn_ks_R(x, y, z), 4);
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
(2*(2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), y) + 
pow(Derivative(bbn_ks_R(x, y, z), x), 2)*Derivative(bbn_ks_Z(x, y, z), y) + 2*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_Z(x, y, z), y) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, y))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), (x, 2), y) - 6*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y))/pow(bbn_ks_R(x, y, z), 4);
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
(2*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, y) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y) + pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_Z(x, y, z), x) + 2*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x, (y, 2)) - 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2))/pow(bbn_ks_R(x, y, z), 4);
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
(2*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
2*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z) + pow(Derivative(bbn_ks_R(x, y, z), z), 2)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_Z(x, y, z), y) + 2*Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), y, (z, 2)) - 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(bbn_ks_R(x, y, z), 4);
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
(6*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) - 
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 3)) + 3*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 3*Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_Z(x, y, z), y))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), (y, 3)) - 6*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 3))/
pow(bbn_ks_R(x, y, z), 4);
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
(2*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) + 
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x, y, z) - 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z))/
pow(bbn_ks_R(x, y, z), 4);
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
(2*(2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_R(x, y, z), z) + 
pow(Derivative(bbn_ks_R(x, y, z), y), 2)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_Z(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), y, z))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), (y, 2), z) - 6*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(bbn_ks_R(x, y, z), 4);
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
(2*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) + 
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x, y, z) - 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z))/
pow(bbn_ks_R(x, y, z), 4);
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
(2*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, y) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y) + pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_Z(x, y, z), x) + 2*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x, (y, 2)) - 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2))/pow(bbn_ks_R(x, y, z), 4);
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
(2*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z) + pow(Derivative(bbn_ks_R(x, y, z), z), 2)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_Z(x, y, z), x) + 2*Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x, (z, 2)) - 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(bbn_ks_R(x, y, z), 4);
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
(6*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) - 
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 3)) + 3*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 3*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_Z(x, y, z), x))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), (x, 3)) - 6*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 3))/
pow(bbn_ks_R(x, y, z), 4);
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
(2*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) + 
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x, y, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, y))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x, y, z) - 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z))/
pow(bbn_ks_R(x, y, z), 4);
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
(2*(2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), z) + 
pow(Derivative(bbn_ks_R(x, y, z), x), 2)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_Z(x, y, z), z) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), x, z))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), (x, 2), z) - 6*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z))/pow(bbn_ks_R(x, y, z), 4);
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
(2*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
2*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z) + pow(Derivative(bbn_ks_R(x, y, z), z), 2)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z) - (bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_Z(x, y, z), y) + 2*Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), y, (z, 2)) - 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2))/pow(bbn_ks_R(x, y, z), 4);
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
(6*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) - 
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 3)) + 3*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 3*Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), (z, 3)) - 6*bbn_ks_Z(x, y, z)*pow(Derivative(bbn_ks_R(x, y, z), z), 3))/
pow(bbn_ks_R(x, y, z), 4);
}
KS_func_def_macro(H)KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// R
// Z
Lambda*M_BH*pow(bbn_ks_R(x, y, z), 3)/(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4));
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
-Lambda*M_BH*(2*pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
pow(bbn_ks_R(x, y, z), 2)/(pow(a_BH, 4)*pow(bbn_ks_Z(x, y, z), 4) + 2*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 4)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 8));
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
-Lambda*M_BH*(2*pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
pow(bbn_ks_R(x, y, z), 2)/(pow(a_BH, 4)*pow(bbn_ks_Z(x, y, z), 4) + 2*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 4)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 8));
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
-Lambda*M_BH*(2*pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 2)/(pow(a_BH, 4)*pow(bbn_ks_Z(x, y, z), 4) + 2*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 4)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 8));
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) - 4*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), y) - 2*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
4*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
pow(bbn_ks_R(x, y, z), 4)) + 8*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y), 2))*
bbn_ks_R(x, y, z)/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3);
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
Lambda*M_BH*(-2*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x, z) + 6*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), z) + 
3*pow(a_BH, 4)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x, z) - 
6*pow(a_BH, 4)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) - 6*pow(a_BH, 4)*bbn_ks_R(x, y, z)*
pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x) + 
6*pow(a_BH, 4)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) - 2*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 6)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, z) - 2*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 6)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), z) + 
2*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x, z) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 5)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), z) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 5)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x) - 24*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 4)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) - 
pow(bbn_ks_R(x, y, z), 9)*Derivative(bbn_ks_R(x, y, z), x, z) + 2*pow(bbn_ks_R(x, y, z), 8)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)/
(pow(a_BH, 6)*pow(bbn_ks_Z(x, y, z), 6) + 3*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 4)*
pow(bbn_ks_Z(x, y, z), 4) + 3*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 8)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 12));
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) - 4*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), z) - 2*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
4*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + pow(Derivative(bbn_ks_R(x, y, z), z), 2))*
pow(bbn_ks_R(x, y, z), 4)) + 8*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z), 2))*
bbn_ks_R(x, y, z)/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3);
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
Lambda*M_BH*(-2*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x, y) + 6*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y) + 
3*pow(a_BH, 4)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x, y) - 
6*pow(a_BH, 4)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) - 6*pow(a_BH, 4)*bbn_ks_R(x, y, z)*
pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x) + 
6*pow(a_BH, 4)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) - 2*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 6)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, y) - 2*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 6)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y) + 
2*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x, y) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 5)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 5)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x) - 24*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 4)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) - 
pow(bbn_ks_R(x, y, z), 9)*Derivative(bbn_ks_R(x, y, z), x, y) + 2*pow(bbn_ks_R(x, y, z), 8)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))*bbn_ks_R(x, y, z)/
(pow(a_BH, 6)*pow(bbn_ks_Z(x, y, z), 6) + 3*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 4)*
pow(bbn_ks_Z(x, y, z), 4) + 3*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 8)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 12));
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) - 4*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), x) - 2*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 
4*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
pow(bbn_ks_R(x, y, z), 4)) + 8*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x), 2))*
bbn_ks_R(x, y, z)/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3);
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
Lambda*M_BH*(-2*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), y, z) + 6*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 
3*pow(a_BH, 4)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y, z) - 
6*pow(a_BH, 4)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) - 6*pow(a_BH, 4)*bbn_ks_R(x, y, z)*
pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) + 
6*pow(a_BH, 4)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) - 2*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 6)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y, z) - 2*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 6)*Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 
2*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y, z) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 5)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 5)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y) - 24*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 4)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) - 
pow(bbn_ks_R(x, y, z), 9)*Derivative(bbn_ks_R(x, y, z), y, z) + 2*pow(bbn_ks_R(x, y, z), 8)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z))*bbn_ks_R(x, y, z)/
(pow(a_BH, 6)*pow(bbn_ks_Z(x, y, z), 6) + 3*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 4)*
pow(bbn_ks_Z(x, y, z), 4) + 3*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 8)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 12));
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x, (z, 2)) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), x, z))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), z) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + pow(Derivative(bbn_ks_R(x, y, z), z), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), x)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, (z, 2)) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, z))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
4*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*pow(Derivative(bbn_ks_R(x, y, z), z), 2)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 4*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, z) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2))*Derivative(bbn_ks_Z(x, y, z), x) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
2*Derivative(bbn_ks_R(x, y, z), x)*pow(Derivative(bbn_ks_Z(x, y, z), z), 2) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z))*pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + pow(Derivative(bbn_ks_R(x, y, z), z), 2))*
pow(bbn_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, z) + 
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), z) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), x))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z)) - 
48*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*pow(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(-2*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), x, y, z) + 6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y, z) + 
6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) - 24*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 3*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), x, y, z) - 6*pow(a_BH, 6)*
pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) - 6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) - 
6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) - 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 18*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 18*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 18*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) + 6*pow(a_BH, 6)*bbn_ks_R(x, y, z)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
6*pow(a_BH, 6)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + 6*pow(a_BH, 6)*bbn_ks_R(x, y, z)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 
12*pow(a_BH, 6)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) - 12*pow(a_BH, 6)*
bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) - 12*pow(a_BH, 6)*
bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x) + 6*pow(a_BH, 6)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) - 4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*
pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_Z(x, y, z), x, y, z) + 4*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 7)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) + 4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 
4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) + 24*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 5*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x, y, z) + 
4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 68*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) - 68*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) - 68*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) - 18*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 5)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 
18*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - 18*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) + pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x) - 186*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 4)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) - 2*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 11)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, y, z) - 2*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 11)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) - 2*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 11)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) - 2*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 11)*Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x, y, z) + 10*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 10)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 10*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 10)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) - 22*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 
22*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - 22*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 60*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) - 60*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 9)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) - 60*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 9)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x) + 186*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 8)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) - 
pow(bbn_ks_R(x, y, z), 14)*Derivative(bbn_ks_R(x, y, z), x, y, z) + 2*
pow(bbn_ks_R(x, y, z), 13)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
2*pow(bbn_ks_R(x, y, z), 13)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + 
2*pow(bbn_ks_R(x, y, z), 13)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 
6*pow(bbn_ks_R(x, y, z), 12)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))/(pow(a_BH, 8)*pow(bbn_ks_Z(x, y, z), 8) + 4*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 4)*pow(bbn_ks_Z(x, y, z), 6) + 6*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 8)*pow(bbn_ks_Z(x, y, z), 4) + 4*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 12)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 16));
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (x, 2), y) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), x, y) + 
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), (x, 2)))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), x) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), y)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_Z(x, y, z), (x, 2))*
Derivative(bbn_ks_Z(x, y, z), y))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), y) + 
2*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*Derivative(bbn_ks_Z(x, y, z), y) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 4*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x, y) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_Z(x, y, z), y) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
4*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) + 2*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), y))*
pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, y) + 
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y) + 2*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))*
pow(bbn_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), y) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)) - 
48*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x), 2)*(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), y)))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (y, 2), z) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), y, z) + 
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), (y, 2)))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), y) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), z)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_Z(x, y, z), (y, 2))*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_R(x, y, z), z) + 
2*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*Derivative(bbn_ks_Z(x, y, z), z) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 4*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y, z) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_Z(x, y, z), z) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
4*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 
2*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y, z) + 
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
2*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), z) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), y))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y)) - 
48*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y), 2)*(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), z)))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (x, 2), z) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), x, z) + 
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), (x, 2)))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), x) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), z)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_Z(x, y, z), (x, 2))*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), z) + 
2*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*Derivative(bbn_ks_Z(x, y, z), z) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 4*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x, z) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_Z(x, y, z), z) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
4*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, z) + 
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), z) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), x))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)) - 
48*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x), 2)*(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), z)))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (y, 2), z) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), y, z) + 
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), (y, 2)))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), y) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), z)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_Z(x, y, z), (y, 2))*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_R(x, y, z), z) + 
2*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*Derivative(bbn_ks_Z(x, y, z), z) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 4*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y, z) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_Z(x, y, z), z) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
4*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2), z) + 
2*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), y, z) + 
Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y, z) + 
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
2*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), z) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), y))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y)) - 
48*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y), 2)*(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), z)))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (x, 3)) - 6*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) - 6*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + 4*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*
pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 
pow(Derivative(bbn_ks_R(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 4))*
Derivative(bbn_ks_R(x, y, z), x) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(9*pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 3)) + 3*
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), (x, 2)))*
pow(bbn_ks_R(x, y, z), 3) - pow(a_BH, 2)*(pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), (x, 3)) + 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2))*
Derivative(bbn_ks_Z(x, y, z), x) + 6*Derivative(bbn_ks_R(x, y, z), x)*
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 3) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 3)) + 3*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), (x, 2)))*pow(bbn_ks_R(x, y, z), 5)) + 24*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), x) + 24*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
pow(bbn_ks_R(x, y, z), 4)) - 48*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x), 3))/
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (z, 3)) - 6*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) - 6*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 4*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*
pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
pow(Derivative(bbn_ks_R(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 4))*
Derivative(bbn_ks_R(x, y, z), z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(9*pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z) + 
pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 3)) + 3*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)))*
pow(bbn_ks_R(x, y, z), 3) - pow(a_BH, 2)*(pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), (z, 3)) + 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2))*
Derivative(bbn_ks_Z(x, y, z), z) + 6*Derivative(bbn_ks_R(x, y, z), z)*
pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), z), 3) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 3)) + 3*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)))*pow(bbn_ks_R(x, y, z), 5)) + 24*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + 24*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + pow(Derivative(bbn_ks_R(x, y, z), z), 2))*
pow(bbn_ks_R(x, y, z), 4)) - 48*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z), 3))/
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y, (z, 2)) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), y, z))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), z) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + pow(Derivative(bbn_ks_R(x, y, z), z), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), y)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y, (z, 2)) + Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y, z))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
4*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*pow(Derivative(bbn_ks_R(x, y, z), z), 2)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 4*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y, z) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2))*Derivative(bbn_ks_Z(x, y, z), y) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
2*Derivative(bbn_ks_R(x, y, z), y)*pow(Derivative(bbn_ks_Z(x, y, z), z), 2) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z))*pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + pow(Derivative(bbn_ks_R(x, y, z), z), 2))*
pow(bbn_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y, z) + 
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
2*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), z) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), y))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z)) - 
48*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*pow(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(-2*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), x, y, z) + 6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y, z) + 
6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) - 24*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 3*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), x, y, z) - 6*pow(a_BH, 6)*
pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) - 6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) - 
6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) - 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 18*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 18*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 18*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) + 6*pow(a_BH, 6)*bbn_ks_R(x, y, z)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
6*pow(a_BH, 6)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + 6*pow(a_BH, 6)*bbn_ks_R(x, y, z)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 
12*pow(a_BH, 6)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) - 12*pow(a_BH, 6)*
bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) - 12*pow(a_BH, 6)*
bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x) + 6*pow(a_BH, 6)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) - 4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*
pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_Z(x, y, z), x, y, z) + 4*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 7)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) + 4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 
4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) + 24*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 5*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x, y, z) + 
4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 68*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) - 68*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) - 68*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) - 18*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 5)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 
18*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - 18*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) + pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x) - 186*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 4)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) - 2*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 11)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, y, z) - 2*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 11)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) - 2*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 11)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) - 2*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 11)*Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x, y, z) + 10*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 10)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 10*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 10)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) - 22*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 
22*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - 22*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 60*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) - 60*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 9)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) - 60*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 9)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x) + 186*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 8)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) - 
pow(bbn_ks_R(x, y, z), 14)*Derivative(bbn_ks_R(x, y, z), x, y, z) + 2*
pow(bbn_ks_R(x, y, z), 13)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
2*pow(bbn_ks_R(x, y, z), 13)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + 
2*pow(bbn_ks_R(x, y, z), 13)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 
6*pow(bbn_ks_R(x, y, z), 12)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))/(pow(a_BH, 8)*pow(bbn_ks_Z(x, y, z), 8) + 4*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 4)*pow(bbn_ks_Z(x, y, z), 6) + 6*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 8)*pow(bbn_ks_Z(x, y, z), 4) + 4*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 12)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 16));
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (x, 2), y) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), x, y) + 
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), (x, 2)))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), x) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), y)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_Z(x, y, z), (x, 2))*
Derivative(bbn_ks_Z(x, y, z), y))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), y) + 
2*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*Derivative(bbn_ks_Z(x, y, z), y) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 4*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x, y) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_Z(x, y, z), y) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
4*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) + 2*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), y) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), y) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, y) + 
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), y))*
pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, y) + 
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y) + 2*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))*
pow(bbn_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), y) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)) - 
48*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x), 2)*(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), y)))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (y, 3)) - 6*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) - 6*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + 4*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*
pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
pow(Derivative(bbn_ks_R(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 4))*
Derivative(bbn_ks_R(x, y, z), y) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(9*pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 3)) + 3*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), (y, 2)))*
pow(bbn_ks_R(x, y, z), 3) - pow(a_BH, 2)*(pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), (y, 3)) + 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 6*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2))*
Derivative(bbn_ks_Z(x, y, z), y) + 6*Derivative(bbn_ks_R(x, y, z), y)*
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 3) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 3)) + 3*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), (y, 2)))*pow(bbn_ks_R(x, y, z), 5)) + 24*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_R(x, y, z), y) + 24*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
pow(bbn_ks_R(x, y, z), 4)) - 48*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y), 3))/
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(-2*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), x, y, z) + 6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y, z) + 
6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) - 24*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 3)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 3*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), x, y, z) - 6*pow(a_BH, 6)*
pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) - 6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) - 
6*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) - 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - 6*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 5)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 18*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 18*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 18*pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 2)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) + 6*pow(a_BH, 6)*bbn_ks_R(x, y, z)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
6*pow(a_BH, 6)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + 6*pow(a_BH, 6)*bbn_ks_R(x, y, z)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 
12*pow(a_BH, 6)*bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) - 12*pow(a_BH, 6)*
bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) - 12*pow(a_BH, 6)*
bbn_ks_R(x, y, z)*pow(bbn_ks_Z(x, y, z), 5)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x) + 6*pow(a_BH, 6)*
pow(bbn_ks_Z(x, y, z), 6)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) - 4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*
pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_Z(x, y, z), x, y, z) + 4*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 7)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) + 4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 
4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) + 24*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 7)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 5*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x, y, z) + 
4*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, y) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + 4*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 3)*
Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 68*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) - 68*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) - 68*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 6)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) - 18*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 5)*
pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 
18*pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - 18*
pow(a_BH, 4)*pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) + pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 5)*pow(bbn_ks_Z(x, y, z), 3)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x) - 186*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 4)*pow(bbn_ks_Z(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) - 2*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 11)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, y, z) - 2*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 11)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) - 2*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 11)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) - 2*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 11)*Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x, y, z) + 10*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 10)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y, z) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, z) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, z) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) + 10*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 10*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 10)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 10*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 10)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) - 22*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) - 
22*pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) - 22*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 60*
pow(a_BH, 2)*pow(bbn_ks_R(x, y, z), 9)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) - 60*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 9)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y) - 60*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 9)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x) + 186*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 8)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) - 
pow(bbn_ks_R(x, y, z), 14)*Derivative(bbn_ks_R(x, y, z), x, y, z) + 2*
pow(bbn_ks_R(x, y, z), 13)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y, z) + 
2*pow(bbn_ks_R(x, y, z), 13)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, z) + 
2*pow(bbn_ks_R(x, y, z), 13)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, y) - 
6*pow(bbn_ks_R(x, y, z), 12)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))/(pow(a_BH, 8)*pow(bbn_ks_Z(x, y, z), 8) + 4*
pow(a_BH, 6)*pow(bbn_ks_R(x, y, z), 4)*pow(bbn_ks_Z(x, y, z), 6) + 6*pow(a_BH, 4)*
pow(bbn_ks_R(x, y, z), 8)*pow(bbn_ks_Z(x, y, z), 4) + 4*pow(a_BH, 2)*
pow(bbn_ks_R(x, y, z), 12)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 16));
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x, (y, 2)) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), x, y))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), y) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), x)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, (y, 2)) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x, y))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, y) + 
4*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y) + 2*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 4*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, y) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_Z(x, y, z), x) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, y) + 
2*Derivative(bbn_ks_R(x, y, z), x)*pow(Derivative(bbn_ks_Z(x, y, z), y), 2) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y))*pow(bbn_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y))*pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
pow(bbn_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, y) + 
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y) + 2*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), y) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y)) - 
48*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*pow(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), y), 2))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (x, 2), z) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(2*
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), x, z) + 
(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), (x, 2)))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), x) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), z)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_Z(x, y, z), (x, 2))*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*(2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), z) + 
2*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*Derivative(bbn_ks_Z(x, y, z), z) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 4*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x, z) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_Z(x, y, z), z) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
4*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) - 12*
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2), z) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), x, z) + 
Derivative(bbn_ks_R(x, y, z), (x, 2))*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, z) + 
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 4)) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (x, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (x, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), x), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), x), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (x, 2)) + pow(Derivative(bbn_ks_R(x, y, z), x), 2))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), z) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), x))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x)) - 
48*pow(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x), 2)*(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), z)))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x, (z, 2)) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), x, z))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, z) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), z) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + pow(Derivative(bbn_ks_R(x, y, z), z), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), x)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, (z, 2)) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), x, z))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
4*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*pow(Derivative(bbn_ks_R(x, y, z), z), 2)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 4*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x, z) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2))*Derivative(bbn_ks_Z(x, y, z), x) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
2*Derivative(bbn_ks_R(x, y, z), x)*pow(Derivative(bbn_ks_Z(x, y, z), z), 2) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), x, z))*pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + pow(Derivative(bbn_ks_R(x, y, z), z), 2))*
pow(bbn_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, z) + 
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, z) + 
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), z) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), x))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z)) - 
48*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*pow(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x, (y, 2)) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), x, y))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, y) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 2*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), y) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), x)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x, (y, 2)) + Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), x, y))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, y) + 
4*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), y) + 2*pow(Derivative(bbn_ks_R(x, y, z), y), 2)*
Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 4*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x, y) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2))*Derivative(bbn_ks_Z(x, y, z), x) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), x, y) + 
2*Derivative(bbn_ks_R(x, y, z), x)*pow(Derivative(bbn_ks_Z(x, y, z), y), 2) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x)*
Derivative(bbn_ks_Z(x, y, z), y))*pow(bbn_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), x, y))*pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (y, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (y, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), y), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), y), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (y, 2)) + pow(Derivative(bbn_ks_R(x, y, z), y), 2))*
pow(bbn_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x, y) + 
Derivative(bbn_ks_Z(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
2*Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_Z(x, y, z), y) + 2*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), x))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x)*
Derivative(bbn_ks_R(x, y, z), y) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), x, y) + 
Derivative(bbn_ks_R(x, y, z), x)*Derivative(bbn_ks_R(x, y, z), y))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), x) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), x) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*
Derivative(bbn_ks_R(x, y, z), y) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), x))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y)) - 
48*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), x) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), x) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), x))*pow(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), y), 2))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
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
Lambda*M_BH*(pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 3)*
pow(bbn_ks_R(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y, (z, 2)) - 2*
pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*
((pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*Derivative(bbn_ks_R(x, y, z), y, z))*
bbn_ks_R(x, y, z) - 2*pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 2)*(2*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y, z) + Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 2*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z))*pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), z) + 
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + pow(Derivative(bbn_ks_R(x, y, z), z), 2))*
pow(bbn_ks_R(x, y, z), 4))*Derivative(bbn_ks_R(x, y, z), y)) - 2*pow(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4), 2)*(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y, (z, 2)) + Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_Z(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), y, z))*pow(bbn_ks_R(x, y, z), 3) + 3*pow(a_BH, 2)*
(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 
2*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
4*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_Z(x, y, z), z) + 2*pow(Derivative(bbn_ks_R(x, y, z), z), 2)*
Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z) - pow(a_BH, 2)*
(pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + 2*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 4*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y, z) + 2*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2))*Derivative(bbn_ks_Z(x, y, z), y) + 
4*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
2*Derivative(bbn_ks_R(x, y, z), y)*pow(Derivative(bbn_ks_Z(x, y, z), z), 2) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y)*
Derivative(bbn_ks_Z(x, y, z), z))*pow(bbn_ks_R(x, y, z), 2) - 12*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + 
(bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_R(x, y, z), z)*
Derivative(bbn_ks_R(x, y, z), y, z))*pow(bbn_ks_R(x, y, z), 5)) + 8*(pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2) + pow(bbn_ks_R(x, y, z), 4))*((pow(a_BH, 2)*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
(-pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), (z, 2)) + 4*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), z))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), (z, 2)) + 
pow(Derivative(bbn_ks_Z(x, y, z), z), 2))*pow(bbn_ks_R(x, y, z), 2) + 3*pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*pow(Derivative(bbn_ks_R(x, y, z), z), 2) + (bbn_ks_R(x, y, z)*
Derivative(bbn_ks_R(x, y, z), (z, 2)) + pow(Derivative(bbn_ks_R(x, y, z), z), 2))*
pow(bbn_ks_R(x, y, z), 4)) + 2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
(pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y, z) + 
Derivative(bbn_ks_Z(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 2) - pow(a_BH, 2)*(bbn_ks_Z(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
2*Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_Z(x, y, z), z) + 2*
Derivative(bbn_ks_R(x, y, z), z)*Derivative(bbn_ks_Z(x, y, z), y))*bbn_ks_R(x, y, z)*
bbn_ks_Z(x, y, z) + 3*pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y)*
Derivative(bbn_ks_R(x, y, z), z) + (bbn_ks_R(x, y, z)*Derivative(bbn_ks_R(x, y, z), y, z) + 
Derivative(bbn_ks_R(x, y, z), y)*Derivative(bbn_ks_R(x, y, z), z))*
pow(bbn_ks_R(x, y, z), 4))) + 8*(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4))*(2*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), y) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), y) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*
Derivative(bbn_ks_R(x, y, z), z) + (pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z))*
Derivative(bbn_ks_R(x, y, z), y))*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*
Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*
Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), z)) - 
48*(pow(a_BH, 2)*bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), y) - 
pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), y) + 
pow(bbn_ks_R(x, y, z), 4)*Derivative(bbn_ks_R(x, y, z), y))*pow(pow(a_BH, 2)*
bbn_ks_R(x, y, z)*bbn_ks_Z(x, y, z)*Derivative(bbn_ks_Z(x, y, z), z) - pow(a_BH, 2)*
pow(bbn_ks_Z(x, y, z), 2)*Derivative(bbn_ks_R(x, y, z), z) + pow(bbn_ks_R(x, y, z), 4)*
Derivative(bbn_ks_R(x, y, z), z), 2))/pow(pow(a_BH, 2)*pow(bbn_ks_Z(x, y, z), 2) + 
pow(bbn_ks_R(x, y, z), 4), 4);
}

