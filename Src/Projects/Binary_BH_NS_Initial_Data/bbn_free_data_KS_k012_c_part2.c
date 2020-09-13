/* msimplify(...) */
/* msimplify(...) */
/* msimplify(...) */
#include "bbn_free_date_analytic.h"
double bbn_k0_KS_freedata KS_func_args_macro;
double bbn_dk0_D1_KS_freedata KS_func_args_macro;
double bbn_dk0_D0_KS_freedata KS_func_args_macro;
double bbn_dk0_D2_KS_freedata KS_func_args_macro;
double bbn_ddk0_D1D2_KS_freedata KS_func_args_macro;
double bbn_ddk0_D0D1_KS_freedata KS_func_args_macro;
double bbn_ddk0_D2D2_KS_freedata KS_func_args_macro;
double bbn_ddk0_D0D0_KS_freedata KS_func_args_macro;
double bbn_ddk0_D0D2_KS_freedata KS_func_args_macro;
double bbn_ddk0_D1D1_KS_freedata KS_func_args_macro;
double bbn_dddk0_D0D2D0_KS_freedata KS_func_args_macro;
double bbn_dddk0_D0D0D1_KS_freedata KS_func_args_macro;
double bbn_dddk0_D0D1D0_KS_freedata KS_func_args_macro;
double bbn_dddk0_D2D2D2_KS_freedata KS_func_args_macro;
double bbn_dddk0_D1D2D2_KS_freedata KS_func_args_macro;
double bbn_dddk0_D1D1D2_KS_freedata KS_func_args_macro;
double bbn_dddk0_D1D1D0_KS_freedata KS_func_args_macro;
double bbn_dddk0_D1D2D0_KS_freedata KS_func_args_macro;
double bbn_dddk0_D2D2D0_KS_freedata KS_func_args_macro;
double bbn_dddk0_D0D0D0_KS_freedata KS_func_args_macro;
double bbn_dddk0_D0D1D1_KS_freedata KS_func_args_macro;
double bbn_dddk0_D0D0D2_KS_freedata KS_func_args_macro;
double bbn_dddk0_D0D2D2_KS_freedata KS_func_args_macro;
double bbn_dddk0_D0D1D2_KS_freedata KS_func_args_macro;
double bbn_dddk0_D1D1D1_KS_freedata KS_func_args_macro;
double bbn_dddk0_D1D2D1_KS_freedata KS_func_args_macro;
double bbn_dddk0_D2D2D1_KS_freedata KS_func_args_macro;
double bbn_dddk0_D0D2D1_KS_freedata KS_func_args_macro;
double bbn_k1_KS_freedata KS_func_args_macro;
double bbn_dk1_D0_KS_freedata KS_func_args_macro;
double bbn_dk1_D2_KS_freedata KS_func_args_macro;
double bbn_dk1_D1_KS_freedata KS_func_args_macro;
double bbn_ddk1_D0D1_KS_freedata KS_func_args_macro;
double bbn_ddk1_D0D0_KS_freedata KS_func_args_macro;
double bbn_ddk1_D1D1_KS_freedata KS_func_args_macro;
double bbn_ddk1_D2D2_KS_freedata KS_func_args_macro;
double bbn_ddk1_D1D2_KS_freedata KS_func_args_macro;
double bbn_ddk1_D0D2_KS_freedata KS_func_args_macro;
double bbn_dddk1_D0D0D2_KS_freedata KS_func_args_macro;
double bbn_dddk1_D1D2D2_KS_freedata KS_func_args_macro;
double bbn_dddk1_D0D1D1_KS_freedata KS_func_args_macro;
double bbn_dddk1_D1D1D0_KS_freedata KS_func_args_macro;
double bbn_dddk1_D0D1D0_KS_freedata KS_func_args_macro;
double bbn_dddk1_D2D2D2_KS_freedata KS_func_args_macro;
double bbn_dddk1_D1D1D2_KS_freedata KS_func_args_macro;
double bbn_dddk1_D2D2D0_KS_freedata KS_func_args_macro;
double bbn_dddk1_D1D2D0_KS_freedata KS_func_args_macro;
double bbn_dddk1_D0D0D1_KS_freedata KS_func_args_macro;
double bbn_dddk1_D2D2D1_KS_freedata KS_func_args_macro;
double bbn_dddk1_D0D1D2_KS_freedata KS_func_args_macro;
double bbn_dddk1_D0D2D0_KS_freedata KS_func_args_macro;
double bbn_dddk1_D0D2D1_KS_freedata KS_func_args_macro;
double bbn_dddk1_D0D0D0_KS_freedata KS_func_args_macro;
double bbn_dddk1_D1D2D1_KS_freedata KS_func_args_macro;
double bbn_dddk1_D1D1D1_KS_freedata KS_func_args_macro;
double bbn_dddk1_D0D2D2_KS_freedata KS_func_args_macro;
double bbn_k2_KS_freedata KS_func_args_macro;
double bbn_dk2_D2_KS_freedata KS_func_args_macro;
double bbn_dk2_D0_KS_freedata KS_func_args_macro;
double bbn_dk2_D1_KS_freedata KS_func_args_macro;
double bbn_ddk2_D0D2_KS_freedata KS_func_args_macro;
double bbn_ddk2_D1D2_KS_freedata KS_func_args_macro;
double bbn_ddk2_D0D0_KS_freedata KS_func_args_macro;
double bbn_ddk2_D2D2_KS_freedata KS_func_args_macro;
double bbn_ddk2_D0D1_KS_freedata KS_func_args_macro;
double bbn_ddk2_D1D1_KS_freedata KS_func_args_macro;
double bbn_dddk2_D2D2D0_KS_freedata KS_func_args_macro;
double bbn_dddk2_D0D1D0_KS_freedata KS_func_args_macro;
double bbn_dddk2_D2D2D1_KS_freedata KS_func_args_macro;
double bbn_dddk2_D2D2D2_KS_freedata KS_func_args_macro;
double bbn_dddk2_D1D1D0_KS_freedata KS_func_args_macro;
double bbn_dddk2_D0D2D0_KS_freedata KS_func_args_macro;
double bbn_dddk2_D1D2D1_KS_freedata KS_func_args_macro;
double bbn_dddk2_D1D2D0_KS_freedata KS_func_args_macro;
double bbn_dddk2_D0D0D1_KS_freedata KS_func_args_macro;
double bbn_dddk2_D1D1D1_KS_freedata KS_func_args_macro;
double bbn_dddk2_D0D0D2_KS_freedata KS_func_args_macro;
double bbn_dddk2_D0D2D2_KS_freedata KS_func_args_macro;
double bbn_dddk2_D1D2D2_KS_freedata KS_func_args_macro;
double bbn_dddk2_D0D1D1_KS_freedata KS_func_args_macro;
double bbn_dddk2_D0D0D0_KS_freedata KS_func_args_macro;
double bbn_dddk2_D1D1D2_KS_freedata KS_func_args_macro;
double bbn_dddk2_D0D2D1_KS_freedata KS_func_args_macro;
double bbn_dddk2_D0D1D2_KS_freedata KS_func_args_macro;
double bbn_k0_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))/(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2));
}
double bbn_dk0_D1_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
((pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) - 2*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2);
}
double bbn_dk0_D0_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
((pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) - 2*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2);
}
double bbn_dk0_D2_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
((pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) - 2*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2);
}
double bbn_ddk0_D1D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), y)) - 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) + 8*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_ddk0_D0D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), x)) - 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) + 8*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_ddk0_D2D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), z)) - 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) + 8*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_ddk0_D0D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), x)) - 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) + 8*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_ddk0_D0D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), x)) - 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) + 8*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_ddk0_D1D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), y)) - 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) + 8*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_dddk0_D0D2D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (x, 2), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_x(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), z) + 2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 2*
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x, z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x, z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_x(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_R(x, y, z), z) + 
3*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D0D0D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (x, 2), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_x(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_x(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), y) + 2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), y) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 2*
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x, y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x, y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_x(x, y, z), y) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_R(x, y, z), y) + 
3*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*
Derivative(bbn_ks_func_R(x, y, z), y) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), y))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D0D1D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (x, 2), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_x(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_x(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), y) + 2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), y) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 2*
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x, y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x, y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_x(x, y, z), y) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_R(x, y, z), y) + 
3*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*
Derivative(bbn_ks_func_R(x, y, z), y) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), y))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D2D2D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (z, 3)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (z, 3)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (z, 3)) + 3*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 3*Derivative(bbn_ks_func_R(x, y, z), (z, 2))*
Derivative(bbn_ks_func_x(x, y, z), z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (z, 3)) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
3*(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 3*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
3*(a_BH*Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)) + 24*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 3))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D1D2D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y, (z, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y, (z, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), (z, 2))*
Derivative(bbn_ks_func_x(x, y, z), y) + 2*Derivative(bbn_ks_func_x(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, (z, 2)) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), y)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D1D1D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (y, 2), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (y, 2), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_x(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), y, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2), z) + 2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 2*
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), y, z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y, z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2))*Derivative(bbn_ks_func_R(x, y, z), z) + 
3*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D1D1D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, (y, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, (y, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_x(x, y, z), x) + 2*Derivative(bbn_ks_func_x(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (y, 2)) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 
2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D1D2D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, y, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), x, y) + Derivative(bbn_ks_func_x(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_x(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_x(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x, y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x, y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_x(x, y, z), y) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y, z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y, z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z) - 48*
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D2D2D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, (z, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, (z, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), (z, 2))*
Derivative(bbn_ks_func_x(x, y, z), x) + 2*Derivative(bbn_ks_func_x(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (z, 2)) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D0D0D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (x, 3)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 3)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 3)) + 3*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + 3*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_x(x, y, z), x)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 3)) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
3*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 3*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) + 
3*(a_BH*Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_x(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 24*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 3))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D0D1D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, (y, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, (y, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_x(x, y, z), x) + 2*Derivative(bbn_ks_func_x(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (y, 2)) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 
2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D0D0D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (x, 2), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_x(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), z) + 2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 2*
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x, z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x, z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_x(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_R(x, y, z), z) + 
3*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D0D2D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, (z, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, (z, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), (z, 2))*
Derivative(bbn_ks_func_x(x, y, z), x) + 2*Derivative(bbn_ks_func_x(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (z, 2)) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D0D1D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, y, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), x, y) + Derivative(bbn_ks_func_x(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_x(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_x(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x, y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x, y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_x(x, y, z), y) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y, z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y, z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z) - 48*
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D1D1D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (y, 3)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (y, 3)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 3)) + 3*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 3*Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_x(x, y, z), y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 3)) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
3*(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 3*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 
3*(a_BH*Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)) + 24*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 3))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D1D2D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), (y, 2), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (y, 2), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_x(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), y, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2), z) + 2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 2*
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), y, z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y, z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(2*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2))*Derivative(bbn_ks_func_R(x, y, z), z) + 
3*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D2D2D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y, (z, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y, (z, 2)) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), (z, 2))*
Derivative(bbn_ks_func_x(x, y, z), y) + 2*Derivative(bbn_ks_func_x(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, (z, 2)) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
2*(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 2*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z) - 48*(a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), y)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk0_D0D2D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, y, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), x, y) + Derivative(bbn_ks_func_x(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_x(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_x(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) - 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + (a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
(a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x, y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x, y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_x(x, y, z), y) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_x(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_x(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_x(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y, z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y, z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_x(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_x(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + (a_BH*bbn_ks_func_y(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + 3*(a_BH*
bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*bbn_ks_func_x(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (a_BH*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_x(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z) - 48*
(a_BH*bbn_ks_func_y(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_x(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_k1_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(-a_BH*bbn_ks_func_x(x, y, z) + bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))/(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2));
}
double bbn_dk1_D0_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
((pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) + 2*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2);
}
double bbn_dk1_D2_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
((pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) + 2*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2);
}
double bbn_dk1_D1_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
((pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) + 2*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2);
}
double bbn_ddk1_D0D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), x)) + 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + (a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) - 8*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_ddk1_D0D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), x)) + 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) - 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) + 8*(-a_BH*bbn_ks_func_x(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_ddk1_D1D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), y)) + 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) - 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) + 8*(-a_BH*bbn_ks_func_x(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_ddk1_D2D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), z)) + 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) - 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) + 8*(-a_BH*bbn_ks_func_x(x, y, z) + 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_ddk1_D1D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), y)) + 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + (a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) - 8*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_ddk1_D0D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 2)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x)) + 2*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) - 8*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3);
}
double bbn_dddk1_D0D0D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (x, 2), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (x, 2), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_y(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), z) + 2*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) - 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), z) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) - 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_y(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, z) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) - 3*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D1D2D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), y, (z, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y, (z, 2)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), (z, 2))*
Derivative(bbn_ks_func_y(x, y, z), y) + 2*Derivative(bbn_ks_func_y(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, (z, 2)) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), y) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) - 
2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) - 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_y(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) - 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) - 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), y)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D0D1D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, (y, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, (y, 2)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_y(x, y, z), x) + 2*Derivative(bbn_ks_func_y(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (y, 2)) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), x) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) - 
2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_y(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x, y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) - 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D1D1D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, (y, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, (y, 2)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (y, 2)) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_y(x, y, z), x) + 2*Derivative(bbn_ks_func_y(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (y, 2)) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), x) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) - 
2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_y(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x, y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) - 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D0D1D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (x, 2), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (x, 2), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_y(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_y(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), y) + 2*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), y) - 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), y) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) - 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_y(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x, y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, y) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), y) - 3*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*
Derivative(bbn_ks_func_R(x, y, z), y) + 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), y))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D2D2D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (z, 3)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (z, 3)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (z, 3)) + 3*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 3*Derivative(bbn_ks_func_R(x, y, z), (z, 2))*
Derivative(bbn_ks_func_y(x, y, z), z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (z, 3)) + 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) - 
3*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 3*(a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) - 
3*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_y(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)) + 24*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 3))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D1D1D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (y, 2), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (y, 2), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_y(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), y, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2), z) + 2*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) - 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), z) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) - 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_y(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), y, z) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) - 3*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D2D2D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, (z, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, (z, 2)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), (z, 2))*
Derivative(bbn_ks_func_y(x, y, z), x) + 2*Derivative(bbn_ks_func_y(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (z, 2)) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), x) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) - 
2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) - 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_y(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) - 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) - 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D1D2D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, y, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + Derivative(bbn_ks_func_y(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_y(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_y(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) - (-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) - (-
a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) - (-
a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x, y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x, y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_y(x, y, z), y) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_y(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y, z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y, z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_y(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_y(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) - (a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z) + 48*
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D0D0D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (x, 2), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (x, 2), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), y) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_y(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_y(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), y) + 2*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), y) - 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), y) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) - 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_y(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x, y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, y) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), y) - 3*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*
Derivative(bbn_ks_func_R(x, y, z), y) + 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), y))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D2D2D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), y, (z, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y, (z, 2)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, (z, 2)) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), (z, 2))*
Derivative(bbn_ks_func_y(x, y, z), y) + 2*Derivative(bbn_ks_func_y(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, (z, 2)) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), y) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) - 
2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) - 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_y(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) - 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) - 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), y)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D0D1D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, y, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + Derivative(bbn_ks_func_y(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_y(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_y(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) - (-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) - (-
a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) - (-
a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x, y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x, y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_y(x, y, z), y) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_y(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y, z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y, z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_y(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_y(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) - (a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z) + 48*
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D0D2D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (x, 2), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (x, 2), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), z) + 2*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + 2*Derivative(bbn_ks_func_y(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2), z) + 2*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) - 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), z) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) - 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_y(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, z) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) - 3*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D0D2D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, y, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x, y) + Derivative(bbn_ks_func_y(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_y(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_y(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, y, z) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + (a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) - (-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) - (-
a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) - (-
a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x, y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x, y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_y(x, y, z), y) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_y(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y, z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y, z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_y(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_y(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 8*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) - (a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) - 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z) + 48*
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D0D0D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (x, 3)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (x, 3)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 3)) + 3*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + 3*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*
Derivative(bbn_ks_func_y(x, y, z), x)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 3)) + 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) - 
3*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 3*(a_BH*
Derivative(bbn_ks_func_x(x, y, z), x) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) - 
3*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), (x, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (x, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_y(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)) + 24*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), x), 2) + 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 3))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D1D2D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (y, 2), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (y, 2), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2), z) + 2*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 2*Derivative(bbn_ks_func_y(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), y, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2), z) + 2*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) - 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), z) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) - 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_y(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), y, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-2*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), y, z) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_R(x, y, z), z) - 3*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*
Derivative(bbn_ks_func_R(x, y, z), z) + 2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), y) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z) + (-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), z) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D1D1D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), (y, 3)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), (y, 3)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 3)) + 3*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 3*Derivative(bbn_ks_func_R(x, y, z), (y, 2))*
Derivative(bbn_ks_func_y(x, y, z), y)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 3)) + 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) - 
3*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 3*(a_BH*
Derivative(bbn_ks_func_x(x, y, z), y) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) - 
3*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), (y, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (y, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_y(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)) + 24*
(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) - 
(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2) + 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), y) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), y) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 3))/pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_dddk1_D0D2D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_x
// bbn_ks_func_y
(pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 3)*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), x, (z, 2)) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, (z, 2)) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (z, 2)) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 2*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), (z, 2))*
Derivative(bbn_ks_func_y(x, y, z), x) + 2*Derivative(bbn_ks_func_y(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z)) + 2*pow(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2), 2)*((a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, (z, 2)) + (a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), x) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
(a_BH*Derivative(bbn_ks_func_x(x, y, z), x) - bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) - bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) - 
2*(-a_BH*Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) - 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), z) - (-a_BH*Derivative(bbn_ks_func_x(x, y, z), (z, 2)) + 
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_y(x, y, z), (z, 2)) + 
bbn_ks_func_y(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_y(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x) - 2*(-
a_BH*Derivative(bbn_ks_func_x(x, y, z), x, z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x, z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_y(x, y, z), z) + Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_y(x, y, z), x))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)) + 8*(pow(a_BH, 2) + 
pow(bbn_ks_func_R(x, y, z), 2))*(-(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) - 
2*(a_BH*bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*
bbn_ks_func_y(x, y, z))*bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, z) - 3*(a_BH*
bbn_ks_func_x(x, y, z) - bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 
(-a_BH*Derivative(bbn_ks_func_x(x, y, z), x) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), x) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x))*bbn_ks_func_R(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2) + 2*(-a_BH*
Derivative(bbn_ks_func_x(x, y, z), z) + bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_y(x, y, z), z) + bbn_ks_func_y(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z))*bbn_ks_func_R(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z))*
bbn_ks_func_R(x, y, z) + 48*(a_BH*bbn_ks_func_x(x, y, z) - 
bbn_ks_func_R(x, y, z)*bbn_ks_func_y(x, y, z))*pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(pow(a_BH, 2) + pow(bbn_ks_func_R(x, y, z), 2), 4);
}
double bbn_k2_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// bbn_ks_func_R
// bbn_ks_func_z
bbn_ks_func_z(x, y, z)/bbn_ks_func_R(x, y, z);
}
double bbn_dk2_D2_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_z
(bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_z(x, y, z), z) - 
bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(bbn_ks_func_R(x, y, z), 2);
}
double bbn_dk2_D0_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_z
(bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_z(x, y, z), x) - 
bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x))/
pow(bbn_ks_func_R(x, y, z), 2);
}
double bbn_dk2_D1_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_z
(bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_z(x, y, z), y) - 
bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y))/
pow(bbn_ks_func_R(x, y, z), 2);
}
double bbn_ddk2_D0D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(-(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), x))*
bbn_ks_func_R(x, y, z) + pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_z(x, y, z), x, z) + 2*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(bbn_ks_func_R(x, y, z), 3);
}
double bbn_ddk2_D1D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(-(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), y))*
bbn_ks_func_R(x, y, z) + pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_z(x, y, z), y, z) + 2*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(bbn_ks_func_R(x, y, z), 3);
}
double bbn_ddk2_D0D0_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_z
(-(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), x))*
bbn_ks_func_R(x, y, z) + pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_z(x, y, z), (x, 2)) + 2*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2))/pow(bbn_ks_func_R(x, y, z), 3);
}
double bbn_ddk2_D2D2_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_z
(-(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
2*Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), z))*
bbn_ks_func_R(x, y, z) + pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_z(x, y, z), (z, 2)) + 2*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/pow(bbn_ks_func_R(x, y, z), 3);
}
double bbn_ddk2_D0D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(-(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), y) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), x))*
bbn_ks_func_R(x, y, z) + pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_z(x, y, z), x, y) + 2*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y))/
pow(bbn_ks_func_R(x, y, z), 3);
}
double bbn_ddk2_D1D1_KS_freedata KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// bbn_ks_func_R
// bbn_ks_func_z
(-(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
2*Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), y))*
bbn_ks_func_R(x, y, z) + pow(bbn_ks_func_R(x, y, z), 2)*
Derivative(bbn_ks_func_z(x, y, z), (y, 2)) + 2*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))/pow(bbn_ks_func_R(x, y, z), 3);
}
double bbn_dddk2_D2D2D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_z(x, y, z), z) + pow(Derivative(bbn_ks_func_R(x, y, z), z), 2)*
Derivative(bbn_ks_func_z(x, y, z), x))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, (z, 2)) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), (z, 2)) + 
2*Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), x, z) + 
Derivative(bbn_ks_func_R(x, y, z), (z, 2))*Derivative(bbn_ks_func_z(x, y, z), x) + 
2*Derivative(bbn_ks_func_z(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), x, z))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), x, (z, 2)) - 6*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D0D1D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(2*bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_R(x, y, z), y) + 
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_z(x, y, z), y) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_z(x, y, z), x))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2), y) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_z(x, y, z), y) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), (x, 2)) + 
2*Derivative(bbn_ks_func_z(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, y))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), (x, 2), y) - 6*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), y))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D2D2D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
2*Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_z(x, y, z), z) + pow(Derivative(bbn_ks_func_R(x, y, z), z), 2)*
Derivative(bbn_ks_func_z(x, y, z), y))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, (z, 2)) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), (z, 2)) + 
2*Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), (z, 2))*Derivative(bbn_ks_func_z(x, y, z), y) + 
2*Derivative(bbn_ks_func_z(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), y, z))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), y, (z, 2)) - 6*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D2D2D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(6*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), z))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (z, 3)) + 3*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), (z, 2)) + 
3*Derivative(bbn_ks_func_R(x, y, z), (z, 2))*Derivative(bbn_ks_func_z(x, y, z), z))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), (z, 3)) - 6*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), z), 3))/pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D1D1D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_z(x, y, z), y) + pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*
Derivative(bbn_ks_func_z(x, y, z), x))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, (y, 2)) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), (y, 2)) + 
2*Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), (y, 2))*Derivative(bbn_ks_func_z(x, y, z), x) + 
2*Derivative(bbn_ks_func_z(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, y))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), x, (y, 2)) - 6*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D0D2D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(2*bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_R(x, y, z), z) + 
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_z(x, y, z), z) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_z(x, y, z), x))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2), z) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), x, z) + 
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_z(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), (x, 2)) + 
2*Derivative(bbn_ks_func_z(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, z))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), (x, 2), z) - 6*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D1D2D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(2*bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2))*Derivative(bbn_ks_func_R(x, y, z), z) + 
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*Derivative(bbn_ks_func_z(x, y, z), z) + 
2*Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_z(x, y, z), y))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2), z) + 
2*Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), (y, 2))*Derivative(bbn_ks_func_z(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), (y, 2)) + 
2*Derivative(bbn_ks_func_z(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), y, z))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), (y, 2), z) - 6*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D1D2D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_z(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), x))*
bbn_ks_func_R(x, y, z) - (bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y, z) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), x, z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), x, y) + 
Derivative(bbn_ks_func_z(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_z(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
Derivative(bbn_ks_func_z(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), x, y))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), x, y, z) - 6*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z))/pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D0D0D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(2*bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_R(x, y, z), y) + 
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_z(x, y, z), y) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_z(x, y, z), x))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2), y) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_z(x, y, z), y) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), (x, 2)) + 
2*Derivative(bbn_ks_func_z(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, y))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), (x, 2), y) - 6*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), y))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D1D1D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(6*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), y))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 3)) + 3*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), (y, 2)) + 
3*Derivative(bbn_ks_func_R(x, y, z), (y, 2))*Derivative(bbn_ks_func_z(x, y, z), y))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), (y, 3)) - 6*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 3))/pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D0D0D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(2*bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), x, z) + bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_R(x, y, z), z) + 
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_z(x, y, z), z) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_z(x, y, z), x))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2), z) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), x, z) + 
Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_z(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), (x, 2)) + 
2*Derivative(bbn_ks_func_z(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), x, z))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), (x, 2), z) - 6*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D0D2D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_z(x, y, z), z) + pow(Derivative(bbn_ks_func_R(x, y, z), z), 2)*
Derivative(bbn_ks_func_z(x, y, z), x))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, (z, 2)) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), (z, 2)) + 
2*Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), x, z) + 
Derivative(bbn_ks_func_R(x, y, z), (z, 2))*Derivative(bbn_ks_func_z(x, y, z), x) + 
2*Derivative(bbn_ks_func_z(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), x, z))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), x, (z, 2)) - 6*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D1D2D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), (z, 2)) + 2*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
2*Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_z(x, y, z), z) + pow(Derivative(bbn_ks_func_R(x, y, z), z), 2)*
Derivative(bbn_ks_func_z(x, y, z), y))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y, (z, 2)) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), (z, 2)) + 
2*Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), (z, 2))*Derivative(bbn_ks_func_z(x, y, z), y) + 
2*Derivative(bbn_ks_func_z(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), y, z))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), y, (z, 2)) - 6*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*pow(Derivative(bbn_ks_func_R(x, y, z), z), 2))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D0D1D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2)) + 2*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, y) + 
2*Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_z(x, y, z), y) + pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*
Derivative(bbn_ks_func_z(x, y, z), x))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, (y, 2)) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), (y, 2)) + 
2*Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), x, y) + 
Derivative(bbn_ks_func_R(x, y, z), (y, 2))*Derivative(bbn_ks_func_z(x, y, z), x) + 
2*Derivative(bbn_ks_func_z(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, y))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), x, (y, 2)) - 6*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*pow(Derivative(bbn_ks_func_R(x, y, z), y), 2))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D0D0D0_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(6*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 2)) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), x))*
bbn_ks_func_R(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (x, 3)) + 3*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), (x, 2)) + 
3*Derivative(bbn_ks_func_R(x, y, z), (x, 2))*Derivative(bbn_ks_func_z(x, y, z), x))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), (x, 3)) - 6*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), x), 3))/pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D1D1D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(2*bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), (y, 2))*Derivative(bbn_ks_func_R(x, y, z), z) + 
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*Derivative(bbn_ks_func_z(x, y, z), z) + 
2*Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_z(x, y, z), y))*bbn_ks_func_R(x, y, z) - 
(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), (y, 2), z) + 
2*Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), (y, 2))*Derivative(bbn_ks_func_z(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), (y, 2)) + 
2*Derivative(bbn_ks_func_z(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), y, z))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), (y, 2), z) - 6*bbn_ks_func_z(x, y, z)*
pow(Derivative(bbn_ks_func_R(x, y, z), y), 2)*Derivative(bbn_ks_func_R(x, y, z), z))/
pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D0D2D1_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_z(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), x))*
bbn_ks_func_R(x, y, z) - (bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y, z) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), x, z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), x, y) + 
Derivative(bbn_ks_func_z(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_z(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
Derivative(bbn_ks_func_z(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), x, y))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), x, y, z) - 6*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z))/pow(bbn_ks_func_R(x, y, z), 4);
}
double bbn_dddk2_D0D1D2_KS_freedata KS_func_args_macro
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
// bbn_ks_func_R
// bbn_ks_func_z
(2*(bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y, z) + bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_R(x, y, z), x, y) + Derivative(bbn_ks_func_R(x, y, z), x)*
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), z) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), z)*
Derivative(bbn_ks_func_z(x, y, z), y) + Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), x))*
bbn_ks_func_R(x, y, z) - (bbn_ks_func_z(x, y, z)*Derivative(bbn_ks_func_R(x, y, z), x, y, z) + 
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_z(x, y, z), y, z) + 
Derivative(bbn_ks_func_R(x, y, z), y)*Derivative(bbn_ks_func_z(x, y, z), x, z) + 
Derivative(bbn_ks_func_R(x, y, z), z)*Derivative(bbn_ks_func_z(x, y, z), x, y) + 
Derivative(bbn_ks_func_z(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y, z) + 
Derivative(bbn_ks_func_z(x, y, z), y)*Derivative(bbn_ks_func_R(x, y, z), x, z) + 
Derivative(bbn_ks_func_z(x, y, z), z)*Derivative(bbn_ks_func_R(x, y, z), x, y))*
pow(bbn_ks_func_R(x, y, z), 2) + pow(bbn_ks_func_R(x, y, z), 3)*
Derivative(bbn_ks_func_z(x, y, z), x, y, z) - 6*bbn_ks_func_z(x, y, z)*
Derivative(bbn_ks_func_R(x, y, z), x)*Derivative(bbn_ks_func_R(x, y, z), y)*
Derivative(bbn_ks_func_R(x, y, z), z))/pow(bbn_ks_func_R(x, y, z), 4);
}

