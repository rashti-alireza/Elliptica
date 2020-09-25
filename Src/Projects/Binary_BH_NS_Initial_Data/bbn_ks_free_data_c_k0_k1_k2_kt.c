#include "bbn_ks_free_data_analytic.h"
KS_func_def_macro(k0) KS_func_args_macro;
KS_func_def_macro(dk0_D2) KS_func_args_macro;
KS_func_def_macro(dk0_D0) KS_func_args_macro;
KS_func_def_macro(dk0_D1) KS_func_args_macro;
KS_func_def_macro(ddk0_D0D0) KS_func_args_macro;
KS_func_def_macro(ddk0_D1D1) KS_func_args_macro;
KS_func_def_macro(ddk0_D0D1) KS_func_args_macro;
KS_func_def_macro(ddk0_D2D2) KS_func_args_macro;
KS_func_def_macro(ddk0_D1D2) KS_func_args_macro;
KS_func_def_macro(ddk0_D0D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddk0_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddk0_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddk0_D0D0D2) KS_func_args_macro;
KS_func_def_macro(k1) KS_func_args_macro;
KS_func_def_macro(dk1_D0) KS_func_args_macro;
KS_func_def_macro(dk1_D2) KS_func_args_macro;
KS_func_def_macro(dk1_D1) KS_func_args_macro;
KS_func_def_macro(ddk1_D1D1) KS_func_args_macro;
KS_func_def_macro(ddk1_D2D2) KS_func_args_macro;
KS_func_def_macro(ddk1_D0D0) KS_func_args_macro;
KS_func_def_macro(ddk1_D0D1) KS_func_args_macro;
KS_func_def_macro(ddk1_D0D2) KS_func_args_macro;
KS_func_def_macro(ddk1_D1D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddk1_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddk1_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddk1_D1D2D0) KS_func_args_macro;
KS_func_def_macro(k2) KS_func_args_macro;
KS_func_def_macro(dk2_D1) KS_func_args_macro;
KS_func_def_macro(dk2_D0) KS_func_args_macro;
KS_func_def_macro(dk2_D2) KS_func_args_macro;
KS_func_def_macro(ddk2_D0D0) KS_func_args_macro;
KS_func_def_macro(ddk2_D1D2) KS_func_args_macro;
KS_func_def_macro(ddk2_D0D1) KS_func_args_macro;
KS_func_def_macro(ddk2_D1D1) KS_func_args_macro;
KS_func_def_macro(ddk2_D0D2) KS_func_args_macro;
KS_func_def_macro(ddk2_D2D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddk2_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddk2_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddk2_D1D2D2) KS_func_args_macro;
KS_func_def_macro(kt) KS_func_args_macro;
KS_func_def_macro(dkt_D0) KS_func_args_macro;
KS_func_def_macro(dkt_D2) KS_func_args_macro;
KS_func_def_macro(dkt_D1) KS_func_args_macro;
KS_func_def_macro(ddkt_D0D1) KS_func_args_macro;
KS_func_def_macro(ddkt_D2D2) KS_func_args_macro;
KS_func_def_macro(ddkt_D0D2) KS_func_args_macro;
KS_func_def_macro(ddkt_D1D2) KS_func_args_macro;
KS_func_def_macro(ddkt_D1D1) KS_func_args_macro;
KS_func_def_macro(ddkt_D0D0) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddkt_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddkt_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddkt_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddkt_D0D1D0) KS_func_args_macro;
KS_func_def_macro(c) KS_func_args_macro;
KS_func_def_macro(dc_D0) KS_func_args_macro;
KS_func_def_macro(dc_D1) KS_func_args_macro;
KS_func_def_macro(dc_D2) KS_func_args_macro;
KS_func_def_macro(ddc_D1D2) KS_func_args_macro;
KS_func_def_macro(ddc_D0D0) KS_func_args_macro;
KS_func_def_macro(ddc_D0D1) KS_func_args_macro;
KS_func_def_macro(ddc_D0D2) KS_func_args_macro;
KS_func_def_macro(ddc_D1D1) KS_func_args_macro;
KS_func_def_macro(ddc_D2D2) KS_func_args_macro;
KS_func_def_macro(dddc_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddc_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddc_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddc_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddc_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddc_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddc_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddc_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddc_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddc_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddc_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddc_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddc_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddc_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddc_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddc_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddc_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddc_D0D1D0) KS_func_args_macro;
KS_func_def_macro(k0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// bbn_ks_K0
// bbn_ks_K1
// bbn_ks_K2
pow(1 - B2, -0.5)*(-B2*Bx - (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*
sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*
bbn_ks_K1(x, y, z) + (Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*
(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - 
sqrt(1 - B2)))*cos(phiz))*cos(phiy))*bbn_ks_K0(x, y, z) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
bbn_ks_K2(x, y, z))/B2;
}
KS_func_def_macro(dk0_D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), z) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), z) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), z))/B2;
}
KS_func_def_macro(dk0_D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), x) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x))/B2;
}
KS_func_def_macro(dk0_D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), y) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), y) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), y))/B2;
}
KS_func_def_macro(ddk0_D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 2)) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (x, 2)))/B2;
}
KS_func_def_macro(ddk0_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 2)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (y, 2)) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (y, 2)))/B2;
}
KS_func_def_macro(ddk0_D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, y) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x, y))/B2;
}
KS_func_def_macro(ddk0_D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (z, 2)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (z, 2)) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (z, 2)))/B2;
}
KS_func_def_macro(ddk0_D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), y, z) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), y, z) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), y, z))/B2;
}
KS_func_def_macro(ddk0_D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), x, z) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, z) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x, z))/B2;
}
KS_func_def_macro(dddk0_D2D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (z, 2)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, (z, 2)) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (z, 2)))/
B2;
}
KS_func_def_macro(dddk0_D0D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (z, 2)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, (z, 2)) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (z, 2)))/
B2;
}
KS_func_def_macro(dddk0_D1D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 3)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (y, 3)) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (y, 3)))/B2;
}
KS_func_def_macro(dddk0_D0D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), z) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 2), z) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), z))/
B2;
}
KS_func_def_macro(dddk0_D0D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (y, 2)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, (y, 2)) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (y, 2)))/
B2;
}
KS_func_def_macro(dddk0_D0D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y, z) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, y, z) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x, y, z))/B2;
}
KS_func_def_macro(dddk0_D1D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 2), z) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (y, 2), z) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 2), z))/
B2;
}
KS_func_def_macro(dddk0_D0D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 3)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 3)) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (x, 3)))/B2;
}
KS_func_def_macro(dddk0_D1D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), y, (z, 2)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), y, (z, 2)) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), y, (z, 2)))/
B2;
}
KS_func_def_macro(dddk0_D1D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 2), z) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (y, 2), z) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 2), z))/
B2;
}
KS_func_def_macro(dddk0_D2D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), y, (z, 2)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), y, (z, 2)) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), y, (z, 2)))/
B2;
}
KS_func_def_macro(dddk0_D0D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), y) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 2), y) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), y))/
B2;
}
KS_func_def_macro(dddk0_D0D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y, z) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, y, z) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x, y, z))/B2;
}
KS_func_def_macro(dddk0_D1D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y, z) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, y, z) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x, y, z))/B2;
}
KS_func_def_macro(dddk0_D0D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), y) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 2), y) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), y))/
B2;
}
KS_func_def_macro(dddk0_D2D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (z, 3)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (z, 3)) - (Bx*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (-
B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (z, 3)))/B2;
}
KS_func_def_macro(dddk0_D1D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (y, 2)) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, (y, 2)) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (y, 2)))/
B2;
}
KS_func_def_macro(dddk0_D0D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(-(Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - 
B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*sin(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), z) + 
(Bx*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
sin(phiz) + (B2*sqrt(1 - B2) + pow(Bx, 2)*(1 - sqrt(1 - B2)))*
cos(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 2), z) - (Bx*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + 
(-B2*sqrt(1 - B2) + pow(Bx, 2)*(sqrt(1 - B2) - 1))*cos(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), z))/
B2;
}
KS_func_def_macro(k1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// bbn_ks_K0
// bbn_ks_K1
// bbn_ks_K2
pow(1 - B2, -0.5)*(-B2*By + (Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*
sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*
bbn_ks_K1(x, y, z) + (By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*
(sqrt(1 - B2) - 1)*cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - 
sqrt(1 - B2)))*sin(phiz))*cos(phiy))*bbn_ks_K0(x, y, z) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
bbn_ks_K2(x, y, z))/B2;
}
KS_func_def_macro(dk1_D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x))/B2;
}
KS_func_def_macro(dk1_D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), z) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), z) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), z))/B2;
}
KS_func_def_macro(dk1_D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), y) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), y) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), y))/B2;
}
KS_func_def_macro(ddk1_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 2)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (y, 2)) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (y, 2)))/B2;
}
KS_func_def_macro(ddk1_D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (z, 2)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (z, 2)) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (z, 2)))/B2;
}
KS_func_def_macro(ddk1_D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 2)) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (x, 2)))/B2;
}
KS_func_def_macro(ddk1_D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, y) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x, y))/B2;
}
KS_func_def_macro(ddk1_D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, z) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, z) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x, z))/B2;
}
KS_func_def_macro(ddk1_D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), y, z) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), y, z) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), y, z))/B2;
}
KS_func_def_macro(dddk1_D1D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (y, 2)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, (y, 2)) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (y, 2)))/
B2;
}
KS_func_def_macro(dddk1_D0D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (z, 2)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, (z, 2)) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (z, 2)))/
B2;
}
KS_func_def_macro(dddk1_D0D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), y) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 2), y) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), y))/
B2;
}
KS_func_def_macro(dddk1_D2D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), y, (z, 2)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), y, (z, 2)) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), y, (z, 2)))/
B2;
}
KS_func_def_macro(dddk1_D0D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y, z) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, y, z) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x, y, z))/B2;
}
KS_func_def_macro(dddk1_D0D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y, z) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, y, z) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x, y, z))/B2;
}
KS_func_def_macro(dddk1_D0D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (y, 2)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, (y, 2)) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (y, 2)))/
B2;
}
KS_func_def_macro(dddk1_D1D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 2), z) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (y, 2), z) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 2), z))/
B2;
}
KS_func_def_macro(dddk1_D1D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 3)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (y, 3)) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (y, 3)))/B2;
}
KS_func_def_macro(dddk1_D0D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), z) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 2), z) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), z))/
B2;
}
KS_func_def_macro(dddk1_D1D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), y, (z, 2)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), y, (z, 2)) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), y, (z, 2)))/
B2;
}
KS_func_def_macro(dddk1_D2D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (z, 2)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, (z, 2)) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (z, 2)))/
B2;
}
KS_func_def_macro(dddk1_D2D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (z, 3)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (z, 3)) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (z, 3)))/B2;
}
KS_func_def_macro(dddk1_D0D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 3)) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 3)) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), (x, 3)))/B2;
}
KS_func_def_macro(dddk1_D0D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), y) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 2), y) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), y))/
B2;
}
KS_func_def_macro(dddk1_D0D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), z) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (x, 2), z) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), z))/
B2;
}
KS_func_def_macro(dddk1_D1D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 2), z) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), (y, 2), z) - (By*
Bz*(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + 
(-B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*
sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 2), z))/
B2;
}
KS_func_def_macro(dddk1_D1D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*By*(sqrt(1 - B2) - 1)*sin(phiz) + (B2*sqrt(1 - 
B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y, z) + 
(By*Bz*(sqrt(1 - B2) - 1)*sin(phiy) + (-Bx*By*(sqrt(1 - B2) - 1)*
cos(phiz) + (B2*sqrt(1 - B2) + pow(By, 2)*(1 - sqrt(1 - B2)))*
sin(phiz))*cos(phiy))*Derivative(bbn_ks_K0(x, y, z), x, y, z) - (By*Bz*
(sqrt(1 - B2) - 1)*cos(phiy) + (Bx*By*(sqrt(1 - B2) - 1)*cos(phiz) + (-
B2*sqrt(1 - B2) + pow(By, 2)*(sqrt(1 - B2) - 1))*sin(phiz))*sin(phiy))*
Derivative(bbn_ks_K2(x, y, z), x, y, z))/B2;
}
KS_func_def_macro(k2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// bbn_ks_K0
// bbn_ks_K1
// bbn_ks_K2
pow(1 - B2, -0.5)*(-B2*Bz + Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - 
B2) - 1)*bbn_ks_K1(x, y, z) - (Bz*(Bx*cos(phiz) + By*sin(phiz))*
(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + pow(Bz, 2)*(sqrt(1 - 
B2) - 1))*cos(phiy))*bbn_ks_K2(x, y, z) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*sqrt(1 - B2) + 
pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*bbn_ks_K0(x, y, z))/
B2;
}
KS_func_def_macro(dk2_D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), y) - (Bz*(Bx*cos(phiz) + By*sin(phiz))*
(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + pow(Bz, 2)*(sqrt(1 - 
B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), y) - (Bz*(Bx*
cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*sqrt(1 - 
B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*Derivative(bbn_ks_K0(x, y, z), y))/
B2;
}
KS_func_def_macro(dk2_D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), x) - (Bz*(Bx*cos(phiz) + By*sin(phiz))*
(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + pow(Bz, 2)*(sqrt(1 - 
B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), x) - (Bz*(Bx*
cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*sqrt(1 - 
B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*Derivative(bbn_ks_K0(x, y, z), x))/
B2;
}
KS_func_def_macro(dk2_D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), z) - (Bz*(Bx*cos(phiz) + By*sin(phiz))*
(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + pow(Bz, 2)*(sqrt(1 - 
B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), z) - (Bz*(Bx*
cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*sqrt(1 - 
B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*Derivative(bbn_ks_K0(x, y, z), z))/
B2;
}
KS_func_def_macro(ddk2_D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (x, 2)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 2)))/B2;
}
KS_func_def_macro(ddk2_D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), y, z) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), y, z) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), y, z))/B2;
}
KS_func_def_macro(ddk2_D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), x, y) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), x, y) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, y))/B2;
}
KS_func_def_macro(ddk2_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (y, 2)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 2)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (y, 2)))/B2;
}
KS_func_def_macro(ddk2_D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), x, z) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), x, z) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, z))/B2;
}
KS_func_def_macro(ddk2_D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (z, 2)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (z, 2)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (z, 2)))/B2;
}
KS_func_def_macro(dddk2_D0D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (x, 2), y) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), y) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 2), y))/B2;
}
KS_func_def_macro(dddk2_D1D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (y, 2), z) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 2), z) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (y, 2), z))/B2;
}
KS_func_def_macro(dddk2_D0D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), x, (z, 2)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (z, 2)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, (z, 2)))/B2;
}
KS_func_def_macro(dddk2_D0D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), x, y, z) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), x, y, z) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, y, z))/B2;
}
KS_func_def_macro(dddk2_D0D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (x, 2), z) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), z) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 2), z))/B2;
}
KS_func_def_macro(dddk2_D2D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), y, (z, 2)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), y, (z, 2)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), y, (z, 2)))/B2;
}
KS_func_def_macro(dddk2_D1D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (y, 3)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 3)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (y, 3)))/B2;
}
KS_func_def_macro(dddk2_D1D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), x, y, z) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), x, y, z) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, y, z))/B2;
}
KS_func_def_macro(dddk2_D0D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (x, 3)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 3)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 3)))/B2;
}
KS_func_def_macro(dddk2_D1D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), x, (y, 2)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (y, 2)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, (y, 2)))/B2;
}
KS_func_def_macro(dddk2_D0D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), x, (y, 2)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (y, 2)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, (y, 2)))/B2;
}
KS_func_def_macro(dddk2_D2D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (z, 3)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (z, 3)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (z, 3)))/B2;
}
KS_func_def_macro(dddk2_D1D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (y, 2), z) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 2), z) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (y, 2), z))/B2;
}
KS_func_def_macro(dddk2_D0D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (x, 2), z) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), z) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 2), z))/B2;
}
KS_func_def_macro(dddk2_D2D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), x, (z, 2)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (z, 2)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, (z, 2)))/B2;
}
KS_func_def_macro(dddk2_D0D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), (x, 2), y) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), y) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 2), y))/B2;
}
KS_func_def_macro(dddk2_D0D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), x, y, z) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), x, y, z) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, y, z))/B2;
}
KS_func_def_macro(dddk2_D1D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*(Bz*(Bx*sin(phiz) - By*cos(phiz))*(sqrt(1 - B2) - 1)*
Derivative(bbn_ks_K1(x, y, z), y, (z, 2)) - (Bz*(Bx*cos(phiz) + By*
sin(phiz))*(sqrt(1 - B2) - 1)*sin(phiy) + (-B2*sqrt(1 - B2) + 
pow(Bz, 2)*(sqrt(1 - B2) - 1))*cos(phiy))*Derivative(bbn_ks_K2(x, y, z), y, (z, 2)) - 
(Bz*(Bx*cos(phiz) + By*sin(phiz))*(sqrt(1 - B2) - 1)*cos(phiy) + (B2*
sqrt(1 - B2) + pow(Bz, 2)*(1 - sqrt(1 - B2)))*sin(phiy))*
Derivative(bbn_ks_K0(x, y, z), y, (z, 2)))/B2;
}
KS_func_def_macro(kt) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// bbn_ks_K0
// bbn_ks_K1
// bbn_ks_K2
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*bbn_ks_K1(x, y, z) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
bbn_ks_K0(x, y, z) - (Bz*cos(phiy) + (Bx*cos(phiz) + By*sin(phiz))*
sin(phiy))*bbn_ks_K2(x, y, z) + 1);
}
KS_func_def_macro(dkt_D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), x) - (Bz*cos(phiy) + (Bx*cos(phiz) + By*
sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x));
}
KS_func_def_macro(dkt_D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), z) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), z) - (Bz*cos(phiy) + (Bx*cos(phiz) + By*
sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), z));
}
KS_func_def_macro(dkt_D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), y) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), y) - (Bz*cos(phiy) + (Bx*cos(phiz) + By*
sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), y));
}
KS_func_def_macro(ddkt_D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, y) - (Bz*cos(phiy) + (Bx*cos(phiz) + 
By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, y));
}
KS_func_def_macro(ddkt_D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (z, 2)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (z, 2)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (z, 2)));
}
KS_func_def_macro(ddkt_D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, z) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, z) - (Bz*cos(phiy) + (Bx*cos(phiz) + 
By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, z));
}
KS_func_def_macro(ddkt_D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), y, z) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), y, z) - (Bz*cos(phiy) + (Bx*cos(phiz) + 
By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), y, z));
}
KS_func_def_macro(ddkt_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 2)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (y, 2)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 2)));
}
KS_func_def_macro(ddkt_D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 2)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2)));
}
KS_func_def_macro(dddkt_D1D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 3)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (y, 3)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 3)));
}
KS_func_def_macro(dddkt_D1D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), y, (z, 2)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), y, (z, 2)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), y, (z, 2)));
}
KS_func_def_macro(dddkt_D0D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), y) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 2), y) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), y));
}
KS_func_def_macro(dddkt_D0D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), z) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 2), z) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), z));
}
KS_func_def_macro(dddkt_D2D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), y, (z, 2)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), y, (z, 2)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), y, (z, 2)));
}
KS_func_def_macro(dddkt_D0D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (z, 2)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, (z, 2)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (z, 2)));
}
KS_func_def_macro(dddkt_D1D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 2), z) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (y, 2), z) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 2), z));
}
KS_func_def_macro(dddkt_D1D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (y, 2), z) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (y, 2), z) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (y, 2), z));
}
KS_func_def_macro(dddkt_D0D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (y, 2)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, (y, 2)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (y, 2)));
}
KS_func_def_macro(dddkt_D1D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y, z) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, y, z) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, y, z));
}
KS_func_def_macro(dddkt_D1D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (y, 2)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, (y, 2)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (y, 2)));
}
KS_func_def_macro(dddkt_D2D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, (z, 2)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, (z, 2)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, (z, 2)));
}
KS_func_def_macro(dddkt_D2D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (z, 3)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (z, 3)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (z, 3)));
}
KS_func_def_macro(dddkt_D0D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y, z) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, y, z) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, y, z));
}
KS_func_def_macro(dddkt_D0D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), z) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 2), z) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), z));
}
KS_func_def_macro(dddkt_D0D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 3)) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 3)) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 3)));
}
KS_func_def_macro(dddkt_D0D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), x, y, z) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), x, y, z) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), x, y, z));
}
KS_func_def_macro(dddkt_D0D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
pow(1 - B2, -0.5)*((Bx*sin(phiz) - By*cos(phiz))*Derivative(bbn_ks_K1(x, y, z), (x, 2), y) + 
(Bz*sin(phiy) - (Bx*cos(phiz) + By*sin(phiz))*cos(phiy))*
Derivative(bbn_ks_K0(x, y, z), (x, 2), y) - (Bz*cos(phiy) + (Bx*
cos(phiz) + By*sin(phiz))*sin(phiy))*Derivative(bbn_ks_K2(x, y, z), (x, 2), y));
}
KS_func_def_macro(c) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*bbn_ks_rolloff(x, y, z);
}
KS_func_def_macro(dc_D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), x) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), x);
}
KS_func_def_macro(dc_D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), y) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), y);
}
KS_func_def_macro(dc_D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), z) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), z);
}
KS_func_def_macro(ddc_D1D2) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), y, z) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), y, z) + 2*
Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), z) + 
2*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), y);
}
KS_func_def_macro(ddc_D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (x, 2)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (x, 2)) + 4*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), x);
}
KS_func_def_macro(ddc_D0D1) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), x, y) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), x, y) + 2*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), y) + 
2*Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), x);
}
KS_func_def_macro(ddc_D0D2) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), x, z) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), x, z) + 2*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), z) + 
2*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), x);
}
KS_func_def_macro(ddc_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (y, 2)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (y, 2)) + 4*
Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), y);
}
KS_func_def_macro(ddc_D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (z, 2)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (z, 2)) + 4*
Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), z);
}
KS_func_def_macro(dddc_D0D2D2) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), x, (z, 2)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), x, (z, 2)) + 2*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), (z, 2)) + 
4*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), x, z) + 
2*Derivative(bbn_ks_H(x, y, z), (z, 2))*Derivative(bbn_ks_rolloff(x, y, z), x) + 
4*Derivative(bbn_ks_rolloff(x, y, z), z)*Derivative(bbn_ks_H(x, y, z), x, z);
}
KS_func_def_macro(dddc_D2D2D0) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), x, (z, 2)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), x, (z, 2)) + 2*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), (z, 2)) + 
4*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), x, z) + 
2*Derivative(bbn_ks_H(x, y, z), (z, 2))*Derivative(bbn_ks_rolloff(x, y, z), x) + 
4*Derivative(bbn_ks_rolloff(x, y, z), z)*Derivative(bbn_ks_H(x, y, z), x, z);
}
KS_func_def_macro(dddc_D1D2D2) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), y, (z, 2)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), y, (z, 2)) + 2*
Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), (z, 2)) + 
4*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), y, z) + 
2*Derivative(bbn_ks_H(x, y, z), (z, 2))*Derivative(bbn_ks_rolloff(x, y, z), y) + 
4*Derivative(bbn_ks_rolloff(x, y, z), z)*Derivative(bbn_ks_H(x, y, z), y, z);
}
KS_func_def_macro(dddc_D1D1D1) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (y, 3)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (y, 3)) + 6*
Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), (y, 2)) + 
6*Derivative(bbn_ks_H(x, y, z), (y, 2))*Derivative(bbn_ks_rolloff(x, y, z), y);
}
KS_func_def_macro(dddc_D0D0D1) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (x, 2), y) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (x, 2), y) + 4*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), x, y) + 
2*Derivative(bbn_ks_H(x, y, z), (x, 2))*Derivative(bbn_ks_rolloff(x, y, z), y) + 
2*Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), (x, 2)) + 
4*Derivative(bbn_ks_rolloff(x, y, z), x)*Derivative(bbn_ks_H(x, y, z), x, y);
}
KS_func_def_macro(dddc_D1D2D0) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), x, y, z) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), x, y, z) + 2*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), y, z) + 
2*Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), x, z) + 
2*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), x, y) + 
2*Derivative(bbn_ks_rolloff(x, y, z), x)*Derivative(bbn_ks_H(x, y, z), y, z) + 
2*Derivative(bbn_ks_rolloff(x, y, z), y)*Derivative(bbn_ks_H(x, y, z), x, z) + 
2*Derivative(bbn_ks_rolloff(x, y, z), z)*Derivative(bbn_ks_H(x, y, z), x, y);
}
KS_func_def_macro(dddc_D0D1D1) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), x, (y, 2)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), x, (y, 2)) + 2*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), (y, 2)) + 
4*Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), x, y) + 
2*Derivative(bbn_ks_H(x, y, z), (y, 2))*Derivative(bbn_ks_rolloff(x, y, z), x) + 
4*Derivative(bbn_ks_rolloff(x, y, z), y)*Derivative(bbn_ks_H(x, y, z), x, y);
}
KS_func_def_macro(dddc_D0D0D0) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (x, 3)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (x, 3)) + 6*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), (x, 2)) + 
6*Derivative(bbn_ks_H(x, y, z), (x, 2))*Derivative(bbn_ks_rolloff(x, y, z), x);
}
KS_func_def_macro(dddc_D0D2D1) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), x, y, z) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), x, y, z) + 2*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), y, z) + 
2*Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), x, z) + 
2*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), x, y) + 
2*Derivative(bbn_ks_rolloff(x, y, z), x)*Derivative(bbn_ks_H(x, y, z), y, z) + 
2*Derivative(bbn_ks_rolloff(x, y, z), y)*Derivative(bbn_ks_H(x, y, z), x, z) + 
2*Derivative(bbn_ks_rolloff(x, y, z), z)*Derivative(bbn_ks_H(x, y, z), x, y);
}
KS_func_def_macro(dddc_D1D1D0) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), x, (y, 2)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), x, (y, 2)) + 2*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), (y, 2)) + 
4*Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), x, y) + 
2*Derivative(bbn_ks_H(x, y, z), (y, 2))*Derivative(bbn_ks_rolloff(x, y, z), x) + 
4*Derivative(bbn_ks_rolloff(x, y, z), y)*Derivative(bbn_ks_H(x, y, z), x, y);
}
KS_func_def_macro(dddc_D2D2D2) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (z, 3)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (z, 3)) + 6*
Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), (z, 2)) + 
6*Derivative(bbn_ks_H(x, y, z), (z, 2))*Derivative(bbn_ks_rolloff(x, y, z), z);
}
KS_func_def_macro(dddc_D2D2D1) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), y, (z, 2)) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), y, (z, 2)) + 2*
Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), (z, 2)) + 
4*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), y, z) + 
2*Derivative(bbn_ks_H(x, y, z), (z, 2))*Derivative(bbn_ks_rolloff(x, y, z), y) + 
4*Derivative(bbn_ks_rolloff(x, y, z), z)*Derivative(bbn_ks_H(x, y, z), y, z);
}
KS_func_def_macro(dddc_D0D0D2) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (x, 2), z) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (x, 2), z) + 4*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), x, z) + 
2*Derivative(bbn_ks_H(x, y, z), (x, 2))*Derivative(bbn_ks_rolloff(x, y, z), z) + 
2*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), (x, 2)) + 
4*Derivative(bbn_ks_rolloff(x, y, z), x)*Derivative(bbn_ks_H(x, y, z), x, z);
}
KS_func_def_macro(dddc_D1D1D2) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (y, 2), z) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (y, 2), z) + 4*
Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), y, z) + 
2*Derivative(bbn_ks_H(x, y, z), (y, 2))*Derivative(bbn_ks_rolloff(x, y, z), z) + 
2*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), (y, 2)) + 
4*Derivative(bbn_ks_rolloff(x, y, z), y)*Derivative(bbn_ks_H(x, y, z), y, z);
}
KS_func_def_macro(dddc_D0D2D0) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (x, 2), z) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (x, 2), z) + 4*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), x, z) + 
2*Derivative(bbn_ks_H(x, y, z), (x, 2))*Derivative(bbn_ks_rolloff(x, y, z), z) + 
2*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), (x, 2)) + 
4*Derivative(bbn_ks_rolloff(x, y, z), x)*Derivative(bbn_ks_H(x, y, z), x, z);
}
KS_func_def_macro(dddc_D0D1D2) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), x, y, z) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), x, y, z) + 2*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), y, z) + 
2*Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), x, z) + 
2*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), x, y) + 
2*Derivative(bbn_ks_rolloff(x, y, z), x)*Derivative(bbn_ks_H(x, y, z), y, z) + 
2*Derivative(bbn_ks_rolloff(x, y, z), y)*Derivative(bbn_ks_H(x, y, z), x, z) + 
2*Derivative(bbn_ks_rolloff(x, y, z), z)*Derivative(bbn_ks_H(x, y, z), x, y);
}
KS_func_def_macro(dddc_D1D2D1) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (y, 2), z) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (y, 2), z) + 4*
Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), y, z) + 
2*Derivative(bbn_ks_H(x, y, z), (y, 2))*Derivative(bbn_ks_rolloff(x, y, z), z) + 
2*Derivative(bbn_ks_H(x, y, z), z)*Derivative(bbn_ks_rolloff(x, y, z), (y, 2)) + 
4*Derivative(bbn_ks_rolloff(x, y, z), y)*Derivative(bbn_ks_H(x, y, z), y, z);
}
KS_func_def_macro(dddc_D0D1D0) KS_func_args_macro
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
// bbn_ks_H
// bbn_ks_rolloff
2*bbn_ks_H(x, y, z)*Derivative(bbn_ks_rolloff(x, y, z), (x, 2), y) + 2*
bbn_ks_rolloff(x, y, z)*Derivative(bbn_ks_H(x, y, z), (x, 2), y) + 4*
Derivative(bbn_ks_H(x, y, z), x)*Derivative(bbn_ks_rolloff(x, y, z), x, y) + 
2*Derivative(bbn_ks_H(x, y, z), (x, 2))*Derivative(bbn_ks_rolloff(x, y, z), y) + 
2*Derivative(bbn_ks_H(x, y, z), y)*Derivative(bbn_ks_rolloff(x, y, z), (x, 2)) + 
4*Derivative(bbn_ks_rolloff(x, y, z), x)*Derivative(bbn_ks_H(x, y, z), x, y);
}

