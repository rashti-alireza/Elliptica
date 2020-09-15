/* msimplify(X) */
/* msimplify(Y) */
/* msimplify(Z) */
#include "bbn_ks_free_date_analytic.h"
KS_func_def_macro(X) KS_func_args_macro;
KS_func_def_macro(dX_D0) KS_func_args_macro;
KS_func_def_macro(dX_D2) KS_func_args_macro;
KS_func_def_macro(dX_D1) KS_func_args_macro;
KS_func_def_macro(ddX_D1D1) KS_func_args_macro;
KS_func_def_macro(ddX_D0D2) KS_func_args_macro;
KS_func_def_macro(ddX_D0D0) KS_func_args_macro;
KS_func_def_macro(ddX_D1D2) KS_func_args_macro;
KS_func_def_macro(ddX_D0D1) KS_func_args_macro;
KS_func_def_macro(ddX_D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddX_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddX_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddX_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddX_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddX_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddX_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddX_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddX_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddX_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddX_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddX_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddX_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddX_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddX_D1D1D0) KS_func_args_macro;
KS_func_def_macro(Y) KS_func_args_macro;
KS_func_def_macro(dY_D0) KS_func_args_macro;
KS_func_def_macro(dY_D2) KS_func_args_macro;
KS_func_def_macro(dY_D1) KS_func_args_macro;
KS_func_def_macro(ddY_D1D2) KS_func_args_macro;
KS_func_def_macro(ddY_D0D1) KS_func_args_macro;
KS_func_def_macro(ddY_D1D1) KS_func_args_macro;
KS_func_def_macro(ddY_D0D2) KS_func_args_macro;
KS_func_def_macro(ddY_D2D2) KS_func_args_macro;
KS_func_def_macro(ddY_D0D0) KS_func_args_macro;
KS_func_def_macro(dddY_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddY_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddY_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddY_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddY_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddY_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddY_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddY_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddY_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddY_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddY_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddY_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddY_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddY_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddY_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddY_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddY_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddY_D0D1D1) KS_func_args_macro;
KS_func_def_macro(Z) KS_func_args_macro;
KS_func_def_macro(dZ_D2) KS_func_args_macro;
KS_func_def_macro(dZ_D1) KS_func_args_macro;
KS_func_def_macro(dZ_D0) KS_func_args_macro;
KS_func_def_macro(ddZ_D2D2) KS_func_args_macro;
KS_func_def_macro(ddZ_D0D0) KS_func_args_macro;
KS_func_def_macro(ddZ_D0D1) KS_func_args_macro;
KS_func_def_macro(ddZ_D1D2) KS_func_args_macro;
KS_func_def_macro(ddZ_D1D1) KS_func_args_macro;
KS_func_def_macro(ddZ_D0D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D2D0) KS_func_args_macro;
KS_func_def_macro(R) KS_func_args_macro;
KS_func_def_macro(dR_D2) KS_func_args_macro;
KS_func_def_macro(dR_D1) KS_func_args_macro;
KS_func_def_macro(dR_D0) KS_func_args_macro;
KS_func_def_macro(ddR_D1D2) KS_func_args_macro;
KS_func_def_macro(ddR_D1D1) KS_func_args_macro;
KS_func_def_macro(ddR_D2D2) KS_func_args_macro;
KS_func_def_macro(ddR_D0D2) KS_func_args_macro;
KS_func_def_macro(ddR_D0D0) KS_func_args_macro;
KS_func_def_macro(ddR_D0D1) KS_func_args_macro;
KS_func_def_macro(dddR_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddR_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddR_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddR_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddR_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddR_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddR_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddR_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddR_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddR_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddR_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddR_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddR_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddR_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddR_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddR_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddR_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddR_D1D2D0) KS_func_args_macro;
KS_func_def_macro(X) KS_func_args_macro
{
return
/* mcode in progress ... */
((-((-1 + Power(1 - B2,0.5) + B2)*Bx*Bz*x) - (-1 + Power(1 - B2,0.5) + 
B2)*By*Bz*y - Power(1 - B2,0.5)*B2*z + (-1 + Power(1 - B2,0.5) + B2)*
Power(Bx,2)*z + (-1 + Power(1 - B2,0.5) + B2)*Power(By,2)*z)*
Sin(phiy) + Cos(phiy)*((-((-1 + Power(1 - B2,0.5))*(Power(By,2)*x - Bx*
By*y + Bz*(Bz*x - Bx*z))) + B2*((Power(1 - B2,0.5) - Power(By,2) - 
Power(Bz,2))*x + Bx*(By*y + Bz*z)))*Cos(phiz) - (-((-1 + Power(1 - 
B2,0.5) + B2)*Bx*By*x) + (-1 + Power(1 - B2,0.5) + B2)*Power(Bx,2)*y + 
(-1 + Power(1 - B2,0.5))*Bz*(Bz*y - By*z) - B2*(Power(1 - B2,0.5)*y - 
Power(Bz,2)*y + By*Bz*z))*Sin(phiz)))/(B2*(Power(1 - B2,0.5) - 
Power(Bx,2) - Power(By,2) - Power(Bz,2)) - (-1 + Power(1 - B2,0.5))*
(Power(Bx,2) + Power(By,2) + Power(Bz,2)))
;
}
KS_func_def_macro(dX_D0) KS_func_args_macro
{
return
/* mcode in progress ... */
(-((-1 + Power(1 - B2,0.5) + B2)*Bx*Bz*Sin(phiy)) + Cos(phiy)*((B2*
(Power(1 - B2,0.5) - Power(By,2) - Power(Bz,2)) - (-1 + Power(1 - 
B2,0.5))*(Power(By,2) + Power(Bz,2)))*Cos(phiz) + (-1 + Power(1 - 
B2,0.5) + B2)*Bx*By*Sin(phiz)))/(B2*(Power(1 - B2,0.5) - Power(Bx,2) - 
Power(By,2) - Power(Bz,2)) - (-1 + Power(1 - B2,0.5))*(Power(Bx,2) + 
Power(By,2) + Power(Bz,2)))
;
}
KS_func_def_macro(dX_D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(((-1 + Power(1 - B2,0.5))*(Power(Bx,2) + Power(By,2)) + B2*(-Power(1 - 
B2,0.5) + Power(Bx,2) + Power(By,2)))*Sin(phiy) + (-1 + Power(1 - 
B2,0.5) + B2)*Bz*Cos(phiy)*(Bx*Cos(phiz) + By*Sin(phiz)))/(B2*
(Power(1 - B2,0.5) - Power(Bx,2) - Power(By,2) - Power(Bz,2)) - (-1 + 
Power(1 - B2,0.5))*(Power(Bx,2) + Power(By,2) + Power(Bz,2)))
;
}
KS_func_def_macro(dX_D1) KS_func_args_macro
{
return
/* mcode in progress ... */
(-((-1 + Power(1 - B2,0.5) + B2)*By*Bz*Sin(phiy)) + Cos(phiy)*((-1 + 
Power(1 - B2,0.5) + B2)*Bx*By*Cos(phiz) + (B2*(Power(1 - B2,0.5) - 
Power(Bx,2) - Power(Bz,2)) - (-1 + Power(1 - B2,0.5))*(Power(Bx,2) + 
Power(Bz,2)))*Sin(phiz)))/(B2*(Power(1 - B2,0.5) - Power(Bx,2) - 
Power(By,2) - Power(Bz,2)) - (-1 + Power(1 - B2,0.5))*(Power(Bx,2) + 
Power(By,2) + Power(Bz,2)))
;
}
KS_func_def_macro(ddX_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddX_D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddX_D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddX_D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddX_D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddX_D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D2D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D0D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D0D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D0D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D1D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D1D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D1D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D2D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D0D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D0D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D0D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D1D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D2D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D1D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D0D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D0D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D0D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddX_D1D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(Y) KS_func_args_macro
{
return
/* mcode in progress ... */
(((-1 + Power(1 - B2,0.5) + B2)*Bx*By*x - (-1 + Power(1 - B2,0.5) + 
B2)*Power(Bx,2)*y - (-1 + Power(1 - B2,0.5))*Bz*(Bz*y - By*z) + B2*
(Power(1 - B2,0.5)*y - Power(Bz,2)*y + By*Bz*z))*Cos(phiz) + ((-1 + 
Power(1 - B2,0.5))*(Power(By,2)*x - Bx*By*y + Bz*(Bz*x - Bx*z)) - B2*
((Power(1 - B2,0.5) - Power(By,2) - Power(Bz,2))*x + Bx*(By*y + Bz*
z)))*Sin(phiz))/(B2*(Power(1 - B2,0.5) - Power(Bx,2) - Power(By,2) - 
Power(Bz,2)) - (-1 + Power(1 - B2,0.5))*(Power(Bx,2) + Power(By,2) + 
Power(Bz,2)))
;
}
KS_func_def_macro(dY_D0) KS_func_args_macro
{
return
/* mcode in progress ... */
((-1 + Power(1 - B2,0.5) + B2)*Bx*By*Cos(phiz) + ((-1 + Power(1 - 
B2,0.5))*(Power(By,2) + Power(Bz,2)) + B2*(-Power(1 - B2,0.5) + 
Power(By,2) + Power(Bz,2)))*Sin(phiz))/(B2*(Power(1 - B2,0.5) - 
Power(Bx,2) - Power(By,2) - Power(Bz,2)) - (-1 + Power(1 - B2,0.5))*
(Power(Bx,2) + Power(By,2) + Power(Bz,2)))
;
}
KS_func_def_macro(dY_D2) KS_func_args_macro
{
return
/* mcode in progress ... */
((-1 + Power(1 - B2,0.5) + B2)*Bz*(By*Cos(phiz) - Bx*Sin(phiz)))/(B2*
(Power(1 - B2,0.5) - Power(Bx,2) - Power(By,2) - Power(Bz,2)) - (-1 + 
Power(1 - B2,0.5))*(Power(Bx,2) + Power(By,2) + Power(Bz,2)))
;
}
KS_func_def_macro(dY_D1) KS_func_args_macro
{
return
/* mcode in progress ... */
((B2*(Power(1 - B2,0.5) - Power(Bx,2) - Power(Bz,2)) - (-1 + Power(1 - 
B2,0.5))*(Power(Bx,2) + Power(Bz,2)))*Cos(phiz) - (-1 + Power(1 - 
B2,0.5) + B2)*Bx*By*Sin(phiz))/(B2*(Power(1 - B2,0.5) - Power(Bx,2) - 
Power(By,2) - Power(Bz,2)) - (-1 + Power(1 - B2,0.5))*(Power(Bx,2) + 
Power(By,2) + Power(Bz,2)))
;
}
KS_func_def_macro(ddY_D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddY_D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddY_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddY_D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddY_D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddY_D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D1D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D1D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D2D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D0D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D1D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D0D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D0D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D0D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D1D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D0D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D0D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D0D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D2D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D2D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D0D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D1D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D1D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddY_D0D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(Z) KS_func_args_macro
{
return
/* mcode in progress ... */
(((-1 + Power(1 - B2,0.5) + B2)*Bx*Bz*x + (-1 + Power(1 - B2,0.5) + 
B2)*By*Bz*y + Power(1 - B2,0.5)*B2*z - (-1 + Power(1 - B2,0.5) + B2)*
Power(Bx,2)*z - (-1 + Power(1 - B2,0.5) + B2)*Power(By,2)*z)*
Cos(phiy) + Sin(phiy)*((-((-1 + Power(1 - B2,0.5))*(Power(By,2)*x - Bx*
By*y + Bz*(Bz*x - Bx*z))) + B2*((Power(1 - B2,0.5) - Power(By,2) - 
Power(Bz,2))*x + Bx*(By*y + Bz*z)))*Cos(phiz) - (-((-1 + Power(1 - 
B2,0.5) + B2)*Bx*By*x) + (-1 + Power(1 - B2,0.5) + B2)*Power(Bx,2)*y + 
(-1 + Power(1 - B2,0.5))*Bz*(Bz*y - By*z) - B2*(Power(1 - B2,0.5)*y - 
Power(Bz,2)*y + By*Bz*z))*Sin(phiz)))/(B2*(Power(1 - B2,0.5) - 
Power(Bx,2) - Power(By,2) - Power(Bz,2)) - (-1 + Power(1 - B2,0.5))*
(Power(Bx,2) + Power(By,2) + Power(Bz,2)))
;
}
KS_func_def_macro(dZ_D2) KS_func_args_macro
{
return
/* mcode in progress ... */
((B2*(Power(1 - B2,0.5) - Power(Bx,2) - Power(By,2)) - (-1 + Power(1 - 
B2,0.5))*(Power(Bx,2) + Power(By,2)))*Cos(phiy) + (-1 + Power(1 - 
B2,0.5) + B2)*Bz*Sin(phiy)*(Bx*Cos(phiz) + By*Sin(phiz)))/(B2*
(Power(1 - B2,0.5) - Power(Bx,2) - Power(By,2) - Power(Bz,2)) - (-1 + 
Power(1 - B2,0.5))*(Power(Bx,2) + Power(By,2) + Power(Bz,2)))
;
}
KS_func_def_macro(dZ_D1) KS_func_args_macro
{
return
/* mcode in progress ... */
((-1 + Power(1 - B2,0.5) + B2)*By*Bz*Cos(phiy) + Sin(phiy)*((-1 + 
Power(1 - B2,0.5) + B2)*Bx*By*Cos(phiz) - ((-1 + Power(1 - B2,0.5))*
(Power(Bx,2) + Power(Bz,2)) + B2*(-Power(1 - B2,0.5) + Power(Bx,2) + 
Power(Bz,2)))*Sin(phiz)))/(B2*(Power(1 - B2,0.5) - Power(Bx,2) - 
Power(By,2) - Power(Bz,2)) - (-1 + Power(1 - B2,0.5))*(Power(Bx,2) + 
Power(By,2) + Power(Bz,2)))
;
}
KS_func_def_macro(dZ_D0) KS_func_args_macro
{
return
/* mcode in progress ... */
((-1 + Power(1 - B2,0.5) + B2)*Bx*Bz*Cos(phiy) + Sin(phiy)*((B2*
(Power(1 - B2,0.5) - Power(By,2) - Power(Bz,2)) - (-1 + Power(1 - 
B2,0.5))*(Power(By,2) + Power(Bz,2)))*Cos(phiz) + (-1 + Power(1 - 
B2,0.5) + B2)*Bx*By*Sin(phiz)))/(B2*(Power(1 - B2,0.5) - Power(Bx,2) - 
Power(By,2) - Power(Bz,2)) - (-1 + Power(1 - B2,0.5))*(Power(Bx,2) + 
Power(By,2) + Power(Bz,2)))
;
}
KS_func_def_macro(ddZ_D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddZ_D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddZ_D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddZ_D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddZ_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddZ_D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D0D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D1D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D0D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D0D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D1D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D0D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D1D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D0D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D0D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D0D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D1D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D2D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D1D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D0D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D0D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D2D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D2D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(dddZ_D1D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(R) KS_func_args_macro
{
return
/* mcode in progress ... */
-0.5*Power(Pattern(a,Blank(BH)),2) + 0.5*Power(XX(x,y,z),2) + 0.5*
Power(YY(x,y,z),2) + 0.5*Power(ZZ(x,y,z),2) + 0.5*Sqrt(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dR_D2) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),z))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
1.*Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z) + (0.5*(4*Hold(D(ZZ(x,y,z),z))*
Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*(Hold(D(XX(x,y,z),z))*
XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dR_D1) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),y))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
1.*Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z) + (0.5*(4*Hold(D(ZZ(x,y,z),y))*
Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*(Hold(D(XX(x,y,z),y))*
XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dR_D0) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),x))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
1.*Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z) + (0.5*(4*Hold(D(ZZ(x,y,z),x))*
Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*(Hold(D(XX(x,y,z),x))*
XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(ddR_D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),z)) + 1.*Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),z)) + 1.*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + 
1.*Hold(D(XX(x,y,z),y,z))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),y,z))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),y,z))*ZZ(x,y,z) + (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),y,z))*XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(ddR_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Power(Hold(D(XX(x,y,z),y)),2) + 1.*Power(Hold(D(YY(x,y,z),y)),2) + 
1.*Power(Hold(D(ZZ(x,y,z),y)),2) + 1.*Hold(D(XX(x,y,z),List(y,2)))*
XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(y,2)))*YY(x,y,z) + 1.*
Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z) + (0.5*(-4*Hold(D(ZZ(x,y,z),y))*
Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*(Hold(D(XX(x,y,z),y))*
XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*Hold(D(ZZ(x,y,z),y))*
Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*(Hold(D(XX(x,y,z),y))*
XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Power(Hold(D(ZZ(x,y,z),y)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(ddR_D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Power(Hold(D(XX(x,y,z),z)),2) + 1.*Power(Hold(D(YY(x,y,z),z)),2) + 
1.*Power(Hold(D(ZZ(x,y,z),z)),2) + 1.*Hold(D(XX(x,y,z),List(z,2)))*
XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(z,2)))*YY(x,y,z) + 1.*
Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z) + (0.5*(-4*Hold(D(ZZ(x,y,z),z))*
Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*(Hold(D(XX(x,y,z),z))*
XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*Hold(D(ZZ(x,y,z),z))*
Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*(Hold(D(XX(x,y,z),z))*
XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Power(Hold(D(ZZ(x,y,z),z)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(ddR_D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),z)) + 1.*Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),z)) + 1.*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + 
1.*Hold(D(XX(x,y,z),x,z))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),x,z))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),x,z))*ZZ(x,y,z) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),x,z))*XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(ddR_D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Power(Hold(D(XX(x,y,z),x)),2) + 1.*Power(Hold(D(YY(x,y,z),x)),2) + 
1.*Power(Hold(D(ZZ(x,y,z),x)),2) + 1.*Hold(D(XX(x,y,z),List(x,2)))*
XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(x,2)))*YY(x,y,z) + 1.*
Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z) + (0.5*(-4*Hold(D(ZZ(x,y,z),x))*
Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*(Hold(D(XX(x,y,z),x))*
XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*Hold(D(ZZ(x,y,z),x))*
Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*(Hold(D(XX(x,y,z),x))*
XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Power(Hold(D(ZZ(x,y,z),x)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(ddR_D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),y)) + 1.*Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y)) + 1.*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + 
1.*Hold(D(XX(x,y,z),x,y))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),x,y))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),x,y))*ZZ(x,y,z) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + 
Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + 
Hold(D(XX(x,y,z),x,y))*XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D0D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),List(x,2))) + 1.*
Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),List(x,2))) + 1.*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(x,2))) + 2.*
Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),x,z)) + 2.*Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),x,z)) + 2.*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,z)) + 
1.*Hold(D(XX(x,y,z),List(x,2),z))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(x,2),z))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),List(x,2),z))*ZZ(x,y,z) + (0.5*Power(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)),2)*(12*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 6*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) - (0.5*(4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),x)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) - (1.*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),x,z))*XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),List(x,2),z))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),x,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,z))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),List(x,2))) + 
Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),List(x,2))) + 
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(x,2))) + 2*
Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),x,z)) + 2*Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),x,z)) + 2*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,z)) + 
Hold(D(XX(x,y,z),List(x,2),z))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2),z))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D0D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),List(z,2))) + 1.*
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),List(z,2))) + 1.*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(z,2))) + 2.*
Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),x,z)) + 2.*Hold(D(YY(x,y,z),z))*
Hold(D(YY(x,y,z),x,z)) + 2.*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,z)) + 
1.*Hold(D(XX(x,y,z),x,List(z,2)))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),x,List(z,2)))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),x,List(z,2)))*ZZ(x,y,z) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-12*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 6*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Power(Hold(D(ZZ(x,y,z),z)),2)*Power(Pattern(a,Blank(BH)),2) - 4*
Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 
4*Power(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z),2) - 2*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (1.*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),x,z))*XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,List(z,2)))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),x,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,z))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),List(z,2))) + 
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),List(z,2))) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(z,2))) + 2*
Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),x,z)) + 2*Hold(D(YY(x,y,z),z))*
Hold(D(YY(x,y,z),x,z)) + 2*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,z)) + 
Hold(D(XX(x,y,z),x,List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),x,List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x,List(z,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D2D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
3.*Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),List(z,2))) + 3.*
Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),List(z,2))) + 3.*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(z,2))) + 1.*
Hold(D(XX(x,y,z),List(z,3)))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(z,3)))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),List(z,3)))*ZZ(x,y,z) + (0.5*(4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Power(Hold(D(ZZ(x,y,z),z)),2)*Power(Pattern(a,Blank(BH)),2) - 4*
Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 
4*Power(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z),2) - 2*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (1.*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),z)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (1.*(6*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2) + 
2*Hold(D(ZZ(x,y,z),List(z,3)))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 6*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z)) + (3*
Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),List(z,2))) + 3*
Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),List(z,2))) + 3*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(z,2))) + 
Hold(D(XX(x,y,z),List(z,3)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,3)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,3)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2)) + 
(12.*Power(Hold(D(XX(x,y,z),z))*XX(x,y,z)*(-1.*Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z)*(-1.*Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)*(Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)),3))/
(Sqrt(Power(Pattern(a,Blank(BH)),4) - 2*Power(Pattern(a,Blank(BH)),2)*
(Power(XX(x,y,z),2) + Power(YY(x,y,z),2) - Power(ZZ(x,y,z),2)) + 
Power(Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))*
Power(Power(Pattern(a,Blank(BH)),4) + Power(XX(x,y,z),4) + 
Power(YY(x,y,z),4) + 2.*Power(YY(x,y,z),2)*Power(ZZ(x,y,z),2) + 
Power(ZZ(x,y,z),4) + Power(Pattern(a,Blank(BH)),2)*(-2.*
Power(XX(x,y,z),2) - 2.*Power(YY(x,y,z),2) + 2.*Power(ZZ(x,y,z),2)) + 
Power(XX(x,y,z),2)*(2.*Power(YY(x,y,z),2) + 2.*Power(ZZ(x,y,z),2)),2))
;
}
KS_func_def_macro(dddR_D2D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),List(z,2))) + 1.*
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),List(z,2))) + 1.*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(z,2))) + 2.*
Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),x,z)) + 2.*Hold(D(YY(x,y,z),z))*
Hold(D(YY(x,y,z),x,z)) + 2.*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,z)) + 
1.*Hold(D(XX(x,y,z),x,List(z,2)))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),x,List(z,2)))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),x,List(z,2)))*ZZ(x,y,z) + (0.5*(12*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 6*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*Power(4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)),2))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) - (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),z)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) - (1.*(4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),x,z))*XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,List(z,2)))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),x,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,z))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),List(z,2))) + 
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),List(z,2))) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(z,2))) + 2*
Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),x,z)) + 2*Hold(D(YY(x,y,z),z))*
Hold(D(YY(x,y,z),x,z)) + 2*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,z)) + 
Hold(D(XX(x,y,z),x,List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),x,List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x,List(z,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D1D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),List(y,2))) + 1.*
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),List(y,2))) + 1.*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(y,2))) + 2.*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),x,y)) + 2.*Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),x,y)) + 2.*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,y)) + 
1.*Hold(D(XX(x,y,z),x,List(y,2)))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),x,List(y,2)))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),x,List(y,2)))*ZZ(x,y,z) + (0.5*(12*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 6*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*Power(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)),2))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) - (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),y)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) - (1.*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + 
Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + 
Hold(D(XX(x,y,z),x,y))*XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,List(y,2)))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),y)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + Hold(D(XX(x,y,z),x,y))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,y))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),List(y,2))) + 
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),List(y,2))) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(y,2))) + 2*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),x,y)) + 2*Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),x,y)) + 2*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,y)) + 
Hold(D(XX(x,y,z),x,List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),x,List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x,List(y,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D2D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),List(z,2))) + 1.*
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),List(z,2))) + 1.*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(z,2))) + 2.*
Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),y,z)) + 2.*Hold(D(YY(x,y,z),z))*
Hold(D(YY(x,y,z),y,z)) + 2.*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),y,z)) + 
1.*Hold(D(XX(x,y,z),y,List(z,2)))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),y,List(z,2)))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),y,List(z,2)))*ZZ(x,y,z) + (0.5*(12*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 6*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*Power(4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)),2))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) - (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),z)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) - (1.*(4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),y,z))*XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),y,List(z,2)))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),y,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y,z))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),List(z,2))) + 
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),List(z,2))) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(z,2))) + 2*
Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),y,z)) + 2*Hold(D(YY(x,y,z),z))*
Hold(D(YY(x,y,z),y,z)) + 2*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),y,z)) + 
Hold(D(XX(x,y,z),y,List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),y,List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y,List(z,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D0D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),List(y,2))) + 1.*
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),List(y,2))) + 1.*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(y,2))) + 2.*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),x,y)) + 2.*Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),x,y)) + 2.*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,y)) + 
1.*Hold(D(XX(x,y,z),x,List(y,2)))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),x,List(y,2)))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),x,List(y,2)))*ZZ(x,y,z) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-12*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 6*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Power(Hold(D(ZZ(x,y,z),y)),2)*Power(Pattern(a,Blank(BH)),2) - 4*
Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 
4*Power(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z),2) - 2*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (1.*(-4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + 
Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + 
Hold(D(XX(x,y,z),x,y))*XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,List(y,2)))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),y)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + Hold(D(XX(x,y,z),x,y))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,y))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),List(y,2))) + 
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),List(y,2))) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(y,2))) + 2*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),x,y)) + 2*Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),x,y)) + 2*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,y)) + 
Hold(D(XX(x,y,z),x,List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),x,List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x,List(y,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D1D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),List(z,2))) + 1.*
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),List(z,2))) + 1.*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(z,2))) + 2.*
Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),y,z)) + 2.*Hold(D(YY(x,y,z),z))*
Hold(D(YY(x,y,z),y,z)) + 2.*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),y,z)) + 
1.*Hold(D(XX(x,y,z),y,List(z,2)))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),y,List(z,2)))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),y,List(z,2)))*ZZ(x,y,z) + (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-12*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 6*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Power(Hold(D(ZZ(x,y,z),z)),2)*Power(Pattern(a,Blank(BH)),2) - 4*
Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 
4*Power(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z),2) - 2*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (1.*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),y,z))*XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(z,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),y,List(z,2)))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),z)),2) + 
Power(Hold(D(YY(x,y,z),z)),2) + Power(Hold(D(ZZ(x,y,z),z)),2) + 
Hold(D(XX(x,y,z),List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(z,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),y,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y,z))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),List(z,2))) + 
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),List(z,2))) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(z,2))) + 2*
Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),y,z)) + 2*Hold(D(YY(x,y,z),z))*
Hold(D(YY(x,y,z),y,z)) + 2*Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),y,z)) + 
Hold(D(XX(x,y,z),y,List(z,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),y,List(z,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y,List(z,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D0D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
3.*Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),List(x,2))) + 3.*
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),List(x,2))) + 3.*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(x,2))) + 1.*
Hold(D(XX(x,y,z),List(x,3)))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(x,3)))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),List(x,3)))*ZZ(x,y,z) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Power(Hold(D(ZZ(x,y,z),x)),2)*Power(Pattern(a,Blank(BH)),2) - 4*
Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 
4*Power(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z),2) - 2*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (1.*(-4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),x)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (1.*(6*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2) + 
2*Hold(D(ZZ(x,y,z),List(x,3)))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 6*(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z)) + (3*
Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),List(x,2))) + 3*
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),List(x,2))) + 3*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),List(x,2))) + 
Hold(D(XX(x,y,z),List(x,3)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,3)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,3)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2)) + 
(12.*Power(Hold(D(XX(x,y,z),x))*XX(x,y,z)*(-1.*Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)) + 
Hold(D(YY(x,y,z),x))*YY(x,y,z)*(-1.*Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z)*(Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)),3))/
(Sqrt(Power(Pattern(a,Blank(BH)),4) - 2*Power(Pattern(a,Blank(BH)),2)*
(Power(XX(x,y,z),2) + Power(YY(x,y,z),2) - Power(ZZ(x,y,z),2)) + 
Power(Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))*
Power(Power(Pattern(a,Blank(BH)),4) + Power(XX(x,y,z),4) + 
Power(YY(x,y,z),4) + 2.*Power(YY(x,y,z),2)*Power(ZZ(x,y,z),2) + 
Power(ZZ(x,y,z),4) + Power(Pattern(a,Blank(BH)),2)*(-2.*
Power(XX(x,y,z),2) - 2.*Power(YY(x,y,z),2) + 2.*Power(ZZ(x,y,z),2)) + 
Power(XX(x,y,z),2)*(2.*Power(YY(x,y,z),2) + 2.*Power(ZZ(x,y,z),2)),2))
;
}
KS_func_def_macro(dddR_D0D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),List(x,2))) + 1.*
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),List(x,2))) + 1.*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(x,2))) + 2.*
Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),x,y)) + 2.*Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),x,y)) + 2.*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,y)) + 
1.*Hold(D(XX(x,y,z),List(x,2),y))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(x,2),y))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),List(x,2),y))*ZZ(x,y,z) + (0.5*(-12*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 6*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),x)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2) - 
4*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + 
Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z)) - 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + 
Hold(D(XX(x,y,z),x,y))*XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + 
Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + 
Hold(D(XX(x,y,z),x,y))*XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),List(x,2),y))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),y)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + Hold(D(XX(x,y,z),x,y))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,y))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),List(x,2))) + 
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),List(x,2))) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(x,2))) + 2*
Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),x,y)) + 2*Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),x,y)) + 2*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,y)) + 
Hold(D(XX(x,y,z),List(x,2),y))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2),y))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D0D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),List(x,2))) + 1.*
Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),List(x,2))) + 1.*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(x,2))) + 2.*
Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),x,z)) + 2.*Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),x,z)) + 2.*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,z)) + 
1.*Hold(D(XX(x,y,z),List(x,2),z))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(x,2),z))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),List(x,2),z))*ZZ(x,y,z) + (0.5*(-12*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 6*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),x)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) - 
4*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) - 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),x,z))*XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),x,z))*XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),List(x,2),z))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),x,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,z))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),List(x,2))) + 
Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),List(x,2))) + 
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(x,2))) + 2*
Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),x,z)) + 2*Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),x,z)) + 2*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,z)) + 
Hold(D(XX(x,y,z),List(x,2),z))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2),z))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D0D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),List(x,2))) + 1.*
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),List(x,2))) + 1.*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(x,2))) + 2.*
Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),x,y)) + 2.*Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),x,y)) + 2.*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,y)) + 
1.*Hold(D(XX(x,y,z),List(x,2),y))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(x,2),y))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),List(x,2),y))*ZZ(x,y,z) + (0.5*Power(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)),2)*(12*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 6*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) - (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),x)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) - (1.*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + 
Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + 
Hold(D(XX(x,y,z),x,y))*XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(x,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),List(x,2),y))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),x)),2) + 
Power(Hold(D(YY(x,y,z),x)),2) + Power(Hold(D(ZZ(x,y,z),x)),2) + 
Hold(D(XX(x,y,z),List(x,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),y)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + Hold(D(XX(x,y,z),x,y))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,y))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),List(x,2))) + 
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),List(x,2))) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(x,2))) + 2*
Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),x,y)) + 2*Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),x,y)) + 2*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),x,y)) + 
Hold(D(XX(x,y,z),List(x,2),y))*XX(x,y,z) + Hold(D(YY(x,y,z),List(x,2),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(x,2),y))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D0D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),x,y)) + 1.*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),x,z)) + 1.*Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y,z)) + 1.*Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),x,y)) + 
1.*Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),x,z)) + 1.*
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),y,z)) + 1.*Hold(D(ZZ(x,y,z),z))*
Hold(D(ZZ(x,y,z),x,y)) + 1.*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,z)) + 
1.*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y,z)) + 1.*
Hold(D(XX(x,y,z),x,y,z))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),x,y,z))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),x,y,z))*ZZ(x,y,z) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-12*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 6*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + 
Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + 
Hold(D(XX(x,y,z),x,y))*XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),x,z))*XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) - 
4*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 4*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) - 2*
(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),y,z))*XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),y)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + Hold(D(XX(x,y,z),x,y))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,y))*
ZZ(x,y,z)) + 4*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),x,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,z))*
ZZ(x,y,z)) + 4*(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),y,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y,z))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),x,y)) + 
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),x,z)) + Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y,z)) + Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),x,y)) + 
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),x,z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y,z)) + Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,y)) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,z)) + Hold(D(ZZ(x,y,z),x))*
Hold(D(ZZ(x,y,z),y,z)) + Hold(D(XX(x,y,z),x,y,z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),x,y,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,y,z))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D0D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),x,y)) + 1.*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),x,z)) + 1.*Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y,z)) + 1.*Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),x,y)) + 
1.*Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),x,z)) + 1.*
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),y,z)) + 1.*Hold(D(ZZ(x,y,z),z))*
Hold(D(ZZ(x,y,z),x,y)) + 1.*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,z)) + 
1.*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y,z)) + 1.*
Hold(D(XX(x,y,z),x,y,z))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),x,y,z))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),x,y,z))*ZZ(x,y,z) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-12*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 6*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + 
Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + 
Hold(D(XX(x,y,z),x,y))*XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),x,z))*XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) - 
4*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 4*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) - 2*
(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),y,z))*XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),y)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + Hold(D(XX(x,y,z),x,y))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,y))*
ZZ(x,y,z)) + 4*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),x,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,z))*
ZZ(x,y,z)) + 4*(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),y,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y,z))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),x,y)) + 
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),x,z)) + Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y,z)) + Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),x,y)) + 
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),x,z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y,z)) + Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,y)) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,z)) + Hold(D(ZZ(x,y,z),x))*
Hold(D(ZZ(x,y,z),y,z)) + Hold(D(XX(x,y,z),x,y,z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),x,y,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,y,z))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D1D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),List(y,2))) + 1.*
Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),List(y,2))) + 1.*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(y,2))) + 2.*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),y,z)) + 2.*Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),y,z)) + 2.*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),y,z)) + 
1.*Hold(D(XX(x,y,z),List(y,2),z))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(y,2),z))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),List(y,2),z))*ZZ(x,y,z) + (0.5*Power(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)),2)*(12*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 6*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) - (0.5*(4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),y)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) - (1.*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),y,z))*XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),List(y,2),z))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),y,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y,z))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),List(y,2))) + 
Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),List(y,2))) + 
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(y,2))) + 2*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),y,z)) + 2*Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),y,z)) + 2*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),y,z)) + 
Hold(D(XX(x,y,z),List(y,2),z))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2),z))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D1D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),List(y,2))) + 1.*
Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),List(y,2))) + 1.*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(y,2))) + 2.*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),y,z)) + 2.*Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),y,z)) + 2.*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),y,z)) + 
1.*Hold(D(XX(x,y,z),List(y,2),z))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(y,2),z))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),List(y,2),z))*ZZ(x,y,z) + (0.5*(-12*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 6*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),y)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) - 
4*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 4*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) - 2*
(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),y,z))*XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),y,z))*XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2) + 
8*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),List(y,2),z))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 4*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z)) + 8*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),y,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y,z))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),List(y,2))) + 
Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),List(y,2))) + 
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),List(y,2))) + 2*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),y,z)) + 2*Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),y,z)) + 2*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),y,z)) + 
Hold(D(XX(x,y,z),List(y,2),z))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2),z))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2),z))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}
KS_func_def_macro(dddR_D1D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
3.*Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),List(y,2))) + 3.*
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),List(y,2))) + 3.*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(y,2))) + 1.*
Hold(D(XX(x,y,z),List(y,3)))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),List(y,3)))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),List(y,3)))*ZZ(x,y,z) + (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Power(Hold(D(ZZ(x,y,z),y)),2)*Power(Pattern(a,Blank(BH)),2) - 4*
Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 
4*Power(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z),2) - 2*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (1.*(-4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Power(Hold(D(ZZ(x,y,z),y)),2)*Power(Pattern(a,Blank(BH)),2) + 4*
Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 
4*Power(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z),2) + 2*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Power(4*
Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (1.*(6*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(y,2)))*Power(Pattern(a,Blank(BH)),2) + 
2*Hold(D(ZZ(x,y,z),List(y,3)))*Power(Pattern(a,Blank(BH)),2)*
ZZ(x,y,z) + 6*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Power(Hold(D(XX(x,y,z),y)),2) + 
Power(Hold(D(YY(x,y,z),y)),2) + Power(Hold(D(ZZ(x,y,z),y)),2) + 
Hold(D(XX(x,y,z),List(y,2)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,2)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,2)))*ZZ(x,y,z)) + (3*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),List(y,2))) + 3*
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),List(y,2))) + 3*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),List(y,2))) + 
Hold(D(XX(x,y,z),List(y,3)))*XX(x,y,z) + Hold(D(YY(x,y,z),List(y,3)))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),List(y,3)))*ZZ(x,y,z))*(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2)) + 
(12.*Power(Hold(D(XX(x,y,z),y))*XX(x,y,z)*(-1.*Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)) + 
Hold(D(YY(x,y,z),y))*YY(x,y,z)*(-1.*Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z)*(Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)),3))/
(Sqrt(Power(Pattern(a,Blank(BH)),4) - 2*Power(Pattern(a,Blank(BH)),2)*
(Power(XX(x,y,z),2) + Power(YY(x,y,z),2) - Power(ZZ(x,y,z),2)) + 
Power(Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))*
Power(Power(Pattern(a,Blank(BH)),4) + Power(XX(x,y,z),4) + 
Power(YY(x,y,z),4) + 2.*Power(YY(x,y,z),2)*Power(ZZ(x,y,z),2) + 
Power(ZZ(x,y,z),4) + Power(Pattern(a,Blank(BH)),2)*(-2.*
Power(XX(x,y,z),2) - 2.*Power(YY(x,y,z),2) + 2.*Power(ZZ(x,y,z),2)) + 
Power(XX(x,y,z),2)*(2.*Power(YY(x,y,z),2) + 2.*Power(ZZ(x,y,z),2)),2))
;
}
KS_func_def_macro(dddR_D1D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
1.*Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),x,y)) + 1.*
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),x,z)) + 1.*Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y,z)) + 1.*Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),x,y)) + 
1.*Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),x,z)) + 1.*
Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),y,z)) + 1.*Hold(D(ZZ(x,y,z),z))*
Hold(D(ZZ(x,y,z),x,y)) + 1.*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,z)) + 
1.*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y,z)) + 1.*
Hold(D(XX(x,y,z),x,y,z))*XX(x,y,z) + 1.*Hold(D(YY(x,y,z),x,y,z))*
YY(x,y,z) + 1.*Hold(D(ZZ(x,y,z),x,y,z))*ZZ(x,y,z) + (0.5*(-12*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 6*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),2.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + 
Hold(D(YY(x,y,z),y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + 
Hold(D(XX(x,y,z),x,y))*XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),y))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 2*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(-4*
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) - 
4*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 4*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) - 2*
(Hold(D(XX(x,y,z),x))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),x,z))*XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(-4*
Hold(D(ZZ(x,y,z),x))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) - 2*
(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2)))*(4*
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z)) + 2*
(Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*
Hold(D(YY(x,y,z),z)) + Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + 
Hold(D(XX(x,y,z),y,z))*XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),y,z))*ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/
Power(4*Power(Pattern(a,Blank(BH)),2)*Power(ZZ(x,y,z),2) + Power(-
Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2),1.5) + (0.5*(4*
Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,y))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y,z))*Power(Pattern(a,Blank(BH)),2) + 
4*Hold(D(ZZ(x,y,z),x,y,z))*Power(Pattern(a,Blank(BH)),2)*ZZ(x,y,z) + 4*
(Hold(D(XX(x,y,z),z))*XX(x,y,z) + Hold(D(YY(x,y,z),z))*YY(x,y,z) + 
Hold(D(ZZ(x,y,z),z))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),y)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),y)) + Hold(D(XX(x,y,z),x,y))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,y))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,y))*
ZZ(x,y,z)) + 4*(Hold(D(XX(x,y,z),y))*XX(x,y,z) + Hold(D(YY(x,y,z),y))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),y))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),x))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),x))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),x,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),x,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,z))*
ZZ(x,y,z)) + 4*(Hold(D(XX(x,y,z),x))*XX(x,y,z) + Hold(D(YY(x,y,z),x))*
YY(x,y,z) + Hold(D(ZZ(x,y,z),x))*ZZ(x,y,z))*(Hold(D(XX(x,y,z),y))*
Hold(D(XX(x,y,z),z)) + Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),z)) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),z)) + Hold(D(XX(x,y,z),y,z))*
XX(x,y,z) + Hold(D(YY(x,y,z),y,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),y,z))*
ZZ(x,y,z)) + 2*(Hold(D(XX(x,y,z),z))*Hold(D(XX(x,y,z),x,y)) + 
Hold(D(XX(x,y,z),y))*Hold(D(XX(x,y,z),x,z)) + Hold(D(XX(x,y,z),x))*
Hold(D(XX(x,y,z),y,z)) + Hold(D(YY(x,y,z),z))*Hold(D(YY(x,y,z),x,y)) + 
Hold(D(YY(x,y,z),y))*Hold(D(YY(x,y,z),x,z)) + Hold(D(YY(x,y,z),x))*
Hold(D(YY(x,y,z),y,z)) + Hold(D(ZZ(x,y,z),z))*Hold(D(ZZ(x,y,z),x,y)) + 
Hold(D(ZZ(x,y,z),y))*Hold(D(ZZ(x,y,z),x,z)) + Hold(D(ZZ(x,y,z),x))*
Hold(D(ZZ(x,y,z),y,z)) + Hold(D(XX(x,y,z),x,y,z))*XX(x,y,z) + 
Hold(D(YY(x,y,z),x,y,z))*YY(x,y,z) + Hold(D(ZZ(x,y,z),x,y,z))*
ZZ(x,y,z))*(-Power(Pattern(a,Blank(BH)),2) + Power(XX(x,y,z),2) + 
Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2))))/Sqrt(4*Power(Pattern(a,Blank(BH)),2)*
Power(ZZ(x,y,z),2) + Power(-Power(Pattern(a,Blank(BH)),2) + 
Power(XX(x,y,z),2) + Power(YY(x,y,z),2) + Power(ZZ(x,y,z),2),2))
;
}

