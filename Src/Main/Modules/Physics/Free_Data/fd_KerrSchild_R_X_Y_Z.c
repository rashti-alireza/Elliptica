/* msimplify(X) */
/* msimplify(Y) */
/* msimplify(Z) */
#include "fd_KerrSchild_header.h"
KS_func_def_macro(X) KS_func_args_macro;
KS_func_def_macro(dX_D2) KS_func_args_macro;
KS_func_def_macro(dX_D0) KS_func_args_macro;
KS_func_def_macro(dX_D1) KS_func_args_macro;
KS_func_def_macro(ddX_D1D2) KS_func_args_macro;
KS_func_def_macro(ddX_D0D2) KS_func_args_macro;
KS_func_def_macro(ddX_D0D1) KS_func_args_macro;
KS_func_def_macro(ddX_D1D1) KS_func_args_macro;
KS_func_def_macro(ddX_D0D0) KS_func_args_macro;
KS_func_def_macro(ddX_D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddX_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddX_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddX_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddX_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddX_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddX_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddX_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddX_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddX_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddX_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddX_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddX_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddX_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddX_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddX_D0D1D0) KS_func_args_macro;
KS_func_def_macro(Y) KS_func_args_macro;
KS_func_def_macro(dY_D1) KS_func_args_macro;
KS_func_def_macro(dY_D0) KS_func_args_macro;
KS_func_def_macro(dY_D2) KS_func_args_macro;
KS_func_def_macro(ddY_D1D1) KS_func_args_macro;
KS_func_def_macro(ddY_D1D2) KS_func_args_macro;
KS_func_def_macro(ddY_D0D0) KS_func_args_macro;
KS_func_def_macro(ddY_D0D1) KS_func_args_macro;
KS_func_def_macro(ddY_D0D2) KS_func_args_macro;
KS_func_def_macro(ddY_D2D2) KS_func_args_macro;
KS_func_def_macro(dddY_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddY_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddY_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddY_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddY_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddY_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddY_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddY_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddY_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddY_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddY_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddY_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddY_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddY_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddY_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddY_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddY_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddY_D0D0D0) KS_func_args_macro;
KS_func_def_macro(Z) KS_func_args_macro;
KS_func_def_macro(dZ_D0) KS_func_args_macro;
KS_func_def_macro(dZ_D1) KS_func_args_macro;
KS_func_def_macro(dZ_D2) KS_func_args_macro;
KS_func_def_macro(ddZ_D0D1) KS_func_args_macro;
KS_func_def_macro(ddZ_D0D0) KS_func_args_macro;
KS_func_def_macro(ddZ_D0D2) KS_func_args_macro;
KS_func_def_macro(ddZ_D1D1) KS_func_args_macro;
KS_func_def_macro(ddZ_D1D2) KS_func_args_macro;
KS_func_def_macro(ddZ_D2D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddZ_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddZ_D2D2D0) KS_func_args_macro;
KS_func_def_macro(R) KS_func_args_macro;
KS_func_def_macro(dR_D0) KS_func_args_macro;
KS_func_def_macro(dR_D1) KS_func_args_macro;
KS_func_def_macro(dR_D2) KS_func_args_macro;
KS_func_def_macro(ddR_D1D1) KS_func_args_macro;
KS_func_def_macro(ddR_D1D2) KS_func_args_macro;
KS_func_def_macro(ddR_D2D2) KS_func_args_macro;
KS_func_def_macro(ddR_D0D1) KS_func_args_macro;
KS_func_def_macro(ddR_D0D0) KS_func_args_macro;
KS_func_def_macro(ddR_D0D2) KS_func_args_macro;
KS_func_def_macro(dddR_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddR_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddR_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddR_D0D0D1) KS_func_args_macro;
KS_func_def_macro(dddR_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddR_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddR_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddR_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddR_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddR_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddR_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddR_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddR_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddR_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddR_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddR_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddR_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddR_D1D1D0) KS_func_args_macro;
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
KS_func_def_macro(ddX_D1D2) KS_func_args_macro
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
KS_func_def_macro(ddX_D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddX_D1D1) KS_func_args_macro
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
KS_func_def_macro(ddX_D2D2) KS_func_args_macro
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
KS_func_def_macro(dddX_D0D2D2) KS_func_args_macro
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
KS_func_def_macro(dddX_D0D2D0) KS_func_args_macro
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
KS_func_def_macro(dddX_D1D1D0) KS_func_args_macro
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
KS_func_def_macro(dddX_D1D2D2) KS_func_args_macro
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
KS_func_def_macro(dddX_D0D1D1) KS_func_args_macro
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
KS_func_def_macro(dddX_D0D0D2) KS_func_args_macro
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
KS_func_def_macro(dddX_D0D0D0) KS_func_args_macro
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
KS_func_def_macro(dddX_D1D2D1) KS_func_args_macro
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
KS_func_def_macro(ddY_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddY_D1D2) KS_func_args_macro
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
KS_func_def_macro(ddY_D0D1) KS_func_args_macro
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
KS_func_def_macro(dddY_D1D2D2) KS_func_args_macro
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
KS_func_def_macro(dddY_D0D2D2) KS_func_args_macro
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
KS_func_def_macro(dddY_D1D2D1) KS_func_args_macro
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
KS_func_def_macro(dddY_D1D1D2) KS_func_args_macro
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
KS_func_def_macro(dddY_D0D1D1) KS_func_args_macro
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
KS_func_def_macro(dddY_D1D2D0) KS_func_args_macro
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
KS_func_def_macro(dddY_D0D0D1) KS_func_args_macro
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
KS_func_def_macro(ddZ_D0D1) KS_func_args_macro
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
KS_func_def_macro(ddZ_D0D2) KS_func_args_macro
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
KS_func_def_macro(ddZ_D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
0
;
}
KS_func_def_macro(ddZ_D2D2) KS_func_args_macro
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
KS_func_def_macro(dddZ_D0D1D2) KS_func_args_macro
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
KS_func_def_macro(dddZ_D0D1D0) KS_func_args_macro
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
KS_func_def_macro(dddZ_D0D2D0) KS_func_args_macro
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
KS_func_def_macro(dddZ_D0D2D2) KS_func_args_macro
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
KS_func_def_macro(dddZ_D1D1D0) KS_func_args_macro
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
KS_func_def_macro(dddZ_D1D1D1) KS_func_args_macro
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
KS_func_def_macro(dddZ_D0D0D0) KS_func_args_macro
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
KS_func_def_macro(dddZ_D1D2D1) KS_func_args_macro
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
KS_func_def_macro(R) KS_func_args_macro
{
return
/* mcode in progress ... */
0.707106781186548*Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5)
;
}
KS_func_def_macro(dR_D0) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5)
;
}
KS_func_def_macro(dR_D1) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5)
;
}
KS_func_def_macro(dR_D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5)
;
}
KS_func_def_macro(ddR_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D1_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D1_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D1_) )),2) + 1.*
Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D1_) )),2) + 2.*Power(Hold(((fd_ks_dY_D1_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) + 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(ddR_D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2) + 4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)) + (2.*
Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 
2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(ddR_D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D2_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D2_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D2_) )),2) + 1.*
Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D2_) )),2) + 2.*Power(Hold(((fd_ks_dY_D2_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) + 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(ddR_D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 1.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2) + 4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)) + (2.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 
2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(ddR_D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D0_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D0_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D0_) )),2) + 1.*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(ddR_D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2) + 4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)) + (2.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + 
2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D0D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D2_) ))*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D0D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D0_) )),2) + 
1.*Power(Hold(((fd_ks_dY_D0_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D0_) )),2) + 
1.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D2_) )) + Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 
1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D2_) ))*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D2_) ))*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dX_D0_) ))*
Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 
4.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 2.*
Hold((fd_ks_dddX_D0D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D0D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D0D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) - 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) - 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) - 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-12.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) - 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D2_) )) - 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) - 
1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + 
1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-3.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 3.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
3.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (1.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D0_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D0_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D0_) )),2) + 1.*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D2D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(2.121320343559644*Power(-1.*Hold(((fd_ks_dX_D2_) ))*Power(Pattern(a,BH),2)*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dX_D2_) ))*Power(fd_ks_X(x,y,z),3) - 1.*
Hold(((fd_ks_dY_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dY_D2_) ))*Power(fd_ks_X(x,y,z),2)*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z)*Power(fd_ks_Y(x,y,z),2) + 1.*
Hold(((fd_ks_dY_D2_) ))*Power(fd_ks_Y(x,y,z),3) + 1.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
Power(fd_ks_X(x,y,z),2)*fd_ks_Z(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
Power(fd_ks_Y(x,y,z),2)*fd_ks_Z(x,y,z) + 1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z)*
Power(fd_ks_Z(x,y,z),2) + 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z)*
Power(fd_ks_Z(x,y,z),2) + 1.*Hold(((fd_ks_dZ_D2_) ))*Power(fd_ks_Z(x,y,z),3) + 1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z)*Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z)*Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)*Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),3))/
(Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)*Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5)) + 
(0.707106781186548*(3.*Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro )) + 
3.*Hold(((fd_ks_dY_D2_) ))*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro )) + 3.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro )) + 1.*
Hold((fd_ks_dddX_D2D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D2D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D2D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (12.*
Power(Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),3))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5) + (0.5*(12.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D2D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*
(Power(Hold(((fd_ks_dX_D2_) )),2) + Power(Hold(((fd_ks_dY_D2_) )),2) + 
Power(Hold(((fd_ks_dZ_D2_) )),2) + Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 8.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D2_) )),2) + 1.*Power(Hold(((fd_ks_dY_D2_) )),2) + 
1.*Power(Hold(((fd_ks_dZ_D2_) )),2) + 1.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (6.*Hold(((fd_ks_dX_D2_) ))*
Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro )) + 6.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro )) + 6.*Hold(((fd_ks_dZ_D2_) ))*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro )) + 2.*Hold((fd_ks_dddX_D2D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D2D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_dddZ_D2D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (1.*(4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D2_) )),2) + 2.*Power(Hold(((fd_ks_dY_D2_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) + 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) - 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 
4.*Power(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (-2.*
Power(Hold(((fd_ks_dX_D2_) )),2) - 2.*Power(Hold(((fd_ks_dY_D2_) )),2) - 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) - 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) - 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(1.4142135623731*(1.*Power(Hold(((fd_ks_dX_D2_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D2_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D2_) )),2) + 1.*
Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D2_) )),2) + 2.*Power(Hold(((fd_ks_dY_D2_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) + 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Power(Hold(((fd_ks_dX_D2_) )),2) - 1.*
Power(Hold(((fd_ks_dY_D2_) )),2) - 1.*Power(Hold(((fd_ks_dZ_D2_) )),2) - 1.*
Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (2.*
Power(Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) - (0.5*(4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D2_) )),2) + 2.*Power(Hold(((fd_ks_dY_D2_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) + 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D0D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
(2.121320343559644*Power(-1.*Hold(((fd_ks_dX_D0_) ))*Power(Pattern(a,BH),2)*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dX_D0_) ))*Power(fd_ks_X(x,y,z),3) - 1.*
Hold(((fd_ks_dY_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*Power(fd_ks_X(x,y,z),2)*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z)*Power(fd_ks_Y(x,y,z),2) + 1.*
Hold(((fd_ks_dY_D0_) ))*Power(fd_ks_Y(x,y,z),3) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Power(fd_ks_X(x,y,z),2)*fd_ks_Z(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Power(fd_ks_Y(x,y,z),2)*fd_ks_Z(x,y,z) + 1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z)*
Power(fd_ks_Z(x,y,z),2) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z)*
Power(fd_ks_Z(x,y,z),2) + 1.*Hold(((fd_ks_dZ_D0_) ))*Power(fd_ks_Z(x,y,z),3) + 1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z)*Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + 
1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z)*Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z)*Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),3))/
(Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)*Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5)) + 
(0.707106781186548*(3.*Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro )) + 
3.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro )) + 3.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro )) + 1.*
Hold((fd_ks_dddX_D0D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (12.*
Power(Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),3))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5) + (0.5*(12.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*
(Power(Hold(((fd_ks_dX_D0_) )),2) + Power(Hold(((fd_ks_dY_D0_) )),2) + 
Power(Hold(((fd_ks_dZ_D0_) )),2) + Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 8.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D0_) )),2) + 1.*Power(Hold(((fd_ks_dY_D0_) )),2) + 
1.*Power(Hold(((fd_ks_dZ_D0_) )),2) + 1.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (6.*Hold(((fd_ks_dX_D0_) ))*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro )) + 6.*Hold(((fd_ks_dY_D0_) ))*
Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro )) + 6.*Hold(((fd_ks_dZ_D0_) ))*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro )) + 2.*Hold((fd_ks_dddX_D0D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D0D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_dddZ_D0D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (1.*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) - 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 
4.*Power(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (-2.*
Power(Hold(((fd_ks_dX_D0_) )),2) - 2.*Power(Hold(((fd_ks_dY_D0_) )),2) - 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) - 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) - 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(1.4142135623731*(1.*Power(Hold(((fd_ks_dX_D0_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D0_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D0_) )),2) + 1.*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Power(Hold(((fd_ks_dX_D0_) )),2) - 1.*
Power(Hold(((fd_ks_dY_D0_) )),2) - 1.*Power(Hold(((fd_ks_dZ_D0_) )),2) - 1.*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (2.*
Power(Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) - (0.5*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D0D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D0D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D0_) )),2) + 
1.*Power(Hold(((fd_ks_dY_D0_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D0_) )),2) + 
1.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D1_) )) + Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D1_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 
1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dX_D0_) ))*
Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 
4.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 2.*
Hold((fd_ks_dddX_D0D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D0D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D0D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) - 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) - 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) - 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-12.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) - 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D1_) )) - 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) - 
1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D1_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 
1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-3.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 3.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
3.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (1.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D0_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D0_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D0_) )),2) + 1.*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D1D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(-1.*Power(Hold(((fd_ks_dX_D2_) )),2) - 1.*
Power(Hold(((fd_ks_dY_D2_) )),2) - 1.*Power(Hold(((fd_ks_dZ_D2_) )),2) - 1.*
Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (2.*
Power(Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) - (0.5*(4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D2_) )),2) + 2.*Power(Hold(((fd_ks_dY_D2_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) + 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D1D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D1D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D1D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D1D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*
(Power(Hold(((fd_ks_dX_D2_) )),2) + Power(Hold(((fd_ks_dY_D2_) )),2) + 
Power(Hold(((fd_ks_dZ_D2_) )),2) + Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 8.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro )) + 
2.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro )) + 4.*
Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 
2.*Hold((fd_ks_dddX_D1D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D1D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D1D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(-4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) - 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 
4.*Power(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (-2.*
Power(Hold(((fd_ks_dX_D2_) )),2) - 2.*Power(Hold(((fd_ks_dY_D2_) )),2) - 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) - 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) - 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (1.*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-12.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(1.4142135623731*(1.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2) + 4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)) + (2.*
Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 
2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-3.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 3.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
3.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (1.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5)
;
}
KS_func_def_macro(dddR_D1D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D2_) ))*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D1D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D1D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D1D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D1D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D1_) )),2) + 
1.*Power(Hold(((fd_ks_dY_D1_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D1_) )),2) + 
1.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
Hold(((fd_ks_dX_D2_) )) + Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*
Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 
1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D2_) ))*
Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D2_) ))*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dX_D1_) ))*
Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 
4.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 2.*
Hold((fd_ks_dddX_D1D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D1D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D1D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) - 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) - 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) - 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-12.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D1_) )),2) + 2.*Power(Hold(((fd_ks_dY_D1_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) + 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) - 1.*Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) - 1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) - 
1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 
1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-3.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 3.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
3.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (1.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D1_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D1_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D1_) )),2) + 1.*
Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D1_) )),2) + 2.*Power(Hold(((fd_ks_dY_D1_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) + 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D0D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dZ_D1_) ))*
Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D1_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 
1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*
(Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) + Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 
Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 
Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D2_) ))*
Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 
2.*Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dY_D2_) ))*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 
2.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 2.*Hold((fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) - 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) - 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) - 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-12.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + 
1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) - 1.*Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) - 1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) - 
1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-3.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 3.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
3.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (1.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 1.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2) + 4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)) + (2.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 
2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D0D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(-1.*Power(Hold(((fd_ks_dX_D1_) )),2) - 1.*
Power(Hold(((fd_ks_dY_D1_) )),2) - 1.*Power(Hold(((fd_ks_dZ_D1_) )),2) - 1.*
Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (2.*
Power(Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) - (0.5*(4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D1_) )),2) + 2.*Power(Hold(((fd_ks_dY_D1_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) + 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D0D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*
(Power(Hold(((fd_ks_dX_D1_) )),2) + Power(Hold(((fd_ks_dY_D1_) )),2) + 
Power(Hold(((fd_ks_dZ_D1_) )),2) + Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 8.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 1.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro )) + 
2.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro )) + 4.*
Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 
2.*Hold((fd_ks_dddX_D0D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D0D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D0D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(-4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) - 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 
4.*Power(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (-2.*
Power(Hold(((fd_ks_dX_D1_) )),2) - 2.*Power(Hold(((fd_ks_dY_D1_) )),2) - 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) - 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) - 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (1.*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-12.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(1.4142135623731*(1.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 1.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2) + 4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)) + (2.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 
2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-3.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 3.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
3.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (1.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5)
;
}
KS_func_def_macro(dddR_D0D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D1_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 
1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D0D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D0_) )),2) + 
1.*Power(Hold(((fd_ks_dY_D0_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D0_) )),2) + 
1.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D1_) )) + Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D1_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 
1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dX_D0_) ))*
Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 
4.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 2.*
Hold((fd_ks_dddX_D0D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D0D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D0D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) - 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) - 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) - 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-12.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D0_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D0_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D0_) )),2) + 1.*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-3.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 3.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 3.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(1.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) - 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D1_) )) - 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) - 
1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D2D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D2_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D2_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D2_) )),2) + 1.*
Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D2_) )),2) + 2.*Power(Hold(((fd_ks_dY_D2_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) + 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-3.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 3.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 3.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(1.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) - 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D2_) )) - 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) - 
1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + 
1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D0D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D2_) )),2) + 
1.*Power(Hold(((fd_ks_dY_D2_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D2_) )),2) + 
1.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D2_) )) + Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 
1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*
Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dX_D2_) ))*
Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D2_) ))*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 
4.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 2.*
Hold((fd_ks_dddX_D0D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D0D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D0D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D2_) )),2) + 2.*Power(Hold(((fd_ks_dY_D2_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) + 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) - 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) - 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) - 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-12.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5)
;
}
KS_func_def_macro(dddR_D0D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + 
1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D2_) ))*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D0D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D0_) )),2) + 
1.*Power(Hold(((fd_ks_dY_D0_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D0_) )),2) + 
1.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D2_) )) + Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 
1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D2_) ))*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D2_) ))*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dX_D0_) ))*
Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 
4.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 2.*
Hold((fd_ks_dddX_D0D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D0D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D0D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) - 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) - 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) - 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-12.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D0_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D0_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D0_) )),2) + 1.*
Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D0_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D0_) )),2) + 2.*Power(Hold(((fd_ks_dY_D0_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D0_) )),2) + 2.*Hold((fd_ks_ddX_D0D0 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D0 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D0D0 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-3.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 3.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 3.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(1.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) - 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D2_) )) - 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) - 
1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D0D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + 
1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dZ_D1_) ))*
Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D1_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 
1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*
(Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) + Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 
Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 
Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D2_) ))*
Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 
2.*Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dY_D2_) ))*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 
2.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 2.*Hold((fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) - 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) - 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) - 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-12.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 1.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2) + 4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)) + (2.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 
2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-3.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 3.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
3.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (1.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) - 1.*Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) - 1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) - 
1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D0D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(-1.*Power(Hold(((fd_ks_dX_D2_) )),2) - 1.*
Power(Hold(((fd_ks_dY_D2_) )),2) - 1.*Power(Hold(((fd_ks_dZ_D2_) )),2) - 1.*
Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (2.*
Power(Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) - (0.5*(4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D2_) )),2) + 2.*Power(Hold(((fd_ks_dY_D2_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) + 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D0D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*
(Power(Hold(((fd_ks_dX_D2_) )),2) + Power(Hold(((fd_ks_dY_D2_) )),2) + 
Power(Hold(((fd_ks_dZ_D2_) )),2) + Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 8.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro )) + 
2.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro )) + 4.*
Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 
2.*Hold((fd_ks_dddX_D0D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D0D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D0D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(-4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) - 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 
4.*Power(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (-2.*
Power(Hold(((fd_ks_dX_D2_) )),2) - 2.*Power(Hold(((fd_ks_dY_D2_) )),2) - 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) - 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) - 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (1.*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-12.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(1.4142135623731*(1.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2) + 4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)) + (2.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + 
2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-3.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 3.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
3.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (1.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5)
;
}
KS_func_def_macro(dddR_D2D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D2_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D2_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D2_) )),2) + 1.*
Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D2_) )),2) + 2.*Power(Hold(((fd_ks_dY_D2_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) + 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-3.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 3.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 3.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - 
(1.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) - 1.*Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) - 1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) - 
1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 
1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D1D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D1D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D1D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D1D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D2_) )),2) + 
1.*Power(Hold(((fd_ks_dY_D2_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D2_) )),2) + 
1.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
Hold(((fd_ks_dX_D2_) )) + Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*
Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 
1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*
Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dX_D2_) ))*
Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D2_) ))*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 
4.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 2.*
Hold((fd_ks_dddX_D1D2D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D1D2D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D1D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D2_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D2_) )),2) + 2.*Power(Hold(((fd_ks_dY_D2_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D2_) )),2) + 2.*Hold((fd_ks_ddX_D2D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D2D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D2D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) - 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) - 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) - 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-12.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5)
;
}
KS_func_def_macro(dddR_D1D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 
1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dZ_D1_) ))*
Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 1.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D1_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 
1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(1.*Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D2_) )) + Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) + Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*
Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 
1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D2_) ))*
Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro )) + 
2.*Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dY_D2_) ))*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 
2.*Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 2.*Hold((fd_ks_dddX_D0D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_dddY_D0D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D0D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) - 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) - 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) - 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-12.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(-3.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 3.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 3.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(1.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 1.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2) + 4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)) + 
(2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) - 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D2_) )) - 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) )) - 
1.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D0D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D1D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) + 1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 
1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D2_) ))*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D2_) ))*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D1D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D1D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D1D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D2_) ))*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D1D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D1_) )),2) + 
1.*Power(Hold(((fd_ks_dY_D1_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D1_) )),2) + 
1.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D1_) ))*
Hold(((fd_ks_dX_D2_) )) + Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*
Hold(((fd_ks_dX_D2_) )) + 1.*Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 
1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) + 1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D2_) ))*
Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D2_) ))*
Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D2_) ))*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dX_D1_) ))*
Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro )) + 
4.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro )) + 2.*
Hold((fd_ks_dddX_D1D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D1D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D1D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) - 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) - 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) - 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D1_) )),2) + 2.*Power(Hold(((fd_ks_dY_D1_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) + 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-12.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D1_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D1_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D1_) )),2) + 1.*
Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D1_) )),2) + 2.*Power(Hold(((fd_ks_dY_D1_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) + 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-3.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 3.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 3.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - 
(1.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D2_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D2_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) - 1.*Hold(((fd_ks_dY_D1_) ))*
Hold(((fd_ks_dY_D2_) )) - 1.*Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) )) - 
1.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D1_) ))*Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D2_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D1_) ))*Hold(((fd_ks_dX_D2_) )) + 2.*
Hold(((fd_ks_dY_D1_) ))*Hold(((fd_ks_dY_D2_) )) + 2.*Hold(((fd_ks_dZ_D1_) ))*
Hold(((fd_ks_dZ_D2_) )) + 2.*Hold((fd_ks_ddX_D1D2 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D1D2 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D1D2 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D2_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D2_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D2_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D2_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D1D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
(2.121320343559644*Power(-1.*Hold(((fd_ks_dX_D1_) ))*Power(Pattern(a,BH),2)*
fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dX_D1_) ))*Power(fd_ks_X(x,y,z),3) - 1.*
Hold(((fd_ks_dY_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dY_D1_) ))*Power(fd_ks_X(x,y,z),2)*fd_ks_Y(x,y,z) + 1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z)*Power(fd_ks_Y(x,y,z),2) + 1.*
Hold(((fd_ks_dY_D1_) ))*Power(fd_ks_Y(x,y,z),3) + 1.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*
Power(fd_ks_X(x,y,z),2)*fd_ks_Z(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*
Power(fd_ks_Y(x,y,z),2)*fd_ks_Z(x,y,z) + 1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z)*
Power(fd_ks_Z(x,y,z),2) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z)*
Power(fd_ks_Z(x,y,z),2) + 1.*Hold(((fd_ks_dZ_D1_) ))*Power(fd_ks_Z(x,y,z),3) + 1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z)*Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + 
1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z)*Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)*Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),3))/
(Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)*Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5)) + 
(0.707106781186548*(3.*Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro )) + 
3.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro )) + 3.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro )) + 1.*
Hold((fd_ks_dddX_D1D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D1D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D1D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (12.*
Power(Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),3))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5) + (0.5*(12.*
Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D1D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*
(Power(Hold(((fd_ks_dX_D1_) )),2) + Power(Hold(((fd_ks_dY_D1_) )),2) + 
Power(Hold(((fd_ks_dZ_D1_) )),2) + Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 8.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D1_) )),2) + 1.*Power(Hold(((fd_ks_dY_D1_) )),2) + 
1.*Power(Hold(((fd_ks_dZ_D1_) )),2) + 1.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (6.*Hold(((fd_ks_dX_D1_) ))*
Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro )) + 6.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro )) + 6.*Hold(((fd_ks_dZ_D1_) ))*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro )) + 2.*Hold((fd_ks_dddX_D1D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D1D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_dddZ_D1D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (1.*(4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D1_) )),2) + 2.*Power(Hold(((fd_ks_dY_D1_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) + 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) - 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 
4.*Power(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (-2.*
Power(Hold(((fd_ks_dX_D1_) )),2) - 2.*Power(Hold(((fd_ks_dY_D1_) )),2) - 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) - 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) - 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5) + 
(1.4142135623731*(1.*Power(Hold(((fd_ks_dX_D1_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D1_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D1_) )),2) + 1.*
Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D1_) )),2) + 2.*Power(Hold(((fd_ks_dY_D1_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) + 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Power(Hold(((fd_ks_dX_D1_) )),2) - 1.*
Power(Hold(((fd_ks_dY_D1_) )),2) - 1.*Power(Hold(((fd_ks_dZ_D1_) )),2) - 1.*
Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (2.*
Power(Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) - (0.5*(4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D1_) )),2) + 2.*Power(Hold(((fd_ks_dY_D1_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) + 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5)
;
}
KS_func_def_macro(dddR_D1D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
(0.707106781186548*(1.*Power(Hold(((fd_ks_dX_D1_) )),2) + 1.*
Power(Hold(((fd_ks_dY_D1_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D1_) )),2) + 1.*
Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (2.*
Power(Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z)*(-1.*Power(Pattern(a,BH),2) - 
1.*Power(fd_ks_X(x,y,z),2) - 1.*Power(fd_ks_Y(x,y,z),2) - 1.*
Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z)*(1.*
Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 1.*
Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z)*(1.*Power(Pattern(a,BH),2) - 1.*Power(fd_ks_X(x,y,z),2) - 
1.*Power(fd_ks_Y(x,y,z),2) - 1.*Power(fd_ks_Z(x,y,z),2)),2))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D1_) )),2) + 2.*Power(Hold(((fd_ks_dY_D1_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) + 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D0_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D0_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-3.*Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) - 3.*
Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) - 3.*Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z) - 
(1.5*(4.*Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + (0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),2.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) + 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(-1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) - 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D1_) )) - 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) - 
1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) - 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) - (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) - (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(-1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) - 1.*
Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) - 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z) - 
(0.5*(4.*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
2.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5))*(1.*
Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*
Hold(((fd_ks_dY_D1_) )) + 1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 
1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),1.5) + 
(0.707106781186548*(1.*Hold(((fd_ks_dX_D0_) ))*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro )) + 
1.*Hold(((fd_ks_dY_D0_) ))*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro )) + 1.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro )) + 2.*
Hold(((fd_ks_dX_D1_) ))*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D1_) ))*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 
1.*Hold((fd_ks_dddX_D0D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_dddY_D0D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_dddZ_D0D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
8.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_dddZ_D0D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*
fd_ks_Z(x,y,z) + 4.*(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Power(Hold(((fd_ks_dX_D1_) )),2) + 
1.*Power(Hold(((fd_ks_dY_D1_) )),2) + 1.*Power(Hold(((fd_ks_dZ_D1_) )),2) + 
1.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 1.*Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + 4.*(1.*
Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D1_) )) + Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z)) + 4.*(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D0_) ))*
Hold(((fd_ks_dX_D1_) )) + 1.*Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 
1.*Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) )) + 1.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 1.*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 1.*
Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*
Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dY_D0_) ))*
Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dX_D1_) ))*
Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro )) + 4.*Hold(((fd_ks_dY_D1_) ))*Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro )) + 
4.*Hold(((fd_ks_dZ_D1_) ))*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro )) + 2.*
Hold((fd_ks_dddX_D0D1D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*Hold((fd_ks_dddY_D0D1D1 KS_func_pass_args_macro ))*
fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_dddZ_D0D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5) + (0.5*(4.*
Power(Hold(((fd_ks_dZ_D1_) )),2)*Power(Pattern(a,BH),2) + 4.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 
4.*Power(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 1.*Hold(((fd_ks_dY_D1_) ))*
fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z),2) + (2.*
Power(Hold(((fd_ks_dX_D1_) )),2) + 2.*Power(Hold(((fd_ks_dY_D1_) )),2) + 2.*
Power(Hold(((fd_ks_dZ_D1_) )),2) + 2.*Hold((fd_ks_ddX_D1D1 KS_func_pass_args_macro ))*
fd_ks_X(x,y,z) + 2.*Hold((fd_ks_ddY_D1D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*
Hold((fd_ks_ddZ_D1D1 KS_func_pass_args_macro ))*fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) + 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) + 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) + 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) + 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) + 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) + 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-4.*
Hold(((fd_ks_dZ_D0_) ))*Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2) - 
4.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 4.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(1.*Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + 
1.*Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 1.*Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z)) + (-2.*Hold(((fd_ks_dX_D0_) ))*Hold(((fd_ks_dX_D1_) )) - 2.*
Hold(((fd_ks_dY_D0_) ))*Hold(((fd_ks_dY_D1_) )) - 2.*Hold(((fd_ks_dZ_D0_) ))*
Hold(((fd_ks_dZ_D1_) )) - 2.*Hold((fd_ks_ddX_D0D1 KS_func_pass_args_macro ))*fd_ks_X(x,y,z) - 2.*
Hold((fd_ks_ddY_D0D1 KS_func_pass_args_macro ))*fd_ks_Y(x,y,z) - 2.*Hold((fd_ks_ddZ_D0D1 KS_func_pass_args_macro ))*
fd_ks_Z(x,y,z))*(-Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*Hold(((fd_ks_dZ_D1_) ))*
Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*(Hold(((fd_ks_dX_D1_) ))*
fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + Hold(((fd_ks_dZ_D1_) ))*
fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/Power(4*
Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),1.5) + (0.5*(-12.*
Hold(((fd_ks_dZ_D0_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 6.*
(Hold(((fd_ks_dX_D0_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D0_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D0_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(-4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) - 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2)))*(4.*
Hold(((fd_ks_dZ_D1_) ))*Power(Pattern(a,BH),2)*fd_ks_Z(x,y,z) + 2.*
(Hold(((fd_ks_dX_D1_) ))*fd_ks_X(x,y,z) + Hold(((fd_ks_dY_D1_) ))*fd_ks_Y(x,y,z) + 
Hold(((fd_ks_dZ_D1_) ))*fd_ks_Z(x,y,z))*(-1.*Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2))))/
Power(4*Power(Pattern(a,BH),2)*Power(fd_ks_Z(x,y,z),2) + Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),2.5)))/Power(-
Power(Pattern(a,BH),2) + Power(fd_ks_X(x,y,z),2) + 
Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2) + Power(4*Power(Pattern(a,BH),2)*
Power(fd_ks_Z(x,y,z),2) + Power(-Power(Pattern(a,BH),2) + 
Power(fd_ks_X(x,y,z),2) + Power(fd_ks_Y(x,y,z),2) + Power(fd_ks_Z(x,y,z),2),2),0.5),0.5)
;
}

