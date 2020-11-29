#include "frda_ks_free_data_analytic.h"
KS_func_def_macro(rolloff) KS_func_args_macro;
KS_func_def_macro(drolloff_D0) KS_func_args_macro;
KS_func_def_macro(drolloff_D2) KS_func_args_macro;
KS_func_def_macro(drolloff_D1) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D2D2) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D1D2) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D0D1) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D0D2) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D0D0) KS_func_args_macro;
KS_func_def_macro(ddrolloff_D1D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D1D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D1D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D2D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D2D2D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D0D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D1D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D2D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D1D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D2D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D2D2D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D2D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D2D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D2D2D0) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D2D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D1D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D1D1D1) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D0D2) KS_func_args_macro;
KS_func_def_macro(dddrolloff_D0D0D1) KS_func_args_macro;
KS_func_def_macro(rolloff) KS_func_args_macro
{
return
/* mcode in progress ... */
exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/pow(r0, 4));
}
KS_func_def_macro(drolloff_D0) KS_func_args_macro
{
return
/* mcode in progress ... */
-4.0*x*pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.0)*exp(-pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 2.0)/pow(r0, 4))/pow(r0, 4);
}
KS_func_def_macro(drolloff_D2) KS_func_args_macro
{
return
/* mcode in progress ... */
-4.0*z*pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.0)*exp(-pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 2.0)/pow(r0, 4))/pow(r0, 4);
}
KS_func_def_macro(drolloff_D1) KS_func_args_macro
{
return
/* mcode in progress ... */
-4.0*y*pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.0)*exp(-pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 2.0)/pow(r0, 4))/pow(r0, 4);
}
KS_func_def_macro(ddrolloff_D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
(-pow(r0, 4)*(8.0*pow(z, 2) + 4.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 1.0)) + 16.0*pow(z, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 8);
}
KS_func_def_macro(ddrolloff_D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
y*z*(-8.0*pow(r0, 4) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 8);
}
KS_func_def_macro(ddrolloff_D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
x*y*(-8.0*pow(r0, 4) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 8);
}
KS_func_def_macro(ddrolloff_D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
x*z*(-8.0*pow(r0, 4) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 8);
}
KS_func_def_macro(ddrolloff_D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
(-pow(r0, 4)*(8.0*pow(x, 2) + 4.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 1.0)) + 16.0*pow(x, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 8);
}
KS_func_def_macro(ddrolloff_D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
(-pow(r0, 4)*(8.0*pow(y, 2) + 4.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 1.0)) + 16.0*pow(y, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 8);
}
KS_func_def_macro(dddrolloff_D0D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
y*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(x, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(x, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D1D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
z*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(y, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(y, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D0D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
z*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(x, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(x, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D2D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
y*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(z, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(z, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D0D0D0) KS_func_args_macro
{
return
/* mcode in progress ... */
x*(-24.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(x, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 48.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(x, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D0D1D2) KS_func_args_macro
{
return
/* mcode in progress ... */
x*y*z*(96.0*pow(r0, 4)*pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.0) - 
64.0*pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 2.0)/pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D0D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
x*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(z, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(z, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D1D1D0) KS_func_args_macro
{
return
/* mcode in progress ... */
x*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(y, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(y, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D1D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
z*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(y, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(y, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D2D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
z*(-24.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(z, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 48.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(z, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D1D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
x*y*z*(96.0*pow(r0, 4)*pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.0) - 
64.0*pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 2.0)/pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D1D2D2) KS_func_args_macro
{
return
/* mcode in progress ... */
y*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(z, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(z, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D2D2D0) KS_func_args_macro
{
return
/* mcode in progress ... */
x*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(z, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(z, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D0D2D1) KS_func_args_macro
{
return
/* mcode in progress ... */
x*y*z*(96.0*pow(r0, 4)*pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.0) - 
64.0*pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 2.0)/pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D0D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
x*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(y, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(y, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D1D1D1) KS_func_args_macro
{
return
/* mcode in progress ... */
y*(-24.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(y, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 48.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(y, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D0D0D2) KS_func_args_macro
{
return
/* mcode in progress ... */
z*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(x, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(x, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}
KS_func_def_macro(dddrolloff_D0D0D1) KS_func_args_macro
{
return
/* mcode in progress ... */
y*(-8.0*pow(r0, 8) + pow(r0, 4)*(96.0*pow(x, 2)*pow(pow(x, 2) + 
pow(y, 2) + pow(z, 2), 1.0) + 16.0*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 2.0)) - 64.0*pow(x, 2)*pow(pow(x, 2) + pow(y, 2) + 
pow(z, 2), 3.0))*exp(-pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.0)/
pow(r0, 4))/pow(r0, 12);
}

