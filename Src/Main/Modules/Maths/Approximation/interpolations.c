/*
// Alireza Rashti
// July 2018
*/

#include "interpolations.h"

/* initilizing an Interpolation_T struct with calloc */
Interpolation_T *init_interpolation(void)
{
  Interpolation_T *s = calloc(1,sizeof(*s));
  pointerEr(s);
  
  s->point = calloc(1,sizeof(*s->point));
  pointerEr(s->point);
  
  return s;
}

/* interpolation function */
double interpolation(Interpolation_T *const interp_s)
{
  UNUSED(interp_s);
  return 0;
}

/*
Interpolation_Func_T *WhichInterpolation(const Patch_T *const patch)
{
   if (strstr_i(GetParameterS("Interpolation_Method"),"Spectral"))
    intp = WhichSpectralInterpolation(f);
  else
    abortEr(INCOMPLETE_FUNC);
 
}
*/

/* given field, based on the f->point find the value of interpolation.
// NOTE: this assumes field is "3-D". for 2-D field refer to "interpolation2"
// ->retutn value: interpolation value.
*/
/*
double interpolation(const Field_T *const f)
{
#error "this function is so slow each time needs to check many ifs! "\
          "its better to come up with better ideas"
  double intp = 0;
  
  if (strstr_i(GetParameterS("Interpolation_Method"),"Spectral"))
    intp = spectral_interpolation(f);
  else
    abortEr(INCOMPLETE_FUNC);
  
  return intp;
}
*/
/* finding interpolation for spectral method for interpolation, 
// based on given f->point.
// ->return value: interpolation value
*/
/*
static double spectral_interpolation(const Field_T *const f)
{
  Interpolation_Func_T *interpolation_f = WhichSpectralInterpolation(f);
  
  return interpolation_f(f);
}
*/
/* finding the best fit spectral interpolation function 
// based on collocation, basis of patch and whether the interested point
// is on plane or in volume. if it is on a plane one can take advantage of
// 2-D interpolation then. 
// ->return value: function for interpolation.
*/
/*
static Interpolation_Func_T *WhichSpectralInterpolation(const Field_T *const f)
{
  
}
*/
