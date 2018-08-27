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
  Interpolation_Func_T *const interp_func = 
    ChooseInterpolationFunction(interp_s);
  
  return interp_func(interp_s);
}

/* choosing interpolation function based on input file and patch properties.
// ->return value: selected interpolation method function.
*/
static Interpolation_Func_T *ChooseInterpolationFunction(const Interpolation_T *const interp_s)
{
  Interpolation_Func_T *func = 0;
  unsigned i;
  Patch_T *const patch = interp_s->point->patch;
  
  if (strstr_i(GetParameterS("Interpolation_Method"),"Spectral"))
  {
    for (i = 0; i < 3; ++i)
    {
      if (patch->basis[i] == Chebyshev_Tn_BASIS 
          && patch->collocation[i] == Chebyshev_Extrema)
      {
        func = interpolation_Chebyshev_Tn;
      }
      else
         abortEr(INCOMPLETE_FUNC);

    }/* end of for (i = 0; i < 3; ++i) */
  }
  else
    abortEr(INCOMPLETE_FUNC);
    
  return func;
}

/* calculating interpolation using Chebyshev spectral method and
// Chebyshev Extrema collocation points.
// ->return value: interpolation value
*/
static double interpolation_Chebyshev_Tn(Interpolation_T *const interp_s)
{
  double interp_v = 0;
  
  if (interp_s->point->sameX && interp_s->point->sameY && interp_s->point->sameZ)
  {
    abortEr("No interpolation is needed,"
      " since the points are collocated.\n");
  }
  else if (interp_s->point->sameX && interp_s->point->sameY)
  {
    interp_v = interpolation_Chebyshev_Tn_Z(interp_s);
  }
  else if (interp_s->point->sameX && interp_s->point->sameZ)
  {
    interp_v = interpolation_Chebyshev_Tn_Y(interp_s);
  }
  else if (interp_s->point->sameY && interp_s->point->sameZ)
  {
    interp_v = interpolation_Chebyshev_Tn_X(interp_s);
  }
  else if (interp_s->point->sameX)
  {
    interp_v = interpolation_Chebyshev_Tn_YZ(interp_s);
  }
  else if (interp_s->point->sameY)
  {
    interp_v = interpolation_Chebyshev_Tn_XZ(interp_s);
  }
  else if (interp_s->point->sameZ)
  {
    interp_v = interpolation_Chebyshev_Tn_XY(interp_s);
  }
  else/* in all directions */
  {
    interp_v = interpolation_Chebyshev_Tn_XYZ(interp_s);
  }
  
  return interp_v;
}

/* interpolation in X direction keeping the other two directions constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_X(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double X = interp_s->point->x[0];/* X coord of the interesting point */
  const unsigned *const n = interp_s->point->patch->n;
  const unsigned N = n[0];/* number of points in the specified direction */
  const unsigned ind = interp_s->point->ind;/* index of interesting point */
  const double *coeffs = 0;/* coeffs of expansion of field in Tn */
  double interp_v = 0;
  unsigned x,y,z;
  
  make_coeffs_1d(field,0);
  coeffs = field->v2;
  
  IJK(ind,n,&x,&y,&z);
  
  for (x = 1; x < N-1; ++x)
  {
    interp_v += coeffs[L(n,x,y,z)]*Cheb_Tn((int)x,X);
  }
  interp_v *= 2;
  
  interp_v += coeffs[L(n,0,y,z)] + coeffs[L(n,N-1,y,z)]*Cheb_Tn((int)N-1,X);
  
  return interp_v;
}

/* interpolation in Y direction keeping the other two directions constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_Y(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double Y = interp_s->point->x[1];/* Y coord of the interesting point */
  const unsigned *const n = interp_s->point->patch->n;
  const unsigned N = n[1];/* number of points in the specified direction */
  const unsigned ind = interp_s->point->ind;/* index of interesting point */
  const double *coeffs = 0;/* coeffs of expansion of field in Tn */
  double interp_v = 0;
  unsigned x,y,z;
  
  make_coeffs_1d(field,1);
  coeffs = field->v2;
  
  IJK(ind,n,&x,&y,&z);
  
  for (y = 1; y < N-1; ++y)
  {
    interp_v += coeffs[L(n,x,y,z)]*Cheb_Tn((int)y,Y);
  }
  interp_v *= 2;
  
  interp_v += coeffs[L(n,x,0,z)] + coeffs[L(n,x,N-1,z)]*Cheb_Tn((int)N-1,Y);
  
  return interp_v;
}

/* interpolation in Z direction keeping the other two directions constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_Z(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double Z = interp_s->point->x[2];/* Z coord of the interesting point */
  const unsigned *const n = interp_s->point->patch->n;
  const unsigned N = n[2];/* number of points in the specified direction */
  const unsigned ind = interp_s->point->ind;/* index of interesting point */
  const double *coeffs = 0;/* coeffs of expansion of field in Tn */
  double interp_v = 0;
  unsigned x,y,z;
  
  make_coeffs_1d(field,2);
  coeffs = field->v2;
  
  IJK(ind,n,&x,&y,&z);
  
  for (z = 1; z < N-1; ++z)
  {
    interp_v += coeffs[L(n,x,y,z)]*Cheb_Tn((int)z,Z);
  }
  interp_v *= 2;
  
  interp_v += coeffs[L(n,x,y,0)] + coeffs[L(n,x,y,N-1)]*Cheb_Tn((int)N-1,Z);
  
  return interp_v;
}

/* interpolation in X&Y directions keeping the other direction constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_XY(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double X = interp_s->point->x[0];/* X coord of the interesting point */
  const double Y = interp_s->point->x[1];/* Y coord of the interesting point */
  const unsigned *const n = interp_s->point->patch->n;
  const unsigned Nx = n[0];/* number of points in the specified direction */
  const unsigned Ny = n[1];/* number of points in the specified direction */
  const unsigned ind = interp_s->point->ind;/* index of interesting point */
  const double *coeffs = 0;/* coeffs of expansion of field in Tn */
  enum MAX_E {MAX = 5};
  double sum[MAX];
  double interp_v = 0;
  unsigned x,y,z,i;
  
  make_coeffs_2d(field,0,1);
  coeffs = field->v2;
  
  IJK(ind,n,&x,&y,&z);
  
  for (i = 0; i < MAX; ++i)
    sum[i] = 0;
  
  sum[0] = coeffs[L(n,0,0,z)]+coeffs[L(n,0,Ny-1,z)]*Cheb_Tn((int)Ny-1,Y)+
                coeffs[L(n,Nx-1,0,z)]*Cheb_Tn((int)Nx-1,X) 			+
                coeffs[L(n,Nx-1,Ny-1,z)]*Cheb_Tn((int)Nx-1,X)*Cheb_Tn((int)Ny-1,Y);
  
  for (y = 1; y < Ny-1; ++y)
  {
    sum[1] += coeffs[L(n,0,y,z)]*Cheb_Tn((int)y,Y);
  }
  
  for (x = 1; x < Nx-1; ++x)
  {
    sum[1] += coeffs[L(n,x,0,z)]*Cheb_Tn((int)x,X);
  }
  sum[1] *= 2;
  
  for (y = 1; y < Ny-1; ++y)
  {
    sum[2] += coeffs[L(n,Nx-1,y,z)]*Cheb_Tn((int)y,Y);
  }
  sum[2] *= 2*Cheb_Tn((int)Nx-1,X);
  
  for (x = 1; x < Nx-1; ++x)
  {
    sum[3] += coeffs[L(n,x,Ny-1,z)]*Cheb_Tn((int)x,X);
  }
  sum[3] *= 2*Cheb_Tn((int)Ny-1,Y);
  
  for (x = 0; x < Nx-1; ++x)
  {
    for (y = 0; y < Ny-1; ++y)
      sum[4] += coeffs[L(n,x,Ny-1,z)]*Cheb_Tn((int)x,X)*Cheb_Tn((int)y,Y);
  }
  sum[4] *= 4;
  
  for (i = 0; i < MAX; ++i)
    interp_v += sum[i];
    
  return interp_v;
}

/* interpolation in X&Z directions keeping the other direction constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_XZ(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double X = interp_s->point->x[0];/* X coord of the interesting point */
  const double Z = interp_s->point->x[2];/* Z coord of the interesting point */
  const unsigned *const n = interp_s->point->patch->n;
  const unsigned Nx = n[0];/* number of points in the specified direction */
  const unsigned Nz = n[2];/* number of points in the specified direction */
  const unsigned ind = interp_s->point->ind;/* index of interesting point */
  const double *coeffs = 0;/* coeffs of expansion of field in Tn */
  enum MAX_E {MAX = 5};
  double sum[MAX];
  double interp_v = 0;
  unsigned x,y,z,i;
  
  make_coeffs_2d(field,0,2);
  coeffs = field->v2;
  
  IJK(ind,n,&x,&y,&z);
  
  for (i = 0; i < MAX; ++i)
    sum[i] = 0;
  
  sum[0] = coeffs[L(n,0,y,0)]+coeffs[L(n,0,y,Nz-1)]*Cheb_Tn((int)Nz-1,Z)+
                coeffs[L(n,Nx-1,y,0)]*Cheb_Tn((int)Nx-1,X) 			+
                coeffs[L(n,Nx-1,y,Nz-1)]*Cheb_Tn((int)Nx-1,X)*Cheb_Tn((int)Nz-1,Z);
  
  for (z = 1; z < Nz-1; ++z)
  {
    sum[1] += coeffs[L(n,0,y,z)]*Cheb_Tn((int)z,Z);
  }
  
  for (x = 1; x < Nx-1; ++x)
  {
    sum[1] += coeffs[L(n,x,y,0)]*Cheb_Tn((int)x,X);
  }
  sum[1] *= 2;
  
  for (z = 1; z < Nz-1; ++z)
  {
    sum[2] += coeffs[L(n,Nx-1,y,z)]*Cheb_Tn((int)z,Z);
  }
  sum[2] *= 2*Cheb_Tn((int)Nx-1,X);
  
  for (x = 1; x < Nx-1; ++x)
  {
    sum[3] += coeffs[L(n,x,y,Nz-1)]*Cheb_Tn((int)x,X);
  }
  sum[3] *= 2*Cheb_Tn((int)Nz-1,Z);
  
  for (x = 0; x < Nx-1; ++x)
  {
    for (z = 0; z < Nz-1; ++z)
      sum[4] += coeffs[L(n,x,y,Nz-1)]*Cheb_Tn((int)x,X)*Cheb_Tn((int)z,Z);
  }
  sum[4] *= 4;
  
  for (i = 0; i < MAX; ++i)
    interp_v += sum[i];
    
  return interp_v;
}

/* interpolation in Y&Z direction keeping the other direction constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_YZ(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double Y = interp_s->point->x[1];/* Y coord of the interesting point */
  const double Z = interp_s->point->x[2];/* Z coord of the interesting point */
  const unsigned *const n = interp_s->point->patch->n;
  const unsigned Ny = n[1];/* number of points in the specified direction */
  const unsigned Nz = n[2];/* number of points in the specified direction */
  const unsigned ind = interp_s->point->ind;/* index of interesting point */
  const double *coeffs = 0;/* coeffs of expansion of field in Tn */
  enum MAX_E {MAX = 5};
  double sum[MAX];
  double interp_v = 0;
  unsigned x,y,z,i;
  
  make_coeffs_2d(field,1,2);
  coeffs = field->v2;
  
  IJK(ind,n,&x,&y,&z);
  
  for (i = 0; i < MAX; ++i)
    sum[i] = 0;
  
  sum[0] = coeffs[L(n,x,0,0)]+coeffs[L(n,x,0,Nz-1)]*Cheb_Tn((int)Nz-1,Z)+
                coeffs[L(n,x,Ny-1,0)]*Cheb_Tn((int)Ny-1,Y) 			+
                coeffs[L(n,x,Ny-1,Nz-1)]*Cheb_Tn((int)Ny-1,Y)*Cheb_Tn((int)Nz-1,Z);
  
  for (z = 1; z < Nz-1; ++z)
  {
    sum[1] += coeffs[L(n,x,0,z)]*Cheb_Tn((int)z,Z);
  }
  
  for (y = 1; y < Ny-1; ++y)
  {
    sum[1] += coeffs[L(n,x,y,0)]*Cheb_Tn((int)y,Y);
  }
  sum[1] *= 2;
  
  for (z = 1; z < Nz-1; ++z)
  {
    sum[2] += coeffs[L(n,x,Ny-1,z)]*Cheb_Tn((int)z,Z);
  }
  sum[2] *= 2*Cheb_Tn((int)Ny-1,Y);
  
  for (y = 1; y < Ny-1; ++y)
  {
    sum[3] += coeffs[L(n,x,y,Nz-1)]*Cheb_Tn((int)y,Y);
  }
  sum[3] *= 2*Cheb_Tn((int)Nz-1,Z);
  
  for (y = 0; y < Ny-1; ++y)
  {
    for (z = 0; z < Nz-1; ++z)
      sum[4] += coeffs[L(n,x,y,Nz-1)]*Cheb_Tn((int)y,Y)*Cheb_Tn((int)z,Z);
  }
  sum[4] *= 4;
  
  for (i = 0; i < MAX; ++i)
    interp_v += sum[i];
    
  return interp_v;
}

/* interpolation in X&Y&Z directions.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_XYZ(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double X = interp_s->point->x[0];/* Y coord of the interesting point */
  const double Y = interp_s->point->x[1];/* Y coord of the interesting point */
  const double Z = interp_s->point->x[2];/* Z coord of the interesting point */
  const unsigned *const n = interp_s->point->patch->n;
  const unsigned Nx = n[0];/* number of points in the specified direction */
  const unsigned Ny = n[1];/* number of points in the specified direction */
  const unsigned Nz = n[2];/* number of points in the specified direction */
  const double *coeffs = 0;/* coeffs of expansion of field in Tn */
  enum MAX_E {MAX = 21};
  double sum[MAX];
  double interp_v = 0;
  unsigned x,y,z,i;
  
  make_coeffs_3d(field);
  coeffs = field->v2;
  
  for (i = 0; i < MAX; ++i)
    sum[i] = 0;

  sum[0] += coeffs[0]+coeffs[L(n,0,0,Nz-1)]*Cheb_Tn((int)Nz-1,Z);
  
  for (z = 1; z < Nz-1; ++z)
    sum[1] += coeffs[L(n,0,0,z)]*Cheb_Tn((int)z,Z);
  sum[1] *= 2;
  
  sum[2] = coeffs[L(n,0,Ny-1,0)]+coeffs[L(n,0,Ny-1,Nz-1)]*Cheb_Tn((int)Nz-1,z);
  sum[2] *= Cheb_Tn((int)Ny-1,Y);
  
  for (z = 1; z < Nz-1; ++z)
    sum[3] += coeffs[L(n,0,Ny-1,z)]*Cheb_Tn((int)z,Z);
  sum[3] *= 2;

  for (y = 1; y < Ny-1; ++y)
    sum[4] += coeffs[L(n,0,y,0)]*Cheb_Tn((int)y,Y);
  sum[4] *= 2;
  
  for (y = 1; y < Ny-1; ++y)
    sum[5] += coeffs[L(n,0,y,Nz-1)]*Cheb_Tn((int)y,Y);
  sum[5] *= 2*Cheb_Tn((int)Nz-1,Z);
  
  for (y = 1; y < Ny-1; ++y)
  {
    for (z = 1; z < Nz-1; ++z)
      sum[6] += coeffs[L(n,0,y,z)]*Cheb_Tn((int)y,Y)*Cheb_Tn((int)z,Z);
  }
  sum[6] *= 2;
  
  sum[7] = coeffs[L(n,Nx-1,Ny-1,0)]+coeffs[L(n,Nx-1,Ny-1,Nz-1)];
  sum[7] *= Cheb_Tn((int)Nx-1,X)*Cheb_Tn((int)Ny-1,Y);
  
  for (z = 1; z < Nz-1; ++z)
    sum[8] += coeffs[L(n,Nx-1,Ny-1,z)]*Cheb_Tn((int)z,Z);
  sum[8] *= 2*Cheb_Tn((int)Ny-1,Y);

  for (y = 1; y < Ny-1; ++y)
    sum[9] += coeffs[L(n,Nx-1,y,0)]*Cheb_Tn((int)y,Y);
  sum[9] *= 2*Cheb_Tn((int)Nx-1,X);
  
  for (y = 1; y < Ny-1; ++y)
    sum[10] += coeffs[L(n,Nx-1,y,Nz-1)]*Cheb_Tn((int)y,Y);
  sum[10] *= 2*Cheb_Tn((int)Nx-1,X)*Cheb_Tn((int)Nz-1,Z);
  
  for (y = 1; y < Ny-1; ++y)
  {
    for (z = 1; z < Nz-1; ++z)
      sum[11] += coeffs[L(n,Nx-1,y,z)]*Cheb_Tn((int)y,Y)*Cheb_Tn((int)z,Z);
  }
  sum[11] *= 4*Cheb_Tn((int)Nx-1,X);
  
  for (x = 0; x < Nx-1; ++x)
    sum[12] += coeffs[L(n,x,0,0)]*Cheb_Tn((int)x,X);
  sum[12] *= 2;

  for (x = 0; x < Nx-1; ++x)
    sum[13] += coeffs[L(n,x,0,Nz-1)]*Cheb_Tn((int)x,X);
  sum[13] *= 2*Cheb_Tn((int)Nz-1,z);

  for (x = 0; x < Nx-1; ++x)
    sum[14] += coeffs[L(n,x,Ny-1,0)]*Cheb_Tn((int)x,X);
  sum[14] *= 2*Cheb_Tn((int)Ny-1,Y);
  
  for (x = 0; x < Nx-1; ++x)
    sum[15] += coeffs[L(n,x,Ny-1,Nz-1)]*Cheb_Tn((int)x,X);
  sum[15] *= 2*Cheb_Tn((int)Ny-1,Y)*Cheb_Tn((int)Nz-1,Z);
  
  for (x = 0; x < Nx-1; ++x)
  {
    for (z = 0; z < Ny-1; ++z)
      sum[16] += coeffs[L(n,x,0,z)]*Cheb_Tn((int)z,Z);
    sum[16] *= Cheb_Tn((int)x,X);
  }
  sum[16] *= 2;

  for (x = 0; x < Nx-1; ++x)
  {
    for (z = 0; z < Nz-1; ++z)
      sum[17] += coeffs[L(n,x,Ny-1,z)]*Cheb_Tn((int)z,Z);
    sum[17] *= Cheb_Tn((int)x,X);
  }
  sum[17] *= 2*Cheb_Tn((int)y,Y);
  
  for (x = 0; x < Nx-1; ++x)
  {
    for (y = 0; y < Ny-1; ++y)
      sum[18] += coeffs[L(n,x,y,0)]*Cheb_Tn((int)y,Y);
    sum[18] *= Cheb_Tn((int)x,X);
  }
  sum[18] *= 4;
  
  for (x = 0; x < Nx-1; ++x)
  {
    for (y = 0; y < Ny-1; ++y)
      sum[19] += coeffs[L(n,x,y,Nz-1)]*Cheb_Tn((int)y,Y);
    sum[19] *= Cheb_Tn((int)x,X);
  }
  sum[19] *= 4*Cheb_Tn((int)Nz-1,Z);
  
  for (x = 0; x < Nx-1; ++x)
  {
    for (y = 0; y < Ny-1; ++y)
    {
      for (z = 0; z < Nz-1; ++z)
        sum[20] += coeffs[L(n,x,y,z)]*Cheb_Tn((int)z,Z);
      sum[20] *= Cheb_Tn((int)y,Y);
    }
    sum[20] *= Cheb_Tn((int)x,X);
  }
  sum[20] *= 8;
  
  for (i = 0; i < MAX; ++i)
    interp_v += sum[i];
  
  return interp_v;
}

