/*
// Alireza Rashti
// June 2018
*/

#include "general.h"

/* taking square root of vector v2-v1 which has n double type components.
// ->return value: root square of v2-v1; i.e. sqrt(sum{(v2[i]-v1[i])^2}) */
double root_square(const Uint n, const double *const v2,const double *const v1)
{
  Uint i;
  double sum;
  sum = 0;
  
  if (v2 == 0 && v1 == 0) return 0;
  else if (v1 == 0)    	  
    for(i = 0; i < n; i++)
      sum += Pow2(v2[i]);
  else if (v2 == 0)       
    for(i = 0; i < n; i++)
      sum += Pow2(v1[i]);
  else
    for(i = 0; i < n; i++)
      sum += Pow2(v2[i]-v1[i]);
    
  sum = sqrt(sum);
  
  return sum;
}

/* calculate the L2 norm of v2 - v1, i.e. :
// L2Norm = sqrt(sum{(v2[i]-v1[i])^2}/n)
// ->return value: L2 norm. */
double L2_norm(const Uint n, const double *const v2,const double *const v1)
{
  return root_square(n,v2,v1)/sqrt(n);
}

/* calculate the L1 norm of v2 - v1, i.e. :
// L1Norm = sum{|v2[i]-v1[i]|}/n
// ->return value: L1 norm. */
double L1_norm(const Uint n, const double *const v2,const double *const v1)
{
  Uint i;
  double sum;
  sum = 0;
  
  if (v2 == 0 && v1 == 0) return 0;
  else if (v1 == 0)    	  
    for(i = 0; i < n; i++)
      sum += ABSd(v2[i]);
  else if (v2 == 0)       
    for(i = 0; i < n; i++)
      sum += ABSd(v1[i]);
  else
    for(i = 0; i < n; i++)
      sum += ABSd(v2[i]-v1[i]);
    
  return sum/n;
}

/* taking root square of vector v2-v1 which has n double type components
// and n is long Uint type.
// ->return value: root square of v2-v1 */
long double root_square_long(const long Uint n, const double *const v2, const double *const v1)
{
  Uint long i;
  long double sum;
  sum = 0;
  
  if (v2 == 0 && v1 == 0) return 0;
  else if (v1 == 0)    	  
    for(i = 0; i < n; i++)
      sum += Pow2(v2[i]);
  else if (v2 == 0)       
    for(i = 0; i < n; i++)
      sum += Pow2(v1[i]);
  else
    for(i = 0; i < n; i++)
      sum += Pow2(v2[i]-v1[i]);
    
  sum = sqrtl(sum);
  
  return sum;
}

/* taking dot product of two v1 and v2 vector with n components
// ->return value : v2.v1
*/
double dot(const Uint n, const double *const v2,const double *const v1)
{
  Uint i;
  double d = 0;
  
  for (i = 0; i < n; i++)
    d += v2[i]*v1[i];
  
  return d;
}

/* finding summation of functional derivative of Chebyshev coefficients
// multiply by derivative of Chebyshev bases. used for functional derivatives 
// in Jacobian J at Jx=-F for Newton method.
// if f(x) = sum_0_^{n}{C(i)T_i(q)} it calculates:
// 	     sum_0_^{n}{dC(i)/df(j)*dT_i(q)/dq}
// arguments:
// o. N # number of grid points in the direction
// o. j # functional derivative w.r.t field(j)
//	# in other words the position of field that is being varied 
// o. q # normalized point q i.e. -1 <= q <= 1
//
// ->return value: sum_0_^{n}{dC(i)/df(j)*dT_i(q)/dq}. */
double sum_0_N_dCi_dfj_by_dTi_dq(const Uint N,const Uint j,const double q)
{
  /* some checks */
  if (q < -1.0 || q > 1)
    Error0("q must be in interval [-1,1].");
  if (j >= N)
    Error0("j must be smaller that N");
    
  double sum = 0;
  const double scale = 0.5/(N-1);/* coming when one does Fourier transformation */
  const double a = acos(q);
  const double SIGN[2] = {1,-1};
  
  if (j == 0)
  {
    sum = (N-1)*Cheb_Un((int)N-2,q)
          + 
          2*d_dq_sum_1_N_cos_ixb_cos_ixa((int)N-2,0,a);
  }
  else if (j == N-1)
  {
    sum = (N-1)*Cheb_Un((int)N-2,q)*SIGN[j%2]
          + 
          2*d_dq_sum_1_N_cos_ixb_cos_ixa((int)N-2,M_PI,a);
  }
  else
  {
    const double b = j*M_PI/(N-1);
    sum = 2*(N-1)*Cheb_Un((int)N-2,q)*SIGN[j%2]
          + 
          4*d_dq_sum_1_N_cos_ixb_cos_ixa((int)N-2,b,a);
  }
  
  return sum*scale;
}

/* finding summation of functional derivative of Chebyshev coefficients
// multiply by Chebyshev bases. used for functional derivatives 
// in Jacobian J at Jx=-F for Newton method.
// if f(x) = sum_0_^{n}{C(i)T_i(q)} it calculates:
// 	     sum_0_^{n}{dC(i)/df(j)*T_i(q)},T_i(q) = cos(i*acos(q))
//
// arguments:
// o. N # number of grid points in the direction
// o. j # functional derivative w.r.t. field(j)
//	# in other words the position of field that is being varied 
// o. q # normalized point q i.e. -1 <= q <= 1
//
// ->return value: sum_0_^{n}{dC(i)/df(j)*cos(i*acos(q))}. */
double sum_0_N_dCi_dfj_by_Ti_q(const Uint N,const Uint j,const double q)
{
  /* some checks */
  if (q < -1.0 || q > 1.0)
    Error0("q must be in interval [-1,1].");
  if (j >= N)
    Error0("j must be smaller that N");
    
  double sum = 0;
  const double scale = 0.5/(N-1);/* coming when one does Fourier transformation */
  const double theta = acos(q);
  const double SIGN[2] = {1,-1};
  
  if (j == 0)
  {
    sum = 1 + cos((N-1)*theta) + 2*sum_1_N_cos_ia(N-2,theta);
  }
  else if (j == N-1)
  {
    sum = 1 + SIGN[j%2]*cos((N-1)*theta)
            + sum_1_N_cos_ia(N-2,theta+M_PI)
            + sum_1_N_cos_ia(N-2,theta-M_PI);
  }
  else
  {
    const double alpha = j*M_PI/(N-1);
    sum = 1 + SIGN[j%2]*cos((N-1)*theta)
            + sum_1_N_cos_ia(N-2,theta+alpha)
            + sum_1_N_cos_ia(N-2,theta-alpha);
    sum *= 2;
  }
  
  return sum*scale;
}

/* calculate:  d/dq{sum_1_^{N}(cos(i.b).Ti(q))} = 
// 	       sum_1_^{N}(cos(i.b).dTi(q)/dq)   = 
//             d/dq{sum_1_^{N}(cos(i.b).cos(i.a))}
// note: a and b are angels.
// ->return value: d/dq{sum_1_^{N}(cos(i.b).cos(i.a))} , q = cos(a) */
double d_dq_sum_1_N_cos_ixb_cos_ixa(const int N, const double b,const double a)
{
  double sum = 0;
  
  if (EQL(a,0))
  {
    if (EQL(b,0))
    {
      sum = N*(N+1)*(2*N+1)/6.0;
    }
    else
    {
      sum = (
              Power(Csc(b/2.0),3)
              *
              (Pow2(1 + N)*Sin((0.5 - N)*b) 
              + 
              (-1 + 2*N + 2*Pow2(N))*Sin((0.5 + N)*b) 
              - 
              Pow2(N)*Sin((1.5 + N)*b))
            )/8.0;
    }
  }
  else if (EQL(a,M_PI))
  {
    if (EQL(b,0))
    {
      double k;
      
      if (N%2 == 0)
      {
        k = N/2.;
        sum = -(2*Pow2(k)+k);
      }
      else
      {
        k = (N-1)/2.;
        sum = 2*Pow2(k)+3*k+1;
      }
    }
    else if (EQL(b,M_PI))
    {
      sum = -N*(N+1)*(2*N+1)/6.0;
    }
    else
    {
      const double g1 = a+b;
      const double g2 = a-b;
      sum = (
              Power(Csc(g1/2.0),3)
              *
              (Pow2(1 + N)*Sin((0.5 - N)*g1) 
              + 
              (-1 + 2*N + 2*Pow2(N))*Sin((0.5 + N)*g1) 
              - 
              Pow2(N)*Sin((1.5 + N)*g1))
            )/8.0 
            +
            (
              Power(Csc(g2/2.0),3)
              *
              (Pow2(1 + N)*Sin((0.5 - N)*g2) 
              + 
              (-1 + 2*N + 2*Pow2(N))*Sin((0.5 + N)*g2) 
              - 
              Pow2(N)*Sin((1.5 + N)*g2))
            )/8.0;
      sum /= -2;
    }
  }
  else
  {
    if (EQL(a-b,0))
    {
      const double g1 = a+b;
      sum = (
             Csc(g1/2.)*((1 + 2*N)*Cos((0.5 + N)*g1) 
             - 
             Cot(g1/2.)*Sin((0.5 + N)*g1))
            )/4.;
      sum /= -2*Sin(a);
    }
    else
    {
      const double g1 = a+b;
      const double g2 = a-b;
      sum = (
             Csc(g1/2.)*((1 + 2*N)*Cos((0.5 + N)*g1) 
             - 
             Cot(g1/2.)*Sin((0.5 + N)*g1))
            )/4.
            +
            (
             Csc(g2/2.)*((1 + 2*N)*Cos((0.5 + N)*g2) 
             - 
             Cot(g2/2.)*Sin((0.5 + N)*g2))
            )/4.;
      sum /= -2*Sin(a);      
    }
  }
  
  return sum;
}

/* ->return value: sum_1_^{N}{cos(i*theta)} */
double sum_1_N_cos_ia(const Uint N, const double a)
{
  double sum = 0;
  
  if (EQL(a,0) || EQL(a,2*M_PI))
  {
    sum = N;
  }
  else
  {
    sum = -0.5 + 0.5*sin((N+0.5)*a)/sin(0.5*a);
  }
  
  return sum;
}

/* given two number a and b it returns the number with
// MAXIMUM MAGNITUDE.
// ->return value: maximum magnitude between a and b
*/
double MaxMag_d(const double a,const double b)
{
  return ABSd(a) > ABSd(b) ? ABSd(a) : ABSd(b);
}

/* given a double array and its dimension N, 
// it finds the maximum magnitude of the elements, i.e. L infinity norm
// ->return value: maximum magnitude of a double array (L infinity). */
double L_inf(const Uint N,const double *const v)
{
  if (!v || !N)
    Error0("The given array is empty.");
  
  double max = ABSd(v[0]);
  Uint i;
  
  for (i = 1; i < N; ++i)
  {
    double a = ABSd(v[i]);
    
    //if (GRT(a,max))
    if (a > max)
      max = a;
  }
  
  return max;
}

/* arctang in full range based on sign of x and y since: y/x = tan(phi)
// ->return value: phi in [0,2Pi) */
double arctan(const double y,const double x)
{
  double atg = DBL_MAX;

  if (EQL(x,0))
  {
    if (EQL(y,0))
      atg = 0;
    else if (y > 0)
      atg = M_PI/2;
    else
      atg = 3*M_PI/2;
  }
  else if (x > 0)/* x + */
  {
    if (EQL(y,0))
      atg = 0;
    else if (y > 0)
      atg = atan(y/x);
    else
      atg = 2*M_PI-atan(-y/x);
  }
  else/* x - */
  {
    if (EQL(y,0))
      atg = M_PI;
    else if (y > 0)
      atg = M_PI-atan(y/-x);
    else
      atg = M_PI+atan(-y/-x);
  }
  
  return atg;
}

/* given phi in range [0,2Pi) it finds out the sign of x and y in y/x = tan(phi) */
void arctan_argument_signum(double *const y_sign,double *const x_sign,const double phi)
{
  if (phi < 0 || phi >= 2*M_PI)
    Error0("Bad argument, phi must be in [0,2Pi).\n");
  
  if (GRTEQL(phi,0) && LSSEQL(phi,M_PI/2))
  {
    *y_sign = 1;
    *x_sign = 1;
  }
  else if (GRTEQL(phi,M_PI/2) && LSSEQL(phi,M_PI))
  {
    *y_sign = 1;
    *x_sign = -1;
  }
  else if (GRTEQL(phi,M_PI) && LSSEQL(phi,3*M_PI/2))
  {
    *y_sign = -1;
    *x_sign = -1;
  }
  else if (GRTEQL(phi,3*M_PI/2) && LSS(phi,2*M_PI))
  {
    *y_sign = -1;
    *x_sign = 1;
  }
  else
    Error0(NO_OPTION);
    
}
