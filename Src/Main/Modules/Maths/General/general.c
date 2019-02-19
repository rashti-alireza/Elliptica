/*
// Alireza Rashti
// June 2018
*/

#include "general.h"

/* taking square root of vector v2-v1 which has l double type components.
// ->return value: root mean square of v2-v1.
*/
double rms(const unsigned n, const double *const v2,const double *const v1)
{
  unsigned i;
  double sum;
  sum = 0;
  
  if (v2 == 0 && v1 == 0) return 0;
  else if (v1 == 0)    	  
    for(i = 0; i < n; i++)
      sum += SQR(v2[i]);
  else if (v2 == 0)       
    for(i = 0; i < n; i++)
      sum += SQR(v1[i]);
  else
    for(i = 0; i < n; i++)
      sum += SQR(v2[i]-v1[i]);
    
  sum = sqrt(sum);
  
  return sum;
}

/* taking root means square of vector v2-v1 which has l double type components
// and l is of order of long unsigned.
// ->return value: root mean square of v2-v1
*/
long double rmsL(const long unsigned n, const double *const v2, const double *const v1)
{
  unsigned long i;
  long double sum;
  sum = 0;
  
  if (v2 == 0 && v1 == 0) return 0;
  else if (v1 == 0)    	  
    for(i = 0; i < n; i++)
      sum += SQR(v2[i]);
  else if (v2 == 0)       
    for(i = 0; i < n; i++)
      sum += SQR(v1[i]);
  else
    for(i = 0; i < n; i++)
      sum += SQR(v2[i]-v1[i]);
    
  sum = sqrtl(sum);
  
  return sum;
}

/* taking dot product of two v1 and v2 vector with n components
// ->return value : v2.v1
*/
double dot(const unsigned n, const double *const v2,const double *const v1)
{
  unsigned i;
  double d = 0;
  
  for (i = 0; i < n; i++)
    d += v2[i]*v1[i];
  
  return d;
}

/* taking absolute value of v
// ->return value: absolute(v)
*/
double ABS(const double v)
{
  return v > 0 ? v : -v;
}

/* Chebyshev polynomial of second kind Un(x). x MUST be normalized value.
// ->return value: Un(x)
*/
double Cheb_Un(const int n, const double x)
{
  double u = DBL_MAX;
  
  if (n == 0) 
    u = 1;
  else if (EQL(x,1))
    u = n+1;
  else if (EQL(x,-1)) 
  {
    if (n%2)
      u = -n-1;
    else
      u = n+1;
  }  
  else
  {
    double th = acos(x);
    u = sin((n+1)*th)/sin(th);
  }
  
  return u;
}

/* Chebyshev polynomial of first kind Tn(x). x MUST be normalized value.
// ->return value: Tn(x)
*/
double Cheb_Tn(const int n, const double x)
{
  double t = DBL_MAX;
  
  if (n == 0)
    t = 1;
  else if (EQL(x,1))
    t = 1;
  else if (EQL(x,-1))
  {
    if (n%2)
      t = -1;
    else
      t = 1;
  }
  else
  {
    double th = acos(x);
    t = cos(n*th);
  }
  
  return t;
}

/* second derivative of Cheb_Tn. 
// ->return value: second derivative of Tn
*/
double d2T_dx2(const int n, const double x)
{
  double d = DBL_MAX;
  
  if (n == 0 || n == 1)
    d = 0;
  else if (n == 2)
    d = 4;
  else if (EQL(x,1))
    d = n*n*(n*n-1)/3.0;
  else if (EQL(x,-1))
  {
    if (n%2)
      d = -n*n*(n*n-1)/3.0;
    else
      d = n*n*(n*n-1)/3.0;
  }
  else
  {
    d = n*((n+1)*Cheb_Tn(n,x)-Cheb_Un(n,x))/(x*x-1);
  }
  
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
double sum_0_N_dCi_dfj_by_dTi_dq(const unsigned N,const unsigned j,const double q)
{
  /* some checks */
  if (LSS(q,-1) || GRT(q,1))
    abortEr("q must be in interval [-1,1].");
  if (j >= N)
    abortEr("j must be smaller that N");
    
  double sum = 0;
  const double scale = 0.5/(N-1);/* coming when one does Fourier transformation */
  const double a = acos(q);
  
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
double sum_0_N_dCi_dfj_by_Ti_q(const unsigned N,const unsigned j,const double q)
{
  /* some checks */
  if (LSS(q,-1) || GRT(q,1))
    abortEr("q must be in interval [-1,1].");
  if (j >= N)
    abortEr("j must be smaller that N");
    
  double sum = 0;
  const double scale = 0.5/(N-1);/* coming when one does Fourier transformation */
  const double theta = acos(q);
  
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
              (SQR(1 + N)*Sin((0.5 - N)*b) 
              + 
              (-1 + 2*N + 2*SQR(N))*Sin((0.5 + N)*b) 
              - 
              SQR(N)*Sin((1.5 + N)*b))
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
        sum = -(2*SQR(k)+k);
      }
      else
      {
        k = (N-1)/2.;
        sum = 2*SQR(k)+3*k+1;
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
              (SQR(1 + N)*Sin((0.5 - N)*g1) 
              + 
              (-1 + 2*N + 2*SQR(N))*Sin((0.5 + N)*g1) 
              - 
              SQR(N)*Sin((1.5 + N)*g1))
            )/8.0 
            +
            (
              Power(Csc(g2/2.0),3)
              *
              (SQR(1 + N)*Sin((0.5 - N)*g2) 
              + 
              (-1 + 2*N + 2*SQR(N))*Sin((0.5 + N)*g2) 
              - 
              SQR(N)*Sin((1.5 + N)*g2))
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
double sum_1_N_cos_ia(const unsigned N, const double a)
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
  return ABS(a) > ABS(b) ? ABS(a) : ABS(b);
}
