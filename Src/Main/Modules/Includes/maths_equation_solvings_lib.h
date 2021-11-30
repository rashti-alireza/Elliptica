#ifndef maths_equation_solvings_LIB_H
#define maths_equation_solvings_LIB_H
#include "elliptica_system_lib.h"

#include "maths_linear_algebra_lib.h"

#define MAX_STR_MATH_EQ_SOLVE_LIB (400)

/* NOTE: SPECTRAL_JACOBIAN_MATRIX_FORM and 
// SPECTRAL_JACOBIAN_ANALYTIC_FORM must be mutually exclusive. */
#define SPECTRAL_JACOBIAN_MATRIX_FORM (1)
#define SPECTRAL_JACOBIAN_ANALYTIC_FORM (0)

/* Define macros for analytic calculation of df/du Jacobian */

/* basic trig math */
#define J__Cos(a) cos((a))
#define J__Sin(a) sin((a))
#define J__Csc(a) (1./sin((a)))
#define J__Cot(a) (1./tan((a)))

/* 2*M_PI */
#define J__2M_PI (6.283185307179586)

/* Kronecker Delta */
#define J__KD(i,j)  ( (i)==(j) ? 1. : 0.)

/* eta_i in Elliptica's paper, NOTE: n = patch->n[?] - 1 */
#define J__eta(i,n) ( (i) == 0 || (i) == (n) ? 1. : 2. )


/* (-1)^i */
#define J__sign(i) (((i)%2) ? -1. : 1.)

/* quick theta: 
// NOTE1: assuming Chebyshev Extrema points.
// NOTE2: assuming patch is defined. */
#define J__theta(i,X_axis) ( (i)*M_PI/((patch->n[(X_axis)])-1.) )


/* normalization, NOTE: n = patch->n */
#define J__norm(n) ( 0.5/((n)-1.) )


/* dX/dx */
#define J__dX_dx(patch,ijk,dX_axis,dx_axis) \
  ( (patch)->JacobianT->dX_dx[(dX_axis)][(dx_axis)][(ijk)] )


/* d2X/dxdy */
#define J__d2X_dxdy(patch,ijk,dX_axis,dxdy_axis) \
  ( (patch)->JacobianT->d2X_dxdy[(dX_axis)][(dxdy_axis)][(ijk)] )

/* dN/dX */
#define J__dN_dX(patch,dX_axis) \
  ( (patch)->JacobianT->dN_dX[(dX_axis)] )


/* sum_{n=0}^{N} cos(n lambda) = 
// 0.5 + 0.5*( sin( (N+0.5)*(lambda) ) ) / ( sin( 0.5*(lambda) ) ),
// N0 = N+0.5. */
#define J__sum_0_N_cos_nlambda(N,N0,lambda) \
  ( \
    0.5 + 0.5*( sin( (N0)*(lambda) ) ) / ( sin( 0.5*(lambda) ) ) \
  )

/* d/dlambda sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define J__d_dlambda_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J__2M_PI) ?\
    (0.0):\
    ( \
      J__Csc(0.5*(lambda))*(2.*(N0)*J__Cos((lambda)*(N0)) - \
      J__Cot(0.5*(lambda))*J__Sin((lambda)*(N0)))\
    )*0.25\
  )


double J__v_lambda       = ????
double J__v_half_lambda  = 0.5*(J__v_lambda);
double J__v_N0_lambda[3] = {J__v_N0[0]*J__v_lambda,
                            J__v_N0[1]*J__v_lambda,
                            J__v_N0[2]*J__v_lambda};

double J__v_cos_lambda = cos(J__v_lambda);
double J__v_sin_lambda = sin(J__v_lambda);
double J__v_csc_lambda = 1./J__v_sin_lambda;

double J__v_cos_half_lambda = cos(J__v_half_lambda);
double J__v_sin_half_lambda = sin(J__v_half_lambda);
double J__v_csc_half_lambda = 1./J__v_sin_half_lambda;

double J__v_cos_N0_lambda[3] = {cos(J__v_N0_lambda[0]),
                                cos(J__v_N0_lambda[1]),
                                cos(J__v_N0_lambda[2])};
double J__v_sin_N0_lambda[3] = {sin(J__v_N0_lambda[0]),
                                sin(J__v_N0_lambda[1]),
                                sin(J__v_N0_lambda[2])};
double J__v_csc_N0_lambda[3] = {1./J__v_sin_N0_lambda[0],
                                1./J__v_sin_N0_lambda[1],
                                1./J__v_sin_N0_lambda[2]};

/* d^2/dlambda^2 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define J__d2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J__2M_PI) ?\
    -(Pow3(N)/3.+Pow2(N)/2.+N/6.)/* simplified, don't forget - sign! */ :\
    (\
      J__Csc(0.5*(lambda))*(-4.*(N0)*J__Cos((lambda)*(N0))*J__Cot(0.5*(lambda)) + \
      (-1. - 4.*Pow2(N0) + 2.*Pow2(J__Csc(0.5*(lambda))))*J__Sin((lambda)*(N0)))\
    )*0.125\
  )


double J__v_2pow2_half_csc_lambda = 2.*Pow2(J__v_csc_half_lambda);

/* d^3/dlambda^3 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define J__d3_dlambda3_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J__2M_PI) ?\
    (0.0):\
    (\
      pow(J__Csc(0.5*(lambda),3))*(2*(N0)*\
        (9 - 4*Pow2((N0)) + (3 + 4*Pow2((N0)))*J__Cos((lambda)))*\
        J__Cos((lambda)*(N0)) - (11 - 12*Pow2((N0)) + J__Cos((lambda)) + \
          12*Pow2((N0))*J__Cos((lambda)))*J__Cot(0.5*(lambda))*J__Sin((lambda)*(N0)))\
    )/32.\
  )

/* d^4/dlambda^4 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define J__d4_dlambda4_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J__2M_PI) ?\
    (Pow4(N)*(N/5.+0.5)+Pow3(N)/3.-N/30.)/* simplified */:\
    (\
      pow(J__Csc(0.5*(lambda)),5)*(-16*(N0)*\
        (11 - 4*Pow2((N0)) + J__Cos((lambda)) + \
          4*Pow2((N0))*J__Cos((lambda)))*J__Cos((lambda)*(N0))*J__Sin((lambda)) + \
       (115 - 120*Pow2((N0)) + 48*Pow4((N0)) + \
          (76 + 96*Pow2((N0)) - 64*Pow4((N0)))*J__Cos((lambda)) + \
          (1 + 24*Pow2((N0)) + 16*Pow4((N0)))*J__Cos(2*(lambda)))*\
        J__Sin((lambda)*(N0)))\
    )/256.\
  )


/* -> eta_j{d/dX(df/du)=d/dX (2*sum_0^N (Tn(Xj) Tn(X)) -1 -(-1)^j *T_{N}(X))},
// NOTE: X = cos(th), N = patch->n[?]-1. */
#define J__d2f_dudX(thi,thj,N,j) \
   (J__eta(j,N)*( d_dXi_2xsum_0_N_Tnj_Tni(thi,thj,N) - J__sign(j)*dT_dx((int)(N),cos(thi)) ))

/* -> d2/dX^2(df/du)=d2/dX^2 (2*sum_0^N (Tn(Xj) Tn(X)) -1 -(-1)^j *T_{N}(X)),
// NOTE: X = cos(th)), N = patch->n[?]-1. */
#define J__d3f_dudXdX(thi,thj,N,j) \
   (J__eta(j,N)*( d2_dXi2_2xsum_0_N_Tnj_Tni(thi,thj,N) - J__sign(j)*d2T_dx2((int)(N),cos(thi)) ))

/* normalization * coords jacobian * J__d2f_dudX */
#define J__d2f_dudx(patch, dx_axis, X_axis, ijk, lmn, qi,qj) \
  ( J__norm(patch->n[X_axis])*J__dX_dx(patch,ijk,X_axis,dx_axis)*J__dN_dX(patch,X_axis)*\
    J__d2f_dudX(J__theta(qi,X_axis),J__theta(qj,X_axis),patch->n[X_axis]-1,qj) )

/* normalization * coords jacobian * J__d3f_dudXdX */
#define J__d3f_dudxdy(patch, dx_axis, dy_axis, dxdy_axis,X_axis, ijk, lmn, qi,qj) \
  ( \
    J__norm(patch->n[X_axis])*J__dN_dX(patch,X_axis)*\
    ( \
      J__d2X_dxdy(patch,ijk,X_axis,dxdy_axis)*\
      J__d2f_dudX(J__theta(qi,X_axis),J__theta(qj,X_axis),patch->n[X_axis]-1,qj) +     \
      J__dX_dx(patch,ijk,X_axis,dx_axis)*J__dX_dx(patch,ijk,X_axis,dy_axis)*J__dN_dX(patch,X_axis)* \
      J__d3f_dudXdX(J__theta(qi,X_axis),J__theta(qj,X_axis),patch->n[X_axis]-1,qj)     \
    )\
  )

/* X_i = cos(theta_i), X_j = cos(theta_j), N = patch->n-1  */
#define J__d_dXi_2xsum_0_N_Tnj_Tni(thi,thj,N)
( EQL(thi,0.) ? 
  -2.*J__d2_dlambda2_sum_0_N_cos_nlambda(N,N0,thj) :
  EQL(thi,M_PI) ?
  ( (lambda = thj+M_PI, _sum  = J__d2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda),
     lambda = thj-M_PI, _sum += J__d2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda),_sum)):
  ())
  
{
  double sum = 0.;
  double N0 = N+0.5;

  if (EQL(thi,0.))
  {
    sum = -2.*Jd2_dlambda2_sum_0_N_cos_nlambda(N,N0,thj);
  }
  else if (EQL(thi,M_PI))
  {
    double lambda = thj+M_PI;
    sum = Jd2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda);

    lambda = thj-M_PI;
    sum += Jd2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda);
  }
  else
  {
    double lambda = thi+thj;
    double dthi_dX   = -1./sin(thi);
    
    sum = Jd_dlambda_sum_0_N_cos_nlambda(N,N0,lambda);
    
    lambda = thi-thj;
    sum += Jd_dlambda_sum_0_N_cos_nlambda(N,N0,lambda);
    
    sum *= dthi_dX;
  }
  
  return sum;
}

/* -> d^2/dX^2 2*sum_0^N (Tn(Xj) Tn(X))| X = Xi.
// X = cos(th). */
static double
d2_dXi2_2xsum_0_N_Tnj_Tni(double thi/* X_i = cos(theta_i) */,
                          double thj/* X_i = cos(theta_i) */,
                          Uint N/* the sum upper limit */)
{
  double sum = 0.;
  double N0 = N+0.5;
  
  if (EQL(thi,0.))
  {
    sum = 
      Jd4_dlambda4_sum_0_N_cos_nlambda(N,N0,thj) +
      Jd2_dlambda2_sum_0_N_cos_nlambda(N,N0,thj);
    sum *= 2./3.;
  }
  else if (EQL(thi,M_PI))
  {
    double lambda = thj+M_PI;
    sum = 
      Jd4_dlambda4_sum_0_N_cos_nlambda(N,N0,lambda) +
      Jd2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda);
    
    lambda = thj-M_PI;
    sum += 
      Jd4_dlambda4_sum_0_N_cos_nlambda(N,N0,lambda) +
      Jd2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda);
      
    sum /= 3.;
  }
  else
  {
    double sin_thi   = sin(thi);
    double d2thi_dX2 = -cos(thi)/(Pow3(sin_thi));
    double dthi_dX   = -1./sin_thi;
    double lambda    = thi+thj;
    
    sum = d2thi_dX2*Jd_dlambda_sum_0_N_cos_nlambda(N,N0,lambda) +
          Pow2(dthi_dX)*Jd2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda);
    
    lambda = thi-thj;
    sum += d2thi_dX2*Jd_dlambda_sum_0_N_cos_nlambda(N,N0,lambda) +
           Pow2(dthi_dX)*Jd2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda);
  }
  
  return sum;
}




#if SPECTRAL_JACOBIAN_MATRIX_FORM

#define Header_Jacobian /* nothing needed yet! */

#define Footer_Jacobian /* free and clean stuffs */

/* it compactifies the prepration of Jacobian of derivatives */
#define Init_Jacobian(xNAME) \
  const char *types_##xNAME[] = {#xNAME,0};\
  prepare_Js_jacobian_eq(patch,types_##xNAME);\
  Matrix_T *m_##xNAME = get_j_matrix(patch,#xNAME);\
  fJs_T *f_##xNAME    = get_j_reader(m_##xNAME);

#define Free_Jacobian(xNAME) /* it's not design yet! maybe in future 
                             // I want to remove Jacobian after each 
                             // population. so this is a place holder. */

#define d2f_dxdu_Jacobian(patch,dx_axis,ijk,lmn,xNAME) \
  ( f_##xNAME(m_##xNAME, ijk, lmn) )

#define d3f_dx2du_Jacobian(patch,dx_axis,ijk,lmn,xNAME) \
  ( f_##xNAME(m_##xNAME, ijk, lmn) )


#elif SPECTRAL_JACOBIAN_ANALYTIC_FORM

#define Init_Jacobian(xNAME) /* nothing needed! */

#define Free_Jacobian(xNAME) /* nothing needed! */

#define Header_Jacobian /* set some variables and initialization */
  const Uint J__v_nm1[3] = {patch->n[0]-1,
                            patch->n[1]-1,
                            patch->n[2]-1};

  const double J__v_pi_o_nm1[3] = {M_PI/J__v_nm1[0],
                                   M_PI/J__v_nm1[1], 
                                   M_PI/J__v_nm1[2]};
  const double J__v_norm[3] = {0.5/J__v_nm1[0],
                               0.5/J__v_nm1[1],
                               0.5/J__v_nm1[2]};

  const double J__v_N0[3] = {0.5 + J__v_nm1[0],
                             0.5 + J__v_nm1[1],
                             0.5 + J__v_nm1[2]};

  const double J__v_c1_d2[3] = {-(Pow3(J__v_nm1[0])/3.+Pow2(J__v_nm1[0])/2.+J__v_nm1[0]/6.),
                                -(Pow3(J__v_nm1[1])/3.+Pow2(J__v_nm1[1])/2.+J__v_nm1[1]/6.),
                                -(Pow3(J__v_nm1[2])/3.+Pow2(J__v_nm1[2])/2.+J__v_nm1[2]/6.)};

  const double J__v_c2_d2[3] = {-1. - 4.*Pow2(J__v_N0[0]), 
                                -1. - 4.*Pow2(J__v_N0[1]),
                                -1. - 4.*Pow2(J__v_N0[2])};

  const double J__v_c1_d4[3] = {(Pow4(J__v_nm1[0])*(J__v_nm1[0]/5.+0.5)+Pow3(J__v_nm1[0])/3.-J__v_nm1[0]/30.),
                                (Pow4(J__v_nm1[1])*(J__v_nm1[1]/5.+0.5)+Pow3(J__v_nm1[1])/3.-J__v_nm1[1]/30.),
                                (Pow4(J__v_nm1[2])*(J__v_nm1[2]/5.+0.5)+Pow3(J__v_nm1[2])/3.-J__v_nm1[2]/30.)};

  const double J__v_c2_d4[3] = {4.*Pow2(J__v_nm1[0]),
                                4.*Pow2(J__v_nm1[1]),
                                4.*Pow2(J__v_nm1[2])};

  const double J__v_c3_d4[3] = {115. - 120.*Pow2(J__v_nm1[0]) + 48.*Pow4(J__v_nm1[0]),
                                115. - 120.*Pow2(J__v_nm1[1]) + 48.*Pow4(J__v_nm1[1]),
                                115. - 120.*Pow2(J__v_nm1[2]) + 48.*Pow4(J__v_nm1[2])};
                                
  const double J__v_c4_d4[3] = {(76. + 96.*Pow2(J__v_nm1[0]) - 64.*Pow4(J__v_nm1[0])),
                                (76. + 96.*Pow2(J__v_nm1[1]) - 64.*Pow4(J__v_nm1[1])),
                                (76. + 96.*Pow2(J__v_nm1[2]) - 64.*Pow4(J__v_nm1[2]))};


  const double J__v_c5_d4[3] = {(1. + 24.*Pow2(J__v_nm1[0]) + 16.*Pow4(J__v_nm1[0])),
                                (1. + 24.*Pow2(J__v_nm1[1]) + 16.*Pow4(J__v_nm1[1])),
                                (1. + 24.*Pow2(J__v_nm1[2]) + 16.*Pow4(J__v_nm1[2]))};

  /* set d^n Cheb/dx^n */
  double *dT_dx[3]   = patch->Solving_Man_T->jacobian_workspace->dT_dx;
  double *d2T_dx2[3] = patch->Solving_Man_T->jacobian_workspace->d2T_dx2;
  /* populate dT_dx if empty */
  if (!dT_dx[0] || !dT_dx[1] || !dT_dx[2])
  {
    Free(dT_dx[0]);
    dT_dx[0] = alloc_double(patch->nn);
    
    Free(dT_dx[1]);
    dT_dx[1] = alloc_double(patch->nn);
    
    Free(dT_dx[2]);
    dT_dx[2] = alloc_double(patch->nn);
    
    /* set */
    for (Uint _ijk; _ijk < patch->nn; ++_ijk)
    {
      Uint _ip,_jp,_kp;
      double _x[3];
      
      ijk_to_i_j_k(_ijk,patch->n,&_ip,&_jp,&_kp);
      _x[0] =  cos(_ip*J__v_pi_o_nm1[0]);
      _x[1] =  cos(_jp*J__v_pi_o_nm1[1]);
      _x[2] =  cos(_kp*J__v_pi_o_nm1[2]);
      
      dT_dx[0][_ijk] = dT_dx(int(J__v_nm1[0]),_x[0]);
      dT_dx[1][_ijk] = dT_dx(int(J__v_nm1[1]),_x[1]);
      dT_dx[2][_ijk] = dT_dx(int(J__v_nm1[2]),_x[2]);
      
      d2T_dx2[0][_ijk] = d2T_dx2(int(J__v_nm1[0]),_x[0]);
      d2T_dx2[1][_ijk] = d2T_dx2(int(J__v_nm1[1]),_x[1]);
      d2T_dx2[2][_ijk] = d2T_dx2(int(J__v_nm1[2]),_x[2]);
    }
    
    /* save */
    for (Uint _i; _i < 3; ++_i)
    {
      patch->Solving_Man_T->jacobian_workspace->dT_dx[_i]   = dT_dx[_i];
      patch->Solving_Man_T->jacobian_workspace->d2T_dx2[_i] = d2T_dx2[_i];
    }
  }

#define Footer_Jacobian /* free and clean stuffs */

#define d2f_dxdu_Jacobian(patch,dx_axis,ijk,lmn,xNAME) \
  d2f_dxdu_spectral_Jacobian_analytic(patch,dx_axis,ijk,lmn)

#define d3f_dx2du_Jacobian(patch,dxdy_axis,ijk,lmn,xNAME) \
  d3f_dxdydu_spectral_Jacobian_analytic(patch,dxdy_axis,ijk,lmn)
  

#endif


/* defining some macros to improve the readability and simplicity */

/* macros for jacobian of equations */
#define DDM_SCHUR_JACOBIAN_EQ_DECLARE \
  Patch_T *const patch  = vp1;\
  DDM_Schur_Complement_T *const S = vp2;\
  double **const B = S->B->reg->A;\
  double **E_Trans;\
  const Uint *const node = S->inv;\
  const Uint Ni = S->Oi;/* number of inner mesh nodes */\
  const Uint Nj = S->NS;/* number of inner mesh+outer-boundary + inner-boundary nodes */\
  const Uint K0 = S->NS;/* number of inner mesh+outer-boundary + inner-boundary nodes */\
  const Uint Nk = patch->nn;/* total number of nodes */\
  Uint i,j,k;

/* macro for B part of jacobian */
#define DDM_SCHUR_JACOBIAN_EQ_Bpart_OPEN \
  for (i = 0; i < Ni; ++i)\
  {\
    ijk = node[i];\
    Uint _i,_j,_k; ijk_to_i_j_k(ijk,patch->n,&_i,&_j,&_k);
    
    for (j = 0; j < Nj; ++j)\
    {\
      lmn = node[j];
      Uint _l,_m,_n; ijk_to_i_j_k(lmn,patch->n,&_l,&_m,&_n);
      int _dx_axis,_dy_axis,_dxdy_axis;
      double _sum;/* for J__d_dXi_2xsum_0_N_Tnj_Tni and J__d2_dXi2_2xsum_0_N_Tnj_Tni. */
      /* first order Jacobian, NOTE: don't change d2f_dxdu_spectral_Jacobian_analytic name */
      double d2f_dxdu_spectral_Jacobian_analytic[3];
      
      _dx_axis = 0;
      d2f_dxdu_spectral_Jacobian_analytic[_dx_axis] = 
        J__d2f_dudx(patch,_dx_axis,0,ijk,lmn,_i,_l)*J__KD(_j,_m)*J__KD(_k,_n)+
        J__d2f_dudx(patch,_dx_axis,1,ijk,lmn,_j,_m)*J__KD(_i,_l)*J__KD(_k,_n)+
        J__d2f_dudx(patch,_dx_axis,2,ijk,lmn,_k,_n)*J__KD(_i,_l)*J__KD(_j,_m);

      _dx_axis = 1;
      d2f_dxdu_spectral_Jacobian_analytic[_dx_axis] = 
        J__d2f_dudx(patch,_dx_axis,0,ijk,lmn,_i,_l)*J__KD(_j,_m)*J__KD(_k,_n)+
        J__d2f_dudx(patch,_dx_axis,1,ijk,lmn,_j,_m)*J__KD(_i,_l)*J__KD(_k,_n)+
        J__d2f_dudx(patch,_dx_axis,2,ijk,lmn,_k,_n)*J__KD(_i,_l)*J__KD(_j,_m);
   
      _dx_axis = 2;
      d2f_dxdu_spectral_Jacobian_analytic[_dx_axis] = 
        J__d2f_dudx(patch,_dx_axis,0,ijk,lmn,_i,_l)*J__KD(_j,_m)*J__KD(_k,_n)+
        J__d2f_dudx(patch,_dx_axis,1,ijk,lmn,_j,_m)*J__KD(_i,_l)*J__KD(_k,_n)+
        J__d2f_dudx(patch,_dx_axis,2,ijk,lmn,_k,_n)*J__KD(_i,_l)*J__KD(_j,_m);

      /* second order Jacobian, NOTE: don't change d3f_dxdydu_spectral_Jacobian_analytic name */
      double d3f_dxdydu_spectral_Jacobian_analytic[6];

      _dxdy_axis = 0;
      _dx_axis = 0;
      _dy_axis = 0;
      d3f_dxdydu_spectral_Jacobian_analytic[_dxdy_axis] = 
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,0,ijk,lmn,_i,_l)*J__KD(_j,_m)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,0,ijk,lmn,_i,_l)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,1,ijk,lmn,_j,_m)*J__KD(_i,_l)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,1,ijk,lmn,_j,_m)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l) +
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,2,ijk,lmn,_k,_n)*J__KD(_j,_m)*J__KD(_i,_l) +
          J__d2f_dudx(patch,_dx_axis,2,ijk,lmn,_k,_n)*
            (
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l)
            );
    
      _dxdy_axis = 1;
      _dx_axis = 0;
      _dy_axis = 1;
      d3f_dxdydu_spectral_Jacobian_analytic[_dxdy_axis] = 
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,0,ijk,lmn,_i,_l)*J__KD(_j,_m)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,0,ijk,lmn,_i,_l)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,1,ijk,lmn,_j,_m)*J__KD(_i,_l)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,1,ijk,lmn,_j,_m)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l) +
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,2,ijk,lmn,_k,_n)*J__KD(_j,_m)*J__KD(_i,_l) +
          J__d2f_dudx(patch,_dx_axis,2,ijk,lmn,_k,_n)*
            (
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l)
            );
        
      _dxdy_axis = 2;
      _dx_axis = 0;
      _dy_axis = 2;
      d3f_dxdydu_spectral_Jacobian_analytic[_dxdy_axis] = 
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,0,ijk,lmn,_i,_l)*J__KD(_j,_m)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,0,ijk,lmn,_i,_l)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,1,ijk,lmn,_j,_m)*J__KD(_i,_l)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,1,ijk,lmn,_j,_m)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l) +
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,2,ijk,lmn,_k,_n)*J__KD(_j,_m)*J__KD(_i,_l) +
          J__d2f_dudx(patch,_dx_axis,2,ijk,lmn,_k,_n)*
            (
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l)
            );
        
      _dxdy_axis = 3;
      _dx_axis = 1;
      _dy_axis = 1;
      d3f_dxdydu_spectral_Jacobian_analytic[_dxdy_axis] = 
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,0,ijk,lmn,_i,_l)*J__KD(_j,_m)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,0,ijk,lmn,_i,_l)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,1,ijk,lmn,_j,_m)*J__KD(_i,_l)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,1,ijk,lmn,_j,_m)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l) +
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,2,ijk,lmn,_k,_n)*J__KD(_j,_m)*J__KD(_i,_l) +
          J__d2f_dudx(patch,_dx_axis,2,ijk,lmn,_k,_n)*
            (
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l)
            );
        
      _dxdy_axis = 4;
      _dx_axis = 1;
      _dy_axis = 2;
      d3f_dxdydu_spectral_Jacobian_analytic[_dxdy_axis] = 
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,0,ijk,lmn,_i,_l)*J__KD(_j,_m)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,0,ijk,lmn,_i,_l)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,1,ijk,lmn,_j,_m)*J__KD(_i,_l)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,1,ijk,lmn,_j,_m)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l) +
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,2,ijk,lmn,_k,_n)*J__KD(_j,_m)*J__KD(_i,_l) +
          J__d2f_dudx(patch,_dx_axis,2,ijk,lmn,_k,_n)*
            (
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l)
            );
        
      _dxdy_axis = 5;
      _dx_axis = 2;
      _dy_axis = 2;
      d3f_dxdydu_spectral_Jacobian_analytic[_dxdy_axis] = 
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,0,ijk,lmn,_i,_l)*J__KD(_j,_m)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,0,ijk,lmn,_i,_l)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,1,ijk,lmn,_j,_m)*J__KD(_i,_l)*J__KD(_k,_n) +
          J__d2f_dudx(patch,_dx_axis,1,ijk,lmn,_j,_m)*
            (
              J__KD(_k,_n)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l) +
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,2,ijk,lmn,_k,_n)
            ) +
          J__d3f_dudxdy(patch,_dx_axis,_dy_axis,_dxdy_axis,2,ijk,lmn,_k,_n)*J__KD(_j,_m)*J__KD(_i,_l) +
          J__d2f_dudx(patch,_dx_axis,2,ijk,lmn,_k,_n)*
            (
              J__KD(_i,_l)*J__d2f_dudx(patch,_dy_axis,1,ijk,lmn,_j,_m) +
              J__KD(_j,_m)*J__d2f_dudx(patch,_dy_axis,0,ijk,lmn,_i,_l)
            );


#define DDM_SCHUR_JACOBIAN_EQ_Bpart_CLOSE \
    }/* end of for (i = 0; i < Ni; ++i) */\
  }/* end of for (j = 0; j < Nj; ++j) */

/* macros for E part of jacobian */
#define DDM_SCHUR_JACOBIAN_EQ_Epart_OPEN \
  if (S->NI)/* if there is any interface points then E is needed */\
  {\
    E_Trans = S->E_Trans->reg->A;\
    for (k = K0; k < Nk; ++k)\
    {\
      lmn = node[k];\
      j = k-K0;\
      for (i = 0; i < Ni; ++i)\
      {\
        ijk = node[i];

#define DDM_SCHUR_JACOBIAN_EQ_Epart_CLOSE \
     }/* end of for (i = 0; i < Ni; ++i) */\
    }/* end of for (k = K0; k < Nk; ++k) */\
  }/* end of if (S->NI) */


/* macros for jacobian of boundary condition equations */
#define DDM_SCHUR_JACOBIAN_BC_DECLARE \
  Patch_T *const patch  = vp1;\
  DDM_Schur_Complement_T *const S = vp2;\
  double **const B = S->B->reg->A;\
  double **E_Trans;\
  const Uint *const node = S->inv;\
  const Uint I0 = S->Oi;/* number of inner mesh nodes */\
  const Uint Ni = S->NS;/* number of inner mesh+outer-boundary + inner-boundary nodes */\
  const Uint Nj = S->NS;/* number of inner mesh+outer-boundary + inner-boundary nodes */\
  const Uint K0 = S->NS;/* number of inner mesh+outer-boundary + inner-boundary nodes */\
  const Uint Nk = patch->nn;/* total number of nodes */\
  Uint i,j,k;

/* macro for B part of outer boundary jacobian */
#define DDM_SCHUR_JACOBIAN_BC_Bpart_OPEN \
  for (i = I0; i < Ni; ++i)\
  {\
    ijk = node[i];\
    for (j = 0; j < Nj; ++j)\
    {\
      lmn = node[j];

#define DDM_SCHUR_JACOBIAN_BC_Bpart_CLOSE \
    }/* end of for (i = I0; i < Ni; ++i) */\
  }/* end of for (j = 0; j < Nj; ++j) */

/* macros for E part of jacobian */
#define DDM_SCHUR_JACOBIAN_BC_Epart_OPEN \
  if (S->NI)/* if there is any interface points then E is needed */\
  {\
    E_Trans = S->E_Trans->reg->A;\
    for (k = K0; k < Nk; ++k)\
    {\
      lmn = node[k];\
      j = k-K0;\
      for (i = I0; i < Ni; ++i)\
      {\
        ijk = node[i];

#define DDM_SCHUR_JACOBIAN_BC_Epart_CLOSE \
     }/* end of for (i = I0; i < Ni; ++i) */\
    }/* end of for (k = K0; k < Nk; ++k) */\
  }/* end of if (S->NI) */


/* macro for equation */
#define DDM_SCHUR_EQ_DECLARE \
  Patch_T *const patch = vp1;\
  DDM_Schur_Complement_T *const S = vp2;\
  double *const F = S->f;\
  const Uint *const node  = S->inv;/* inverse map to node */\
  const Uint N = S->Oi;/* number of inner mesh nodes */\
  Uint n;
  
#define DDM_SCHUR_EQ_OPEN \
  for (n = 0; n < N; ++n)\
  {\
    ijk  = node[n];


#define DDM_SCHUR_EQ_CLOSE }

/* macro for boundary condition */
#define DDM_SCHUR_BC_DECLARE \
  Boundary_Condition_T *const bc = vp1;\
  DDM_Schur_Complement_T *const S = vp2;\
  double *const F      = S->f;\
  Uint *const map  = S->map;\
  Patch_T *const patch = bc->patch;\
  const Uint *const node = bc->node;/* nodes at boundary */\
  const Uint N = bc->nn;/* number of nodes at boundary */\
  Uint n;

#define DDM_SCHUR_BC_OPEN \
  for (n = 0; n < N; ++n)\
  {\
    ijk  = node[n];
    
#define DDM_SCHUR_BC_CLOSE }



/* forward declaration structures */
struct FIELD_T;
struct MATRIX_T;


/* typedef function for Solve_Equations_T */
typedef void fFunc_field_update_T(Patch_T *const patch,const char *const name);
typedef void fFunc_source_update_T(Grid_T *const grid,const char *const name);
typedef int  fFunc_stop_criteria_T(Grid_T *const grid,const char *const name);

/* a general prototype to embrace various types of equations */
typedef void *fEquation_T(void *vp1,void *vp2);

/* elements of Jacobian for equations like dfxx_df etc. */
typedef double fJs_T(struct MATRIX_T *const m,const long i,const long j);

/* equation stucture */
typedef struct sEQUATION_T
{
  char name[MAX_STR_MATH_EQ_SOLVE_LIB];
  fEquation_T *eq;/* the equation needed to be satisfied */
}sEquation_T;

/* different quantities giving info abour pairing used in Schur complement */
typedef struct PAIR_T
{
  struct SEWING_T *sewing;/* refers to its sewing */
  double *pg;/* partial g that comping from this pair*/
  SubFace_T *subface;/* differet subfaces due to patch[pn] that
                     // is related to the current patch that equations
                     // are being set up 
                     */
  struct PAIR_T *mirror;/* the pair that is mirror of itself but
                        // from the other patch. */
  Uint patchN;/* patch number which is equal to its sewing number */
  struct/* interpolation points;general coords of points
        // needed for interpolation subfaces */
  {
    double X[3];
  }*ip;
  struct/* normal vector at the subface */
  {
    double N[3];
  }*nv;
  
}Pair_T;

/* boundary information and how different patches are sown.
// this struct is made specially for having a good concurrency.
*/
typedef struct SEWING_T
{
  Pair_T **pair;
  Uint patchN;/* patch number which is equal to its sewing number */
  Uint npair;/* number of pairs */
  /* the following are the quantities that 
  // patch[patchN]->method->SchurC has.
  // it's used for purpose of concurrency and avoing race condition
  // bewteen pairs of different patches. more definition of each quantity
  // refer to SchurC strcut. */
  Uint NS;
  Uint NI;
  Uint Oi;
  Uint *map;
  Uint *inv;
  Uint *Imap;
  Uint *Iinv;
}Sewing_T;

/* ingredients needed for mapping, labeling and etc for
// domain decomposition schur complement method
*/
typedef struct DDM_SCHUR_COMPLEMENT_T
{
  struct PATCH_T *patch;/* refers to its patch itself */
  /* regular means i_j_k_to_ijk(n,i,j,k) */
  Uint *map;/* map: regular -> relabeled. ex: map[2] = 5 */
  Uint *inv;/* inv: relabeled -> regular. ex: inv[5] = 2 */
  Uint *Imap;/* interface point map, if it is given a point
                 // outside of its domain, it returns UINT_MAX. */
  Uint *Iinv;/* interface point inverse map */
  Uint NS;/* Number of subdomain points i.e. 
              // number of inner points + outerboundar points (NO) =>
              // total nodes - NS = number of interface points. Note:
              // outerboundary points excluded from interface points.
              */
  Uint NI;/* total number of interface points, if 0, it means there
              // is no interface for this patch, for example when you
              // only have one single patch, all sides of the patch
              // are outerbounday so no interface with other patches. */
  Uint Oi;/* initial index of outer boundary points at new label.
              // e.g. if NS = 10 and the last 3 points are 
              // outer boundary points then Oi = 7. 
              // furthermore, if there is no any outer boundary points 
              // then Oi = NS. */
  
/* namings:
   |B E||x|   |f|
   |F C||y| = |g|
*/
  double *f;
  double *f_prime;
  double *F_by_f_prime;
  double *g;
  double *x;
  double *y;
  struct MATRIX_T *B;
  struct MATRIX_T *E_Trans;/* NOE: this is TRANSPOSE of E */
  struct MATRIX_T *E_Trans_prime;/* NOTE: it is E' of E_Trnas. */
  struct MATRIX_T **F_by_E_prime;/* F*E' for each patch in regular format */
  struct MATRIX_T **F;
  struct MATRIX_T **C;
  struct MATRIX_T *subS;/* subS = C - F_by_E_prime_reg in ccs format */
  
  Sewing_T **sewing;/* sewing[patch_number] */
  Uint nsewing;/* number of sewings which is = number of patches */
  Uint np;/* total number of patches */
  Uint *NS_p;/* SchurC->NS for each patch p */
  Uint NS_total;/* summation of all NS_p */
  Uint *NI_p;/* SchurC->NI for each patch p */
  Uint NI_total;/* summation of all NI_p */
  
}DDM_Schur_Complement_T;

/* solving management */
typedef struct SOLVING_MAN_T
{
  struct PATCH_T *patch;/* refers to its patch itself */
  char **field_name;/* field to be solved */
  Uint nf;/* number of fields */
  Uint cf;/* current field; index of the field is being solved */
  double Frms;/* the current residual(rms) of F in, Jx=-F for this field 
              // at this patch. note: it's initialized to DBL_MAX. */
  fEquation_T **field_eq;/* the equation needed to be satisfied */
  fEquation_T **bc_eq;/* the B.C. needed to be satisfied */
  fEquation_T **jacobian_field_eq;/* jacobian for field equations */
  fEquation_T **jacobian_bc_eq;/* jacobian for B.C. equations */
  struct/* spectral jacobian elements for numeric populations */
  {
    char type[MAX_STR_MATH_EQ_SOLVE_LIB];
    struct MATRIX_T *J;/* spectral Jacobian */
  }**jacobian;
  Uint nj;/* number of jacobian */
  
  struct/* workspace for spectral jacobian */
  {
    double *dT_dx[3];/* save dCheb_Tn(n[?],x)/dx|ijk, 
                     // where n[?] = patch->n[?]-1. */
    double *d2T_dx2[3];/* save d2Cheb_Tn(n[?],x)/dx2|ijk, 
                       // where n[?] = patch->n[?]-1. */
  }jacobian_workspace[1];
  
  struct/* various method to solve */
  {
    /* type of method */
    Uint Schur_Complement: 1;/*1 if schur complement, 0 otherwise*/
    DDM_Schur_Complement_T *SchurC;
  }method[1];
  
  struct/* settings and options for solver */
  {
    double relaxation_factor;/* relaxation factor in relaxation scheme */
    double Frms_i;/* the very beginning Frms (see Frms above for definition)
                 // which this field has at the its entrance to the solver. */
    double *HFrms;/* history of all Frms start form 0 to NFrms */
    double *last_sol;/* it is back up of last solution, 
                     // so in case the residula goes up, it uses this value. */
    Uint NHFrms;/* number of HFrms */
    int solver_step;/* number of steps have been taken by solver till now. starting from 0 */
    int umfpack_size;/* 0 = di, otherwise dl (default is 0) */
    double umfpack_refine;/* max iter. refinement step, default is the default of UMFPACK which is 2 */
  }settings[1];
}Solving_Man_T;

/* equation solver */
typedef int fEquation_Solver_T(void *vp);

/* general function for variation of various kind of interpolation with 
// respect to the field. */
typedef double fdInterp_dfs_T(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);


/* boundary condition struct */
typedef struct BOUNDARY_CONDITION_T
{
  Patch_T *patch;/* patch that has this boundary */
  SubFace_T *subface;/* the subface located at interesting boundary */
  struct FIELD_T *field;/* the field this B.C.s to be imposed */
  Uint cn;/* collection number */
  Uint *node;/* nodes's index at the boundary, i.e node[i] = node number used in the patch */
  Uint nn;/* number of nodes */
}Boundary_Condition_T;


enum ROOT_FINDER_enum
{
  ROOT_FINDER_UNDEF/* undefined */,
  ROOT_FINDER_OK/* root was found successfully */,
  ROOT_FINDER_EXTREMA/* it stuck in an extrema */,
  ROOT_FINDER_MAX_ITER/* exceeds from maximum number of iteration */,
  ROOT_FINDER_NO_IMPROVEMENT/* it could not improve it more */,
  ROOT_FINDER_INTERRUPTED/* it was interrupted by a condition by user */,
  ROOT_FINDER_NAN/* the residual gets nan */
};

/* struct for root finder routine */
typedef struct ROOT_FINDER_T
{
  const char *type;/* type of root finder */
  const char *description;/* if might give some description for the root finder */
  double residual;/* residual of the function from zero */
  double tolerance;/* tolerance for f(x) = 0, 
                   // if |f^{iter+1}(x)-f^{iter}(x)| < tol, the root finder stops */
  Uint n;/* number of variables (or equations) that make f = 0, 
             // e.g in {f1(x1,x2) = 0,f2(x1,x2) = 0, n is 2 */
  Uint MaxIter;/* maximum iteration */
  Uint eq_number;/* current equation number, there are cases ,e.g. PDE, 
                     // that the equations are the same but they are 
                     // at different point, this could help to populate the
                     // root->f with one function but the funcation is evaluated
                     // at different points. */
  const double *x_gss;/* initial guess */
  double *x_sol;/* solution of f(x) = 0 */
  Uint FD_Left : 1;/* if 1 it uses finite difference with Left side stencil */
  Uint FD_Right: 1;/* if 1 it uses finite difference with Right side stencil */
  void *params;/* parameters needed for evaluation of f(x) */ 
  /* f(x1,x2,...) = 0, params is supposed to refere to whatever is needed for evaluation of f */
  // note: since it might be systems of equations like {f1=0,f2=0,...} I used pointer to pointer function */
  double (**f)(void *params,const double *const x);
  /* df/dx^{dir}, params is the parameters are used for evalution of df_dx,
  // x is the dependent variables and dir is the direction of derivative */
  double (**df_dx)(void *params,const double *const x,const Uint dir);
  double *(*root_finder_func)(struct ROOT_FINDER_T *const root);
  enum ROOT_FINDER_enum exit_status;/* exit status of root finder */
  int interrupt;/* if interrupt != 0, the root finder is interrupted.
                // for example, this controls if during search of root, 
                // root finder exceeds the domain of function. note, this 
                // must be set by the user at the equation function f(x). */
  int verbose;/* if 1, prints every step of root finding */
  double a_bisect;/* note: f(x) must change sign for x in [a,b]. */
  double b_bisect;/* note: f(x) must change sign for x in [a,b]. */
}Root_Finder_T;

/* solve equation struct that is passed to the solver.
// it may contain various functions and parameters to control and
// execute different tasks. */
typedef struct SOLVE_EQUATIONS_T
{
  Grid_T *grid;/* default grid, if no specific grid for particular
               // field has been specified, it uses this grid which
               // is given at the time of initialization. */
  const char *field_name;/* the name of the field that is being solved now */
  const char *solving_order;/* field name separated with comma to be solved,
                      // e.g. "phi,psi' means solve for first phi 
                      // and then psi. */
  const char *residual_suffix;/* residual field name suffix */                    
  double relaxation_factor;/* in relaxation scheme we have :
                           // X_new = A*X'+(A-1)X_old, where
                           // A is the relaxation_factor, 
                           // X' is the solution found by the solver. 
                           // this factor can be set for each field separately
                           // and if no info is given, it is equal to 1, which
                           // means no relaxation. */
  
  int umfpack_size;/* (0 = di) otherwise long, default is 0 */
  double umfpack_refine;/* max iter. refinement step, default is the default of UMFPACK which is 2 */
  /* some fields need their own grid, called sgrid (Special GRID) here. 
  // e.g. phi in Euler's equations is solved only in NS not the whole grid */
  struct
  {
    char *name;/* name of the field with special grid, e.g. phi */
    Grid_T *sgrid;/* e.g. the grid composed of NS patches for phi */
  }**Sgrid;/* the end of this struct determined by Null */
  
  /* instructions for updating field and its derivative according to the
  // field name particulare task for updating is done. note, if 
  // it has not been assigned it won't be execute. */
  fFunc_field_update_T(*FieldUpdate);
  
  /* instructions for updating the sources after the field has been 
  // solved on the whole grid and its derivative according to the given
  // name of the field. 
  // note, if it has not been assigned it won't be executed.*/
  fFunc_source_update_T(*SourceUpdate);
  
  /* this is the function specifies the stop criteria of the solver
  // if 1 it means continue, 0 means stop. if no function defined 
  // the default function is made using 
  // Solving_Residual and Solving_Max_Number_of_Solver_Step parameter */
  fFunc_stop_criteria_T(*StopCriteria);
}Solve_Equations_T;

void calculate_equation_residual(Solve_Equations_T *const SolveEqs);
char **get_solving_field_name(const char *const solving_order,Uint *const nf);
void print_root_finder_exit_status(const Root_Finder_T *const root);
Root_Finder_T *init_root_finder(const Uint n);
double *execute_root_finder(Root_Finder_T *const root);
void plan_root_finder(Root_Finder_T *const root);
void free_root_finder(Root_Finder_T *root);
int solve_eqs(Solve_Equations_T *const SolveEqs);
void free_solve_equations(Solve_Equations_T *solve);
double get_relaxation_factor_solve_equations(Solve_Equations_T *const solve);
Solve_Equations_T *init_solve_equations(Grid_T *const grid);
Grid_T *get_grid_solve_equations(Solve_Equations_T *const solve);
void add_special_grid_solve_equations(Grid_T *const grid,const char *const name, Solve_Equations_T *const solve);
void make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void test_make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void test_dfs_df_spectral_vs_FiniteDiff(Grid_T *const grid);
void test_dfs_df_Spectral_vs_analytic(Grid_T *const grid);
void test_dInterp_a_df(Grid_T *const grid);
void *init_eq(void);
void add_eq(sEquation_T ***const data_base, fEquation_T *const eq,const char *const name);

void initialize_solving_man(Grid_T *const grid,
                            sEquation_T **const field_eq,
                            sEquation_T **const bc_eq,
                            sEquation_T **const jacobian_field_eq,
                            sEquation_T **const jacobian_bc_eq,
                            const char *const prefix);

struct MATRIX_T *get_j_matrix(const Patch_T *const patch,const char *type);
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
double read_matrix_entry_ccs(struct MATRIX_T *const m, const long r,const long c);
fJs_T *get_j_reader(const struct MATRIX_T *const m);
void test_Jacobian_of_equations(Solve_Equations_T *const SolveEqs);
void test_root_finders(Grid_T *const grid);
fdInterp_dfs_T *get_dInterp_df(const Patch_T *const patch,const SubFace_T *const sf,const char *const dir);
Sewing_T *alloc_sewing(void);
void free_db_eqs(sEquation_T **db);
void free_patch_SolMan_jacobian(Patch_T *const patch);
void free_patch_SolMan_method_Schur(Patch_T *const patch);
void move_dfdu_jacobian_patch(Patch_T *const patch2,Patch_T *const patch1);

double
  d2f_dxdu_spectral_Jacobian_analytic(Patch_T *const patch,
                                      const Uint dx_axis, 
                                      const Uint ijk,const Uint lmn);
  
double
  d3f_dxdydu_spectral_Jacobian_analytic(Patch_T *const patch,
                                        const int dxdy_axis,
                                        const Uint ijk,const Uint lmn);

#endif



