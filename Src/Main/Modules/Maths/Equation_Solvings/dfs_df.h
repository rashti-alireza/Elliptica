#include "core_lib.h"
#include "utilities_lib.h"
#include "macros_lib.h"
#include "manifold_lib.h"
#include "error_handling_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "maths_linear_algebra_lib.h"
#include "fields_lib.h"
#include "maths_equation_solvings_lib.h"

#define _MAX_STR_ (400)

/* Define macros for analytic calculation of df/du Jacobian */

/* a short hand notation */
#define JW (patch->solving_man->jacobian_workspace)

/* basic trig math */
#define JCos(a) cos((a))
#define JSin(a) sin((a))
#define JCsc(a) (1./sin((a)))
#define JCot(a) (1./tan((a)))

/* 2*M_PI */
#define J2M_PI (6.283185307179586)

/* Kronecker Delta */
#define JKD(i,j)  ( (i)==(j) ? 1. : 0.)

/* eta_i in Elliptica's paper, NOTE: n = patch->n[?] - 1 */
#define Jeta(i,n) ( (i) == 0 || (i) == (n) ? 1. : 2. )

/* (-1)^i */
#define Jsign(i) (((i)%2) ? -1. : 1.)

/* quick theta: 
// NOTE1: assuming Chebyshev Extrema points.
// NOTE2: assuming patch is defined. */
#define Jtheta(i,X_axis) ( (i)*M_PI/((patch->n[(X_axis)])-1.) )

/* normalization, NOTE: n = patch->n */
#define Jnorm(n) ( 0.5/((n)-1.) )

/* dX/dx */
#define JdX_dx(patch,ijk,dX_axis,dx_axis) \
  ( (patch)->JacobianT->dX_dx[(dX_axis)][(dx_axis)][(ijk)] )

/* d2X/dxdy */
#define Jd2X_dxdy(patch,ijk,dX_axis,dxdy_axis) \
  ( (patch)->JacobianT->d2X_dxdy[(dX_axis)][(dxdy_axis)][(ijk)] )

/* dN/dX */
#define JdN_dX(patch,dX_axis) \
  ( (patch)->JacobianT->dN_dX[(dX_axis)] )

/* sum_{n=0}^{N} cos(n lambda) = 
// 0.5 + 0.5*( sin( (N+0.5)*(lambda) ) ) / ( sin( 0.5*(lambda) ) ),
// N0 = N+0.5. */
#define Jsum_0_N_cos_nlambda(N,N0,lambda) \
  ( \
    0.5 + 0.5*( sin( (N0)*(lambda) ) ) / ( sin( 0.5*(lambda) ) ) \
  )

/* d/dlambda sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define Jd_dlambda_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J2M_PI) ?\
    (0.0):\
    ( \
      JCsc(0.5*(lambda))*(2.*(N0)*JCos((lambda)*(N0)) - \
      JCot(0.5*(lambda))*JSin((lambda)*(N0)))\
    )*0.25\
  )

/* d^2/dlambda^2 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define Jd2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J2M_PI) ?\
    -(Pow3(N)/3.+Pow2(N)/2.+N/6.)/* simplified, don't forget - sign! */ :\
    (\
      JCsc(0.5*(lambda))*(-4.*(N0)*JCos((lambda)*(N0))*JCot(0.5*(lambda)) + \
      (-1. - 4.*Pow2(N0) + 2.*Pow2(JCsc(0.5*(lambda))))*JSin((lambda)*(N0)))\
    )*0.125\
  )

/* d^3/dlambda^3 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define Jd3_dlambda3_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J2M_PI) ?\
    (0.0):\
    (\
      pow(JCsc(0.5*(lambda),3))*(2*(N0)*\
        (9 - 4*Pow2((N0)) + (3 + 4*Pow2((N0)))*JCos((lambda)))*\
        JCos((lambda)*(N0)) - (11 - 12*Pow2((N0)) + JCos((lambda)) + \
          12*Pow2((N0))*JCos((lambda)))*JCot(0.5*(lambda))*JSin((lambda)*(N0)))\
    )/32.\
  )

/* d^4/dlambda^4 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define Jd4_dlambda4_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J2M_PI) ?\
    (Pow4(N)*(N/5.+0.5)+Pow3(N)/3.-N/30.)/* simplified */:\
    (\
      pow(JCsc(0.5*(lambda)),5)*(-16*(N0)*\
        (11 - 4*Pow2((N0)) + JCos((lambda)) + \
          4*Pow2((N0))*JCos((lambda)))*JCos((lambda)*(N0))*JSin((lambda)) + \
       (115 - 120*Pow2((N0)) + 48*Pow4((N0)) + \
          (76 + 96*Pow2((N0)) - 64*Pow4((N0)))*JCos((lambda)) + \
          (1 + 24*Pow2((N0)) + 16*Pow4((N0)))*JCos(2*(lambda)))*\
        JSin((lambda)*(N0)))\
    )/256.\
  )

/* -> eta_j{d/dX(df/du)=d/dX (2*sum_0^N (Tn(Xj) Tn(X)) -1 -(-1)^j *T_{N}(X))},
// NOTE: X = cos(th), N = patch->n[?]-1. */
#define Jd2f_dudX(thi,thj,N,j) \
   (Jeta(j,N)*( d_dXi_2xsum_0_N_Tnj_Tni(thi,thj,N) - Jsign(j)*dCheb_Tn_dx((int)(N),cos(thi)) ))

/* -> d2/dX^2(df/du)=d2/dX^2 (2*sum_0^N (Tn(Xj) Tn(X)) -1 -(-1)^j *T_{N}(X)),
// NOTE: X = cos(th)), N = patch->n[?]-1. */
#define Jd3f_dudXdX(thi,thj,N,j) \
   (Jeta(j,N)*( d2_dXi2_2xsum_0_N_Tnj_Tni(thi,thj,N) - Jsign(j)*d2Cheb_Tn_dx2((int)(N),cos(thi)) ))

/* normalization * coords jacobian * Jd2f_dudX */
#define Jd2f_dudx(patch, dx_axis, X_axis, ijk, lmn, qi,qj) \
  ( Jnorm(patch->n[X_axis])*JdX_dx(patch,ijk,X_axis,dx_axis)*JdN_dX(patch,X_axis)*\
    Jd2f_dudX(Jtheta(qi,X_axis),Jtheta(qj,X_axis),patch->n[X_axis]-1,qj) )

/* normalization * coords jacobian * Jd3f_dudXdX */
#define Jd3f_dudxdy(patch, dx_axis, dy_axis, dxdy_axis,X_axis, ijk, lmn, qi,qj) \
  ( \
    Jnorm(patch->n[X_axis])*JdN_dX(patch,X_axis)*\
    ( \
      Jd2X_dxdy(patch,ijk,X_axis,dxdy_axis)*\
      Jd2f_dudX(Jtheta(qi,X_axis),Jtheta(qj,X_axis),patch->n[X_axis]-1,qj) +     \
      JdX_dx(patch,ijk,X_axis,dx_axis)*JdX_dx(patch,ijk,X_axis,dy_axis)*JdN_dX(patch,X_axis)* \
      Jd3f_dudXdX(Jtheta(qi,X_axis),Jtheta(qj,X_axis),patch->n[X_axis]-1,qj)     \
    )\
  )


/* quick theta (optimized): 
// NOTE1: assuming Chebyshev Extrema points.
// NOTE2: assuming patch is defined. */
#define Jtheta_opt(i,X_axis) ( (i)*(JW->pi_o_nm1[X_axis]) )

/* normalization, optimized */
#define Jnorm_opt(X) ( JW->norm[X] )

/* d/dlambda sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. (optimized) */
#define Jd_dlambda_sum_0_N_cos_nlambda_opt(X_axis,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J2M_PI) ?\
    (0.0):\
    ( \
      csc_half_lambda*(2.*(N0)*JCos((lambda)*(N0)) - \
      cot_half_lambda*JSin((lambda)*(N0)))\
    )*0.25\
  )

/* d^2/dlambda^2 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define Jd2_dlambda2_sum_0_N_cos_nlambda_opt(X_axis,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J2M_PI) ?\
    JW->c1_d2[X_axis]:\
    (\
      csc_half_lambda*(-4.*(N0)*JCos((lambda)*(N0))*cot_half_lambda + \
      (JW->c2_d2[X_axis] + 2.*Pow2(csc_half_lambda))*JSin((lambda)*(N0)))\
    )*0.125\
  )

/* d^4/dlambda^4 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5, (optimized). */
#define Jd4_dlambda4_sum_0_N_cos_nlambda_opt(X_axis,N0,lambda) \
  ( EQL((lambda),0.) || EQL((lambda),J2M_PI) ?\
    JW->c1_d4[X_axis]:\
    (\
      Pow2(csc_half_lambda)*Pow3(csc_half_lambda)*(-16*(N0)*\
        ( 11 - JW->c2_d4[X_axis] + cos_lambda + \
          JW->c2_d4[X_axis]*cos_lambda)*JCos((lambda)*(N0))*JSin((lambda)) + \
       (  JW->c3_d4[X_axis] + \
          JW->c4_d4[X_axis]*cos_lambda + \
          JW->c5_d4[X_axis]*JCos(2*(lambda)))*\
        JSin((lambda)*(N0)))\
    )/256.\
  )

/* -> eta_j{d/dX(df/du)=d/dX (2*sum_0^N (Tn(Xj) Tn(X)) -1 -(-1)^j *T_{N}(X))},
// NOTE: X = cos(th), N = patch->n[?]-1 (optimized). */
#define Jd2f_dudX_opt(thi,thj,N,patch,j,X_axis,ijk) \
   (Jeta(j,N)*( d_dXi_2xsum_0_N_Tnj_Tni_opt(thi,thj,patch,X_axis) - Jsign(j)*(JW->dT_dx[X_axis][ijk]) ))

/* -> d2/dX^2(df/du)=d2/dX^2 (2*sum_0^N (Tn(Xj) Tn(X)) -1 -(-1)^j *T_{N}(X)),
// NOTE: X = cos(th)), N = patch->n[?]-1, (optimized). */
#define Jd3f_dudXdX_opt(thi,thj,N,patch,j,X_axis,ijk) \
   (Jeta(j,N)*( d2_dXi2_2xsum_0_N_Tnj_Tni_opt(thi,thj,patch,X_axis) - Jsign(j)*(JW->d2T_dx2[X_axis][ijk]) ))

/* normalization * coords jacobian * Jd2f_dudX (optimized) */
#define Jd2f_dudx_opt(patch, dx_axis, X_axis, ijk, qi,qj) \
  ( Jnorm_opt(X_axis)*JdX_dx(patch,ijk,X_axis,dx_axis)*JdN_dX(patch,X_axis)*\
    Jd2f_dudX_opt(Jtheta_opt(qi,X_axis),Jtheta_opt(qj,X_axis),(JW->nm1[X_axis]),patch,qj,X_axis,ijk) )

/* normalization * coords jacobian * Jd3f_dudXdX */
#define Jd3f_dudxdy_opt(patch, dx_axis, dy_axis, dxdy_axis,X_axis, ijk, qi,qj) \
  ( \
    Jnorm_opt(X_axis)*JdN_dX(patch,X_axis)*\
    ( \
      Jd2X_dxdy(patch,ijk,X_axis,dxdy_axis)*\
      Jd2f_dudX_opt(Jtheta(qi,X_axis),Jtheta(qj,X_axis),(JW->nm1[X_axis]),patch,qj,X_axis,ijk) +     \
      JdX_dx(patch,ijk,X_axis,dx_axis)*JdX_dx(patch,ijk,X_axis,dy_axis)*JdN_dX(patch,X_axis)* \
      Jd3f_dudXdX_opt(Jtheta_opt(qi,X_axis),Jtheta_opt(qj,X_axis),(JW->nm1[X_axis]),patch,qj,X_axis,ijk)     \
    )\
  )



/* Jacobian type */
typedef enum JTYPE_E
{
  T_x/* df_x/df */,
  T_xx/* df_xx/df */,
  T_xy/* df_xy/df */,
  T_xz/* df_xz/df */,
  T_y/* df_y/df */,
  T_yy/* df_yy/df */,
  T_yz/* df_yz/df */,
  T_z/* df_z/df */,
  T_zz/* df_zz/df */,
  T_UNDEF
}JType_E;

/* Jacobain for equation */
typedef void Js_Jacobian_eq_F(double **const J,Patch_T *const patch,JType_E jt_e);

void free_patch_SolMan_jacobian(Patch_T *const patch);
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types);
void make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void test_make_Js_jacobian_eq(Grid_T *const grid, const char * const* types);
void obsolete_fill_jacobian_spectral_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e);
Matrix_T *get_j_matrix(const Patch_T *const patch,const char *type);
double read_matrix_entry_ccs(Matrix_T *const m, const long r,const long c);
fdInterp_dfs_T *get_dInterp_df(const Patch_T *const patch,const SubFace_T *const sf,const char *const dir);
static JType_E str2JType_E(const char *const str);
static char *interpret_type(const char *const type);
static void JType_E2str(const JType_E e,char *const str);
static void make_jacobian_spectral_method(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_spectral_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_spectral_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir);
static void make_jacobian_direct_method(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_direct_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e);
static void fill_jacobian_direct_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir);
static double ChebExtrema_1point(const Uint n, const Uint p);
static double dc_df(const Uint n,const Uint i,const Uint l);
static void read_1st_and_2nd_deriv(const JType_E deriv_dir,JType_E *const deriv_1st,JType_E *const deriv_2nd);
static void JType_E2Dd_T(const JType_E jt_e, Dd_T *const q_dir);
static void write_J_in_disk_ccs(void);
static double J_sizeMb_ccs(const Matrix_T *const m);
fJs_T *get_j_reader(const Matrix_T *const m);
static void coarse_grain_Ap_ccs_matrix(Matrix_T *const m,const int Nslice) __attribute__((unused)) ;
static double dInterp_x_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_y_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_z_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_x_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_y_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_z_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_x_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_y_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_z_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_x_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_y_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_z_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);
static double dInterp_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane);

static double
d2_dXi2_2xsum_0_N_Tnj_Tni(double thi/* X_i = cos(theta_i) */,
                          double thj/* X_i = cos(theta_i) */,
                          Uint N/* the sum upper limit */);
static double
d_dXi_2xsum_0_N_Tnj_Tni(double thi/* X_i = cos(theta_i) */,
                        double thj/* X_i = cos(theta_i) */,
                        Uint N/* the sum upper limit */);

double
  d2f_dxdu_spectral_Jacobian_analytic(Patch_T *const patch,
                                      const Uint dx_axis, 
                                      const Uint ijk,const Uint lmn);
  
double
  d3f_dxdydu_spectral_Jacobian_analytic(Patch_T *const patch,
                                        const int dxdy_axis,
                                        const Uint ijk,const Uint lmn);

void set_Solving_Man_jacobian_workspace(Patch_T *const patch);



static double
d_dXi_2xsum_0_N_Tnj_Tni_opt(const double thi/* X_i = cos(theta_i) */,
                            const double thj/* X_j = cos(theta_j) */,
                            Patch_T *const patch,
                            Uint X_axis);

static double
d2_dXi2_2xsum_0_N_Tnj_Tni_opt(const double thi/* X_i = cos(theta_i) */,
                              const double thj/* X_i = cos(theta_i) */,
                              Patch_T *const patch,
                              Uint X_axis);
