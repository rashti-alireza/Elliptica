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

/* basic trig math */
#define Cos(a) cos((a))
#define Sin(a) sin((a))
#define Csc(a) (1./sin((a)))
#define Cot(a) (1./tan((a)))

/* Kronecker Delta */
#define K__D(i,j)  ( (i)==(j) ? 1. : 0.)

/* (-1)^n */
#define SIGN(n) (((n)%2) ? -1. : 1.)

/* quick theta */
#define THETA(axis,ijk) (patch->node[(ijk)]->theta[(axis)])

/* normalization */
#define NORMALIZATION(n) ( 1./( 2.*((n)-1.) ) )

/* sum_{n=0}^{N} cos(n lambda) = 
// 0.5 + 0.5*( sin( (N+0.5)*(lambda) ) ) / ( sin( 0.5*(lambda) ) ),
// N0 = N+0.5. */
#define sum_0_N_cos_nlambda(N,N0,lambda) \
  ( \
    0.5 + 0.5*( sin( (N0)*(lambda) ) ) / ( sin( 0.5*(lambda) ) ) \
  )

/* d/dlambda sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define d_dlambda_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) ?\
    (0.0):\
    ( \
      Csc(0.5*(lambda))*(2.*(N0)*Cos((lambda)*(N0)) - \
      Cot(0.5*(lambda))*Sin((lambda)*(N0)))\
    )*0.25\
  )

/* d^2/dlambda^2 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define d2_dlambda2_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) ?\
    -(Pow3(N)/3.+Pow2(N)/2.+N/6.)/* simplified, don't forget - sign! */ :\
    (\
      Csc(0.5*(lambda))*(-4.*(N0)*Cos((lambda)*(N0))*Cot(0.5*(lambda)) + \
      (-1. - 4.*Pow2(N0) + 2.*Pow2(Csc(0.5*(lambda))))*Sin((lambda)*(N0)))\
    )*0.125\
  )

/* d^3/dlambda^3 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define d3_dlambda3_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) ?\
    (0.0):\
    (\
      pow(Csc((lambda)/2.,3))*(2*(N0)*\
        (9 - 4*Pow2((N0)) + (3 + 4*Pow2((N0)))*Cos((lambda)))*\
        Cos((lambda)*(N0)) - (11 - 12*Pow2((N0)) + Cos((lambda)) + \
          12*Pow2((N0))*Cos((lambda)))*Cot((lambda)/2.)*Sin((lambda)*(N0)))\
    )/32.\
  )

/* d^4/dlambda^4 sum_{n=0}^{N} cos(n lambda).
// N0 = N+0.5. */
#define d4_dlambda4_sum_0_N_cos_nlambda(N,N0,lambda) \
  ( EQL((lambda),0.) ?\
    (Pow4(N)*(N/5.+0.5)+Pow3(N)/3.-N/30.)/* simplified */:\
    (\
      pow(Csc((lambda)/2.),5)*(-16*(N0)*\
        (11 - 4*Pow2((N0)) + Cos((lambda)) + \
          4*Pow2((N0))*Cos((lambda)))*Cos((lambda)*(N0))*Sin((lambda)) + \
       (115 - 120*Pow2((N0)) + 48*Pow4((N0)) + \
          (76 + 96*Pow2((N0)) - 64*Pow4((N0)))*Cos((lambda)) + \
          (1 + 24*Pow2((N0)) + 16*Pow4((N0)))*Cos(2*(lambda)))*\
        Sin((lambda)*(N0)))\
    )/256.\
  )


/* -> d/dX(df/du)=d/dX (2*sum_0^N (Tn(Xj) Tn(X)) -1 -(-1)^j *T_{N-1}(X)),
// NOTE: X = cos(th) */
#define JACOBIAN_d_dX_df_du(thi,thj,N,j) \
   ( d_dXi_2xsum_0_N_Tnj_Tni(thi,thj,N) - SIGN(j)*dT_dx((int)(N)-1,cos(thi)) )

/* -> d2/dX^2(df/du)=d2/dX^2 (2*sum_0^N (Tn(Xj) Tn(X)) -1 -(-1)^j *T_{N-1}(X)),
// NOTE: X = cos(th)). */
#define JACOBIAN_d2_dX2_df_du(thi,thj,N,j) \
   ( d2_dXi2_2xsum_0_N_Tnj_Tni(thi,thj,N) - SIGN(j)*d2T_dx2((int)(N)-1,cos(thi)) )

/* normalization * coords jacobian * JACOBIAN_d_dX_df_du */
#define JACOBIAN_dX_dx_d_dX_df_du(patch, x_axis, X_axis, ijk, lmn, pn/* 1-d point number */) \
  ( NORMALIZATION(patch->n[X_axis])*dX_dx[ijk][X_axis][x_axis]*dN_dX[X_axis]*\
    JACOBIAN_d_dX_df_du(THETA(X_axis,ijk),THETA(X_axis,lmn),patch->n[X_axis]-1,pn) )

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
static double dT_dx(const int n,const double x);
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


