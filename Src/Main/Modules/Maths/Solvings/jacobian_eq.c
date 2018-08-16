/*
// Alireza Rashti
// August 2018
*/

#include "jacobian_eq.h"
#define MAX_STR_LEN 100

/* making Jacobian for equations at the inner mesh
// types are pointers to string determining the type of jacobian
// e.g. *types[3] = {"J_xx","J_y",0}.Note: the number of
// types is found by null.
*/
void make_jacobian_eq(Grid_T *const grid, char **const types)
{
  Jacobian_eq_F *Jacobian;
  double **J = 0;
  JType_E jt_e;
  unsigned i,p,nn;
  
  /* selecting Jacobian method for making of jacobian equation */
  if (strcmp_i(GetParameterS_E("Making_Jacobian_Eq_Method"),"spectral"))
    Jacobian = make_jacobian_spectral_method;
  else if (strcmp_i(GetParameterS_E("Making_Jacobian_Eq_Method"),"direct"))
    Jacobian = make_jacobian_direct_method;
  else
    abortEr(INCOMPLETE_FUNC);
  
  i = 0;
  while (types[i] != 0)
  {
    jt_e = str2enum(types[i]);
    
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      nn = total_nodes_patch(patch);
      J = alloc_matrix(nn,nn);
      Jacobian(J,patch,jt_e);
      
      free_matrix(J,nn);
    }
  }
}

/* translating string to enum JType_E */
static JType_E str2enum(const char *const str)
{
  JType_E jt_e = T_UNDEF;
  
  if (strcmp_i(str,"J_x"))
    jt_e = T_x;
  else if (strcmp_i(str,"J_xx"))
    jt_e = T_xx;
  else if (strcmp_i(str,"J_y"))
    jt_e = T_y;
  else if (strcmp_i(str,"J_yy"))
    jt_e = T_yy;
  else if (strcmp_i(str,"J_z"))
    jt_e = T_z;
  else if (strcmp_i(str,"J_zz"))
    jt_e = T_zz;
  else
    abortEr(INCOMPLETE_FUNC);
  
  return jt_e;
}

/* making Jacobian equations using direct method */
static void make_jacobian_direct_method(double **const J,Patch_T *const patch,JType_E jt_e)
{
  switch(jt_e)
  {
    case T_x:
      fill_jacobian_direct_method_x(J,patch);
      break;
    case T_xx:
      fill_jacobian_direct_method_xx(J,patch);
      break;
    case T_y:
      fill_jacobian_direct_method_y(J,patch);
      break;
    case T_yy:
      fill_jacobian_direct_method_yy(J,patch);
      break;
    case T_z:
      fill_jacobian_direct_method_z(J,patch);
      break;
    case T_zz:
      fill_jacobian_direct_method_zz(J,patch);
      break;
    default:
      abortEr("No such type for Jacobian defined!\n");
  }
}

/* making Jacobian using direct method in direction x 
// d(df(i,j,k)/dx)/df(l,m,n) = (d(f+df)/dx-df/dx)/df
*/
static void fill_jacobian_direct_method_x(double **const J, Patch_T *const patch)
{
  const unsigned num_thread = (unsigned)omp_get_num_threads();
  Field_T **j;
  const double EPS = 1E-7;
  const double CONST = 1.0;
  const unsigned nn = total_nodes_patch(patch);
  unsigned lmn,tn;
  char name[MAX_STR_LEN];
  
  j = malloc(num_thread*sizeof(*j));
  pointerEr(j);
  
  /* adding a field for each thread */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    j[tn] = add_field(name,"(3dim)",patch,YES);
    
    for (lmn = 0; lmn < nn; ++lmn)
      j[tn]->v[lmn] = CONST;
  }
  
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Field_T *Jf = j[omp_get_thread_num()];
    double *J_x = 0;
    unsigned ijk;
    
    Jf->v[lmn] += EPS;   
    J_x = Partial_Derivative(Jf,"x");
    /* since it was added v2 and info in J_x we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_info(Jf);
    free_v2(Jf);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_x[ijk]/EPS;
      
    free(J_x);
    Jf->v[lmn] -= EPS;
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    Field_T *f = patch->pool[Ind(name)];
    remove_field(f);
  }
  free(j);
}

/* making Jacobian using direct method in direction xx 
// d(d^2f(i,j,k)/dx^2)/df(l,m,n) = (d^2(f+df)/dx^2-d^2f/dx^2)/df
*/
static void fill_jacobian_direct_method_xx(double **const J, Patch_T *const patch)
{
  const unsigned num_thread = (unsigned)omp_get_num_threads();
  Field_T **j;
  const double EPS = 1E-7;
  const double CONST = 1.0;
  const unsigned nn = total_nodes_patch(patch);
  unsigned lmn,tn;
  char name[MAX_STR_LEN];
  
  j = malloc(num_thread*sizeof(*j));
  pointerEr(j);
  
  /* adding a field for each thread */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    j[tn] = add_field(name,"(3dim)",patch,YES);
    
    for (lmn = 0; lmn < nn; ++lmn)
      j[tn]->v[lmn] = CONST;
  }
  
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Field_T *Jf = j[omp_get_thread_num()];
    double *J_xx = 0;
    unsigned ijk;
    
    Jf->v[lmn] += EPS;   
    J_xx = Partial_Derivative(Jf,"xx");
    /* since it was added v2 and info in J_xx we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_info(Jf);
    free_v2(Jf);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_xx[ijk]/EPS;
      
    free(J_xx);
    Jf->v[lmn] -= EPS;
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    Field_T *f = patch->pool[Ind(name)];
    remove_field(f);
  }
  free(j);
}

/* making Jacobian using direct method in direction y
// d(df(i,j,k)/dy)/df(l,m,n) = (d(f+df)/dy-df/dy)/df
*/
static void fill_jacobian_direct_method_y(double **const J, Patch_T *const patch)
{
  const unsigned num_thread = (unsigned)omp_get_num_threads();
  Field_T **j;
  const double EPS = 1E-7;
  const double CONST = 1.0;
  const unsigned nn = total_nodes_patch(patch);
  unsigned lmn,tn;
  char name[MAX_STR_LEN];
  
  j = malloc(num_thread*sizeof(*j));
  pointerEr(j);
  
  /* adding a field for each thread */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    j[tn] = add_field(name,"(3dim)",patch,YES);
    
    for (lmn = 0; lmn < nn; ++lmn)
      j[tn]->v[lmn] = CONST;
  }
  
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Field_T *Jf = j[omp_get_thread_num()];
    double *J_y = 0;
    unsigned ijk;
    
    Jf->v[lmn] += EPS;   
    J_y = Partial_Derivative(Jf,"y");
    /* since it was added v2 and info in J_y we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_info(Jf);
    free_v2(Jf);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_y[ijk]/EPS;
      
    free(J_y);
    Jf->v[lmn] -= EPS;
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    Field_T *f = patch->pool[Ind(name)];
    remove_field(f);
  }
  free(j);
}

/* making Jacobian using direct method in direction yy
// d(d^2f(i,j,k)/dy^2)/df(l,m,n) = (d^2(f+df)/dy^2-d^2f/dy^2)/df
*/
static void fill_jacobian_direct_method_yy(double **const J, Patch_T *const patch)
{
  const unsigned num_thread = (unsigned)omp_get_num_threads();
  Field_T **j;
  const double EPS = 1E-7;
  const double CONST = 1.0;
  const unsigned nn = total_nodes_patch(patch);
  unsigned lmn,tn;
  char name[MAX_STR_LEN];
  
  j = malloc(num_thread*sizeof(*j));
  pointerEr(j);
  
  /* adding a field for each thread */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    j[tn] = add_field(name,"(3dim)",patch,YES);
    
    for (lmn = 0; lmn < nn; ++lmn)
      j[tn]->v[lmn] = CONST;
  }
  
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Field_T *Jf = j[omp_get_thread_num()];
    double *J_yy = 0;
    unsigned ijk;
    
    Jf->v[lmn] += EPS;   
    J_yy = Partial_Derivative(Jf,"yy");
    /* since it was added v2 and info in J_yy we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_info(Jf);
    free_v2(Jf);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_yy[ijk]/EPS;
      
    free(J_yy);
    Jf->v[lmn] -= EPS;
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    Field_T *f = patch->pool[Ind(name)];
    remove_field(f);
  }
  free(j);
}

/* making Jacobian using direct method in direction z
// d(df(i,j,k)/dz)/df(l,m,n) = (d(f+df)/dz-df/dz)/df
*/
static void fill_jacobian_direct_method_z(double **const J, Patch_T *const patch)
{
  const unsigned num_thread = (unsigned)omp_get_num_threads();
  Field_T **j;
  const double EPS = 1E-7;
  const double CONST = 1.0;
  const unsigned nn = total_nodes_patch(patch);
  unsigned lmn,tn;
  char name[MAX_STR_LEN];
  
  j = malloc(num_thread*sizeof(*j));
  pointerEr(j);
  
  /* adding a field for each thread */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    j[tn] = add_field(name,"(3dim)",patch,YES);
    
    for (lmn = 0; lmn < nn; ++lmn)
      j[tn]->v[lmn] = CONST;
  }
  
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Field_T *Jf = j[omp_get_thread_num()];
    double *J_z = 0;
    unsigned ijk;
    
    Jf->v[lmn] += EPS;   
    J_z = Partial_Derivative(Jf,"z");
    /* since it was added v2 and info in J_z we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_info(Jf);
    free_v2(Jf);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_z[ijk]/EPS;
      
    free(J_z);
    Jf->v[lmn] -= EPS;
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    Field_T *f = patch->pool[Ind(name)];
    remove_field(f);
  }
  free(j);
}

/* making Jacobian using direct method in direction x 
// d(d^2f(i,j,k)/dz^2)/df(l,m,n) = (d^2(f+df)/dz^2-d^2f/dz^2)/df
*/
static void fill_jacobian_direct_method_zz(double **const J, Patch_T *const patch)
{
  const unsigned num_thread = (unsigned)omp_get_num_threads();
  Field_T **j;
  const double EPS = 1E-7;
  const double CONST = 1.0;
  const unsigned nn = total_nodes_patch(patch);
  unsigned lmn,tn;
  char name[MAX_STR_LEN];
  
  j = malloc(num_thread*sizeof(*j));
  pointerEr(j);
  
  /* adding a field for each thread */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    j[tn] = add_field(name,"(3dim)",patch,YES);
    
    for (lmn = 0; lmn < nn; ++lmn)
      j[tn]->v[lmn] = CONST;
  }
  
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Field_T *Jf = j[omp_get_thread_num()];
    double *J_zz = 0;
    unsigned ijk;
    
    Jf->v[lmn] += EPS;   
    J_zz = Partial_Derivative(Jf,"zz");
    /* since it was added v2 and info in J_zz we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_info(Jf);
    free_v2(Jf);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_zz[ijk]/EPS;
      
    free(J_zz);
    Jf->v[lmn] -= EPS;
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    Field_T *f = patch->pool[Ind(name)];
    remove_field(f);
  }
  free(j);
}

/* making Jacobian equations using spectral method */
static void make_jacobian_spectral_method(double **const J,Patch_T *const patch,JType_E jt_e)
{
  switch(jt_e)
  {
    case T_x:
      fill_jacobian_spectral_method_x(J,patch);
      break;
    case T_xx:
      fill_jacobian_spectral_method_xx(J,patch);
      break;
    case T_y:
      fill_jacobian_spectral_method_y(J,patch);
      break;
    case T_yy:
      fill_jacobian_spectral_method_yy(J,patch);
      break;
    case T_z:
      fill_jacobian_spectral_method_z(J,patch);
      break;
    case T_zz:
      fill_jacobian_spectral_method_zz(J,patch);
      break;
    default:
      abortEr("No such type for Jacobian defined!\n");
  }
}

/* making Jacobian using spectral method in direction x 
// d(df(i,j,k)/dx)/df(l,m,n) = j(N_i,x) * \sum_{ip,jp,kp} dc(ip,jp,kp)/df(l,m,n)*dT(i,j,k)/dN_i
*/
static void fill_jacobian_spectral_method_x(double **const J, Patch_T *const patch)
{
  const unsigned nn = total_nodes_patch(patch);
  const unsigned *const N = patch->n;
  unsigned ijk;
  
  OpenMP_1d_Pragma(omp parallel for)
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double cj0 = dq2_dq1(patch,_N0_,_x_,ijk);/* coordinate jacobian */
    double cj1 = dq2_dq1(patch,_N1_,_x_,ijk);/* coordinate jacobian */
    double cj2 = dq2_dq1(patch,_N2_,_x_,ijk);/* coordinate jacobian */
    double x,y,z;
    unsigned lmn;
    unsigned i,j,k;
    
    IJK(ijk,N,&i,&j,&k);
    x = ChebExtrema_1point(N[0],i);
    y = ChebExtrema_1point(N[1],j);
    z = ChebExtrema_1point(N[2],k);
    
    for (lmn = 0; lmn < nn; ++lmn)
    {
      double dc_df;
      double j0,j1,j2;
      unsigned l,m,n,ip,jp,kp;
      
      IJK(lmn,N,&l,&m,&n);
      j0 = 0;
      j1 = 0;
      j2 = 0;
      
      if (m == j && n == k)
      {
        if (!EQL(cj0,0))
        {
          for (ip = 1; ip < N[0]-1; ++ip)
          {
            dc_df = dc0_df(N[0],ip,l);
            j0 += dc_df*dT_dx((int)ip,x);
          }
          j0 *= cj0;
        }
      }
      
      if (l == i && n == k)
      {
        if (!EQL(cj1,0))
        {
          for (jp = 1; jp < N[1]-1; ++jp)
          {
            dc_df = dc1_df(N[1],jp,m);
            j1 += dc_df*dT_dx((int)jp,y);
          }
          j1 *= cj1;
        }
      }
      
      if (l == i && m == j)
      {
        if (!EQL(cj2,0))
        {
          for (kp = 1; kp < N[2]-1; ++kp)
          {
            dc_df = dc2_df(N[2],kp,n);
            j2 += dc_df*dT_dx((int)kp,z);
          }
          j2 *= cj2;
        }
      }
      
      J[ijk][lmn] = j0+j1+j2;  
    }/* end of for (lmn = 0; lmn < nn; ++lmn) */
    
  }/* end of for (ijk = 0; ijk < nn; ++ijk) */
  
}

/* making Jacobian using spectral method in direction y
// d(df(i,j,k)/dy)/df(l,m,n) = j(N_i,y) * \sum_{ip,jp,kp} dc(ip,jp,kp)/df(l,m,n)*dT(i,j,k)/dN_i
*/
static void fill_jacobian_spectral_method_y(double **const J,Patch_T *const patch)
{
  const unsigned nn = total_nodes_patch(patch);
  const unsigned *const N = patch->n;
  unsigned ijk;
  
  OpenMP_1d_Pragma(omp parallel for)
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double cj0 = dq2_dq1(patch,_N0_,_y_,ijk);/* coordinate jacobian */
    double cj1 = dq2_dq1(patch,_N1_,_y_,ijk);/* coordinate jacobian */
    double cj2 = dq2_dq1(patch,_N2_,_y_,ijk);/* coordinate jacobian */
    double x,y,z;
    unsigned lmn;
    unsigned i,j,k;
    
    IJK(ijk,N,&i,&j,&k);
    x = ChebExtrema_1point(N[0],i);
    y = ChebExtrema_1point(N[1],j);
    z = ChebExtrema_1point(N[2],k);
    
    for (lmn = 0; lmn < nn; ++lmn)
    {
      double dc_df;
      double j0,j1,j2;
      unsigned l,m,n,ip,jp,kp;
      
      IJK(lmn,N,&l,&m,&n);
      j0 = 0;
      j1 = 0;
      j2 = 0;
      
      if (m == j && n == k)
      {
        if (!EQL(cj0,0))
        {
          for (ip = 1; ip < N[0]-1; ++ip)
          {
            dc_df = dc0_df(N[0],ip,l);
            j0 += dc_df*dT_dx((int)ip,x);
          }
          j0 *= cj0;
        }
      }
      
      if (l == i && n == k)
      {
        if (!EQL(cj1,0))
        {
          for (jp = 1; jp < N[1]-1; ++jp)
          {
            dc_df = dc1_df(N[1],jp,m);
            j1 += dc_df*dT_dx((int)jp,y);
          }
          j1 *= cj1;
        }
      }
      
      if (l == i && m == j)
      {
        if (!EQL(cj2,0))
        {
          for (kp = 1; kp < N[2]-1; ++kp)
          {
            dc_df = dc2_df(N[2],kp,n);
            j2 += dc_df*dT_dx((int)kp,z);
          }
          j2 *= cj2;
        }
      }
      
      J[ijk][lmn] = j0+j1+j2;  
    }/* end of for (lmn = 0; lmn < nn; ++lmn) */
    
  }/* end of for (ijk = 0; ijk < nn; ++ijk) */
  
}

/* making Jacobian using spectral method in direction z
// d(df(i,j,k)/dz)/df(l,m,n) = j(N_i,z) * \sum_{ip,jp,kp} dc(ip,jp,kp)/df(l,m,n)*dT(i,j,k)/dN_i
*/
static void fill_jacobian_spectral_method_z(double **const J,Patch_T *const patch)
{
  const unsigned nn = total_nodes_patch(patch);
  const unsigned *const N = patch->n;
  unsigned ijk;
  
  OpenMP_1d_Pragma(omp parallel for)
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double cj0 = dq2_dq1(patch,_N0_,_z_,ijk);/* coordinate jacobian */
    double cj1 = dq2_dq1(patch,_N1_,_z_,ijk);/* coordinate jacobian */
    double cj2 = dq2_dq1(patch,_N2_,_z_,ijk);/* coordinate jacobian */
    double x,y,z;
    unsigned lmn;
    unsigned i,j,k;
    
    IJK(ijk,N,&i,&j,&k);
    x = ChebExtrema_1point(N[0],i);
    y = ChebExtrema_1point(N[1],j);
    z = ChebExtrema_1point(N[2],k);
    
    for (lmn = 0; lmn < nn; ++lmn)
    {
      double dc_df;
      double j0,j1,j2;
      unsigned l,m,n,ip,jp,kp;
      
      IJK(lmn,N,&l,&m,&n);
      j0 = 0;
      j1 = 0;
      j2 = 0;
      
      if (m == j && n == k)
      {
        if (!EQL(cj0,0))
        {
          for (ip = 1; ip < N[0]-1; ++ip)
          {
            dc_df = dc0_df(N[0],ip,l);
            j0 += dc_df*dT_dx((int)ip,x);
          }
          j0 *= cj0;
        }
      }
      
      if (l == i && n == k)
      {
        if (!EQL(cj1,0))
        {
          for (jp = 1; jp < N[1]-1; ++jp)
          {
            dc_df = dc1_df(N[1],jp,m);
            j1 += dc_df*dT_dx((int)jp,y);
          }
          j1 *= cj1;
        }
      }
      
      if (l == i && m == j)
      {
        if (!EQL(cj2,0))
        {
          for (kp = 1; kp < N[2]-1; ++kp)
          {
            dc_df = dc2_df(N[2],kp,n);
            j2 += dc_df*dT_dx((int)kp,z);
          }
          j2 *= cj2;
        }
      }
      
      J[ijk][lmn] = j0+j1+j2;  
    }/* end of for (lmn = 0; lmn < nn; ++lmn) */
    
  }/* end of for (ijk = 0; ijk < nn; ++ijk) */
  
}

/* making Jacobian using spectral method in direction xx
// d(d^2 f(i,j,k)/dx^2)/df(l,m,n) = d (d(df(i,j,k)/dx)/df(l,m,n))/dx
*/
static void fill_jacobian_spectral_method_xx(double **const J, Patch_T *const patch)
{
  const unsigned nn = total_nodes_patch(patch);
  const unsigned num_thread = (unsigned)omp_get_num_threads();
  Field_T **j_x;
  unsigned lmn,tn;
  char name[MAX_STR_LEN];
  
  j_x = malloc(num_thread*sizeof(*j_x));
  pointerEr(j_x);
  
  /* adding a field for each thread */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j_x%d",tn);
    j_x[tn] = add_field(name,"(3dim)",patch,YES);
  }
  
  fill_jacobian_spectral_method_x(J,patch);/* -> J_x */
  
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    unsigned ijk;
    double *J_xx = 0;
    Field_T *J_x = j_x[omp_get_thread_num()];
   
    for (ijk = 0; ijk < nn; ++ijk)
      J_x->v[ijk] = J[ijk][lmn];
      
    J_xx = Partial_Derivative(J_x,"x");
    /* since it was added v2 and info in J_x we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_info(J_x);
    free_v2(J_x);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_xx[ijk];/* -> J_xx */
    
    free(J_xx);
   
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j_x%d",tn);
    Field_T *f = patch->pool[Ind(name)];
    remove_field(f);
  }
  free(j_x);
}

/* making Jacobian using spectral method in direction yy
// d(d^2 f(i,j,k)/dy^2)/df(l,m,n) = d (d(df(i,j,k)/dy)/df(l,m,n))/dy
*/
static void fill_jacobian_spectral_method_yy(double **const J, Patch_T *const patch)
{
  const unsigned nn = total_nodes_patch(patch);
  const unsigned num_thread = (unsigned)omp_get_num_threads();
  Field_T **j_y;
  unsigned lmn,tn;
  char name[MAX_STR_LEN];
  
  j_y = malloc(num_thread*sizeof(*j_y));
  pointerEr(j_y);
  
  /* adding a field for each thread */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j_y%d",tn);
    j_y[tn] = add_field(name,"(3dim)",patch,YES);
  }
  
  fill_jacobian_spectral_method_y(J,patch);/* -> J_y */
  
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    unsigned ijk;
    double *J_yy = 0;
    Field_T *J_y = j_y[omp_get_thread_num()];
   
    for (ijk = 0; ijk < nn; ++ijk)
      J_y->v[ijk] = J[ijk][lmn];
      
    J_yy = Partial_Derivative(J_y,"y");
    /* since it was added v2 and info in J_y we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_info(J_y);
    free_v2(J_y);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_yy[ijk];/* -> J_yy */
    
    free(J_yy);
   
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j_y%d",tn);
    Field_T *f = patch->pool[Ind(name)];
    remove_field(f);
  }
  free(j_y);
}

/* making Jacobian using spectral method in direction yy
// d(d^2 f(i,j,k)/dy^2)/df(l,m,n) = d (d(df(i,j,k)/dy)/df(l,m,n))/dy
*/
static void fill_jacobian_spectral_method_zz(double **const J, Patch_T *const patch)
{
  const unsigned nn = total_nodes_patch(patch);
  const unsigned num_thread = (unsigned)omp_get_num_threads();
  Field_T **j_z;
  unsigned lmn,tn;
  char name[MAX_STR_LEN];
  
  j_z = malloc(num_thread*sizeof(*j_z));
  pointerEr(j_z);
  
  /* adding a field for each thread */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j_z%d",tn);
    j_z[tn] = add_field(name,"(3dim)",patch,YES);
  }
  
  fill_jacobian_spectral_method_z(J,patch);/* -> J_z */
  
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    unsigned ijk;
    double *J_zz = 0;
    Field_T *J_z = j_z[omp_get_thread_num()];
   
    for (ijk = 0; ijk < nn; ++ijk)
      J_z->v[ijk] = J[ijk][lmn];
      
    J_zz = Partial_Derivative(J_z,"z");
    /* since it was added v2 and info in J_y we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_info(J_z);
    free_v2(J_z);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_zz[ijk];/* -> J_zz */
    
    free(J_zz);
   
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j_z%d",tn);
    Field_T *f = patch->pool[Ind(name)];
    remove_field(f);
  }
  free(j_z);
}

/* dc/df where c is coefficients of expansion in direction 0
// ->return value: dc(i)/df(l)
*/
static double dc0_df(const unsigned n0,const unsigned i,const unsigned l)
{
  double dc_df = 0;
  
  if (l == 0)
    dc_df = 1;
  else if (l == n0-1)
    dc_df = SIGN[i%2];
  else
  {
    double xi = ChebExtrema_1point(n0,i);
    dc_df = 2 * Cheb_Tn((int)l,xi);
  }
  
  return dc_df/(2*(n0-1));
}

/* dc/df where c is coefficients of expansion in direction 1
// ->return value: dc(j)/df(m)
*/
static double dc1_df(const unsigned n1,const unsigned j,const unsigned m)
{
  double dc_df = 0;
  
  if (m == 0)
    dc_df = 1;
  else if (m == n1-1)
    dc_df = SIGN[j%2];
  else
  {
    double yj = ChebExtrema_1point(n1,j);
    dc_df = 2 * Cheb_Tn((int)m,yj);
  }
  
  return dc_df/(2*(n1-1));
}

/* dc/df where c is coefficients of expansion in direction 2
// ->return value: dc(k)/df(n)
*/
static double dc2_df(const unsigned n2,const unsigned k,const unsigned n)
{
  double dc_df = 0;
  
  if (n == 0)
    dc_df = 1;
  else if (n == n2-1)
    dc_df = SIGN[k%2];
  else
  {
    double zk = ChebExtrema_1point(n2,k);
    dc_df = 2 * Cheb_Tn((int)n,zk);
  }
  
  return dc_df/(2*(n2-1));
}

/* given number of point in a line and the point number
// it gives the Chebyshev extrema.
// ->return value: Chebyshev extrema at p
*/
static double ChebExtrema_1point(const unsigned n, const unsigned p)
{
  return cos(p*M_PI/(n-1));
}

/* ->return value: d(Cheb_Tn(x))/dx */
static double dT_dx(const int n,const double x)
{
  return n*Cheb_Un(n-1,x);
}
