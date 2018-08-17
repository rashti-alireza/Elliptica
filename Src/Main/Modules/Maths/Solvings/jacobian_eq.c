/*
// Alireza Rashti
// August 2018
*/

#include "jacobian_eq.h"
#define MAX_STR_LEN 400
static const double CONST = 1.0;

/* making Jacobian for equations at the inner mesh
// types are pointers to string determining the type of jacobian
// e.g. *types[3] = {"J_xx","J_y",0}.Note: the number of
// types is found by null.
*/
void make_jacobian_eq(Grid_T *const grid, const char * const* types)
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
    jt_e = str2JType_E(types[i]);
    
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      nn = total_nodes_patch(patch);
      J = alloc_matrix(nn,nn);
      Jacobian(J,patch,jt_e);
      
      free_matrix(J,nn);
    }
    i++;
  }
}

/* testing make_jacobian_eq function:
// it checks the consistensy between direct and spectral methods 
*/
void test_make_jacobian_eq(Grid_T *const grid, const char * const* types)
{
  enum Method_E {Spectral_e = 0,Direct_e,N_Method_E};
  Jacobian_eq_F *Jacobian[N_Method_E] = 
      {make_jacobian_spectral_method,make_jacobian_direct_method};
  double **cmp[N_Method_E];
  double **J = 0;
  const char *path_par = GetParameterS_E("output_directory_path");
  char *path = make_directory(path_par,"Test_Jacobian_Eq");
  char file_name[MAX_STR_LEN];
  char line[MAX_STR_LEN]={'\0'};
  FILE *file = 0;
  double Err = 0;
  JType_E jt_e;
  unsigned i,p,nn,r,c;
  enum Method_E e;
  Flag_T flg = NONE;
  
  i = 0;
  while (types[i] != 0)
  {
    jt_e = str2JType_E(types[i]);
    
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      nn = total_nodes_patch(patch);
      Err = CONST/nn;
      
      for (e = Spectral_e; e < N_Method_E; ++e)
      {
        J = alloc_matrix(nn,nn);
        Jacobian[e](J,patch,jt_e);
        cmp[e] = J;
      }
      
      printf("Testing Jacobian for Equations: patch=%s, type:%5s\t",patch->name,types[i]);
      
      sprintf(file_name,"%s/%s_SepctalDirect.patch%u",path,types[i],patch->pn);
      file = fopen(file_name,"w");
      pointerEr(file);
      fprintf(file,"Row Column J_Spectal J_Direct\n");
      
      for (r = 0; r < nn; ++r)
      {
        for (c = 0; c < nn; ++c)
        {
          if (GRT(ABS(cmp[Spectral_e][r][c]-cmp[Direct_e][r][c]),Err))
            fprintf(file,"%3u %3u   %0.10f    %0.10f\n",
              r,c,cmp[Spectral_e][r][c],cmp[Direct_e][r][c]);
        }
      }
      fclose(file);
      
      for (e = Spectral_e; e < N_Method_E; ++e)
      {
        free_matrix(cmp[e],nn);
      }
      
      flg = NO;
      /* check if the second line is empty so both approach are equal */
      file = fopen(file_name,"r");
      fgets(line,sizeof(line),file);
      if(fgets(line,sizeof(line),file) == 0)
        flg = YES;
        
      if (flg == YES) printf("[+].\n");
      else	      printf("[-].\n");

    }
    i++;
  }
  
  free(path);
}

/* translating string to enum JType_E */
static JType_E str2JType_E(const char *const str)
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
static void make_jacobian_direct_method(double **const J,Patch_T *const patch,const JType_E jt_e)
{
  switch(jt_e)
  {
    case T_x:
      fill_jacobian_direct_method_1stOrder(J,patch,T_x);
      break;
    case T_xx:
      fill_jacobian_direct_method_2ndOrder(J,patch,T_xx);
      break;
    case T_y:
      fill_jacobian_direct_method_1stOrder(J,patch,T_y);
      break;
    case T_yy:
      fill_jacobian_direct_method_2ndOrder(J,patch,T_yy);
      break;
    case T_z:
      fill_jacobian_direct_method_1stOrder(J,patch,T_z);
      break;
    case T_zz:
      fill_jacobian_direct_method_2ndOrder(J,patch,T_zz);
      break;
    default:
      abortEr("No such type for Jacobian defined!\n");
  }
}

/* making Jacobian using direct method in direction $
// d(df(i,j,k)/d$)/df(l,m,n) = (d(f+df)/d$-df/d$)/df
*/
static void fill_jacobian_direct_method_1stOrder(double **const J, Patch_T *const patch,const JType_E jt_e)
{
  Field_T **j_1st_deriv_field = 0;
  Patch_T *temp_patch = 0;
  const unsigned nn = total_nodes_patch(patch);
  const double EPS = CONST/nn;
  unsigned num_thread;
  char name[MAX_STR_LEN],deriv_str[MAX_STR_LEN] ;
  unsigned lmn,tn;
  
  JType_E2str(jt_e,deriv_str);
  
  #pragma omp parallel 
  { 
    #pragma omp single
    {
      num_thread = (unsigned)omp_get_num_threads();
      j_1st_deriv_field = malloc(num_thread*sizeof(*j_1st_deriv_field));
      pointerEr(j_1st_deriv_field);
      temp_patch = malloc(num_thread*sizeof(*temp_patch));
      pointerEr(temp_patch);
  
      /* adding a patch and field for each thread */
      for (tn = 0; tn < num_thread; ++tn)
      {
        sprintf(name,"j_1st_deriv_field%d",tn);
        temp_patch[tn] = make_temp_patch(patch);
        j_1st_deriv_field[tn] = add_field(name,"(3dim)",&temp_patch[tn],YES);
        
        for (lmn = 0; lmn < nn; ++lmn)
          j_1st_deriv_field[tn]->v[lmn] = CONST;
      }
    }
  }
    
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Field_T *Jf = j_1st_deriv_field[omp_get_thread_num()];
    double *J_deriv = 0;
    unsigned ijk;
    
    Jf->v[lmn] += EPS;   
    J_deriv = Partial_Derivative(Jf,deriv_str);
    /* since it was added v2 and info in J_deriv we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_coeffs(Jf);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_deriv[ijk]/EPS;
      
    free(J_deriv);
    Jf->v[lmn] -= EPS;
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j_1st_deriv_field%d",tn);
    Field_T *f = temp_patch[tn].pool[LookUpField(name,&temp_patch[tn])];
    remove_field(f);
    free_temp_patch(&temp_patch[tn]);
  }
  free(j_1st_deriv_field);
}

/* making Jacobian using direct method in direction $&
// d(d^2f(i,j,k)/d$d&)/df(l,m,n) = (d^2(f+df)/d$d&-d^2f/d$d&)/df
*/
static void fill_jacobian_direct_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir)
{
  const unsigned nn = total_nodes_patch(patch);
  const double EPS = CONST/nn;
  unsigned num_thread;
  Field_T **j;
  Patch_T *temp_patch = 0;
  unsigned lmn,tn;
  char name[MAX_STR_LEN];
  char deriv_str[MAX_STR_LEN];
  
  JType_E2str(deriv_dir,deriv_str);

  #pragma omp parallel 
  { 
    #pragma omp single
    {
      num_thread = (unsigned)omp_get_num_threads();
      j = malloc(num_thread*sizeof(*j));
      pointerEr(j);
      temp_patch = malloc(num_thread*sizeof(*temp_patch));
      pointerEr(temp_patch);
      
      /* adding a field for each thread */
      for (tn = 0; tn < num_thread; ++tn)
      {
        sprintf(name,"j%d",tn);
        temp_patch[tn] = make_temp_patch(patch);
        j[tn] = add_field(name,"(3dim)",&temp_patch[tn],YES);
        
        for (lmn = 0; lmn < nn; ++lmn)
          j[tn]->v[lmn] = CONST;
      }
    }
  }
  
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Field_T *Jf = j[omp_get_thread_num()];
    double *J_2nd = 0;
    unsigned ijk;
    
    Jf->v[lmn] += EPS;   
    J_2nd = Partial_Derivative(Jf,deriv_str);
    /* since it was added v2 and info in J_2nd we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_coeffs(Jf);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_2nd[ijk]/EPS;
      
    free(J_2nd);
    Jf->v[lmn] -= EPS;
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j%d",tn);
    Field_T *f = temp_patch[tn].pool[LookUpField(name,&temp_patch[tn])];
    remove_field(f);
    free_temp_patch(&temp_patch[tn]);
  }
  free(j);
}

/* making Jacobian equations using spectral method */
static void make_jacobian_spectral_method(double **const J,Patch_T *const patch,const JType_E jt_e)
{
  switch(jt_e)
  {
    case T_x:
      fill_jacobian_spectral_method_1stOrder(J,patch,T_x);
      break;
    case T_xx:
      fill_jacobian_spectral_method_2ndOrder(J,patch,T_xx);
      break;
    case T_y:
      fill_jacobian_spectral_method_1stOrder(J,patch,T_y);
      break;
    case T_yy:
      fill_jacobian_spectral_method_2ndOrder(J,patch,T_yy);
      break;
    case T_z:
      fill_jacobian_spectral_method_1stOrder(J,patch,T_z);
      break;
    case T_zz:
      fill_jacobian_spectral_method_2ndOrder(J,patch,T_zz);
      break;
    default:
      abortEr("No such type for Jacobian defined!\n");
  }
}

/* making Jacobian using spectral method in direction $
// d(df(i,j,k)/d$)/df(l,m,n) = j(N_i,$) * 2 \sum_{ip,jp,kp} dc(ip,jp,kp)/df(l,m,n)*dT(i,j,k)/dN_i
*/
static void fill_jacobian_spectral_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e)
{
  const unsigned nn = total_nodes_patch(patch);
  const unsigned *const N = patch->n;
  Dd_T q_dir;
  unsigned ijk;
  
  JType_E2Dd_T(jt_e,&q_dir);
  
  OpenMP_1d_Pragma(omp parallel for)
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double cj0 = dq2_dq1(patch,_N0_,q_dir,ijk);/* coordinate jacobian */
    double cj1 = dq2_dq1(patch,_N1_,q_dir,ijk);/* coordinate jacobian */
    double cj2 = dq2_dq1(patch,_N2_,q_dir,ijk);/* coordinate jacobian */
    double x,y,z;
    unsigned lmn;
    unsigned i,j,k;
    
    IJK(ijk,N,&i,&j,&k);
    x = ChebExtrema_1point(N[0],i);
    y = ChebExtrema_1point(N[1],j);
    z = ChebExtrema_1point(N[2],k);
    
    for (lmn = 0; lmn < nn; ++lmn)
    {
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
            j0 += dc_df(N[0],ip,l)*dT_dx((int)ip,x);
          }
          j0 *= 2*cj0;
        }
      }
      
      if (l == i && n == k)
      {
        if (!EQL(cj1,0))
        {
          for (jp = 1; jp < N[1]-1; ++jp)
          {
            j1 += dc_df(N[1],jp,m)*dT_dx((int)jp,y);
          }
          j1 *= 2*cj1;
        }
      }
      
      if (l == i && m == j)
      {
        if (!EQL(cj2,0))
        {
          for (kp = 1; kp < N[2]-1; ++kp)
          {
            j2 += dc_df(N[2],kp,n)*dT_dx((int)kp,z);
          }
          j2 *= 2*cj2;
        }
      }
      
      J[ijk][lmn] = j0+j1+j2;  
    }/* end of for (lmn = 0; lmn < nn; ++lmn) */
    
  }/* end of for (ijk = 0; ijk < nn; ++ijk) */
  
}

/* reading 1st order derivative for given JType_E */
static void JType_E2Dd_T(const JType_E jt_e, Dd_T *const q_dir)
{
  switch(jt_e)
  {
    case T_x:
      *q_dir = _x_;
      break;
    case T_y:
      *q_dir = _y_;
      break;
    case T_z:
      *q_dir = _z_;
      break;
    default:
      abortEr(INCOMPLETE_FUNC);
  }
}

/* making Jacobian using spectral method in direction @&
// d(d^2 f(i,j,k)/d@d&)/df(l,m,n) = d (d(df(i,j,k)/d@)/df(l,m,n))/d&
*/
static void fill_jacobian_spectral_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir)
{
  const unsigned nn = total_nodes_patch(patch);
  unsigned num_thread;
  Field_T **j_1st_deriv_field = 0;
  Patch_T *temp_patch = 0;
  JType_E deriv_1st,deriv_2nd;
  char name[MAX_STR_LEN],deriv_2nd_s[MAX_STR_LEN];
  unsigned lmn,tn;
  
  read_1st_and_2nd_deriv(deriv_dir,&deriv_1st,&deriv_2nd);
  JType_E2str(deriv_2nd,deriv_2nd_s);
  
  #pragma omp parallel 
  { 
    #pragma omp single
    {
      num_thread = (unsigned)omp_get_num_threads();
      j_1st_deriv_field = malloc(num_thread*sizeof(*j_1st_deriv_field));
      pointerEr(j_1st_deriv_field);
      temp_patch = malloc(num_thread*sizeof(*temp_patch));
      pointerEr(temp_patch);
      
      /* adding a patch and field for each thread to make them thread safe*/
      for (tn = 0; tn < num_thread; ++tn)
      {
        sprintf(name,"j_1st_deriv_field%d",tn);
        temp_patch[tn] = make_temp_patch(patch);
        j_1st_deriv_field[tn] = add_field(name,"(3dim)",&temp_patch[tn],YES);
      }
    }
  }
  
  fill_jacobian_spectral_method_1stOrder(J,patch,deriv_1st);/* -> J = d(df/d@)/df */
  
  OpenMP_1d_Pragma(omp parallel for)
  for (lmn = 0; lmn < nn; ++lmn)
  {
    unsigned ijk;
    double *J_2nd_deriv_value = 0;
    Field_T *J_1st_deriv_field = j_1st_deriv_field[omp_get_thread_num()];
   
    for (ijk = 0; ijk < nn; ++ijk)
      J_1st_deriv_field->v[ijk] = J[ijk][lmn];
      
    J_2nd_deriv_value = Partial_Derivative(J_1st_deriv_field,deriv_2nd_s);
    /* since it was added v2 and info in J_1st_deriv_field we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_coeffs(J_1st_deriv_field);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = J_2nd_deriv_value[ijk];/* -> J_2nd_deriv_value */
    
    free(J_2nd_deriv_value);
   
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  for (tn = 0; tn < num_thread; ++tn)
  {
    sprintf(name,"j_1st_deriv_field%d",tn);
    Field_T *f = temp_patch[tn].pool[LookUpField(name,&temp_patch[tn])];
    remove_field(f);
    free_temp_patch(&temp_patch[tn]);
  }
  free(j_1st_deriv_field);
}

/* decomposing deriv_dir and finding 1st and 2nd direction for derivative */
static void read_1st_and_2nd_deriv(const JType_E deriv_dir,JType_E *const deriv_1st,JType_E *const deriv_2nd)
{
  switch(deriv_dir)
  {
    case T_xx:
      *deriv_1st = T_x;
      *deriv_2nd = T_x;
      break;
    case T_yy:
      *deriv_1st = T_y;
      *deriv_2nd = T_y;
      break;
    case T_zz:
      *deriv_1st = T_z;
      *deriv_2nd = T_z;
      break;
    default:
      abortEr(INCOMPLETE_FUNC);
  }
}

/* dc/df where c is coefficients of expansion in a direction with n nodes
// ->return value: dc(i)/df(l)
*/
static double dc_df(const unsigned n,const unsigned i,const unsigned l)
{
  double dcdf = 0;
  
  if (l == 0)
    dcdf = 1;
  else if (l == n-1)
    dcdf = SIGN[i%2];
  else
  {
    double xi = ChebExtrema_1point(n,i);
    dcdf = 2 * Cheb_Tn((int)l,xi);
  }
  
  return dcdf/(2*(n-1));
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

/* converitng enum to str */
static void JType_E2str(const JType_E e,char *const str)
{
  switch(e)
  {
    case T_x:
      sprintf(str,"x");
      break;
    case T_y:
      sprintf(str,"y");
      break;
    case T_z:
      sprintf(str,"z");
      break;
    case T_xx:
      sprintf(str,"x,x");
      break;
    case T_yy:
      sprintf(str,"y,y");
      break;
    case T_zz:
      sprintf(str,"z,z");
      break;
    default:
      abortEr(INCOMPLETE_FUNC);
  }
}
