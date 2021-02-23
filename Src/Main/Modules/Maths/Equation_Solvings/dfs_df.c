/*
// Alireza Rashti
// August 2018
*/

#include "dfs_df.h"
#define MAX_STR_LEN (400)
#define MAX_J_SIZE (10000) /* 10Gb */

static const double CONST = 1.0;

/* preparing J_* used in equations for each patch
// according to given types. note: end of types pointer
// must be marked with null pointer, e.g. *types[3] = {"dfxx_df","dfy_df",0}.
*/
void prepare_Js_jacobian_eq(Patch_T *const patch,const char * const *types)
{
  Js_Jacobian_eq_F *Jacobian = 0;
  Solving_Man_T *const sol_man = patch->solving_man;
  const Uint nn = patch->nn;
  /* default value of J if it gets larger than 10 Mb it will 
  // be written in a file unless user assigne other value. */
  double max_j_size = MAX_J_SIZE;
  Matrix_T *J = 0;
  char *jtype = 0;
  JType_E jt_e = T_UNDEF;
  Uint i;
  
  if (get_parameter("Maximum_Size_of_J_Kept_in_Mb"))
    max_j_size = PgetdEZ("Maximum_Size_of_J_Kept_in_Mb");
  
  /* selecting Jacobian method for making of jacobian equation */
  if (strcmp_i(Pgets("dF/du_for_Newton_Method"),"Spectral"))
    Jacobian = make_jacobian_spectral_method;
  else if (strcmp_i(Pgets("dF/du_for_Newton_Method"),"Finite_Difference"))
    Jacobian = make_jacobian_direct_method;
  else
    Error0(INCOMPLETE_FUNC);
  
  i = 0;
  while (types[i] != 0)
  {
    jtype = interpret_type(types[i]);
    
    /* check if this type has been already made then skip this */
    Flag_T flg = NONE;
    Uint c;
    for (c = 0; c < sol_man->nj; ++c)
      if (strcmp_i(sol_man->jacobian[c]->type,jtype))
      {
        flg = FOUND;
        break;
      }
    
    if (flg == FOUND)
    {
      i++;
      free(jtype);
      continue;
    }
      
    /* making Jacobian elements */
    jt_e = str2JType_E(jtype);
    J = alloc_matrix(REG_SF,nn,nn);
    Jacobian(J->reg->A,patch,jt_e);
    
    /* saving jacobian elements */
    c = sol_man->nj;
    sol_man->jacobian = 
      realloc(sol_man->jacobian,(c+1)*sizeof(*sol_man->jacobian));
    IsNull(sol_man->jacobian);
    sol_man->jacobian[c] = calloc(1,sizeof(*sol_man->jacobian[c]));
    IsNull(sol_man->jacobian[c]);
    sprintf(sol_man->jacobian[c]->type,jtype);
    sol_man->jacobian[c]->J = cast_matrix_ccs(J);
    free_matrix(J);
    
    if (GRT(J_sizeMb_ccs(sol_man->jacobian[i]->J),max_j_size))
      write_J_in_disk_ccs();
    
    free(jtype);
    
    ++i;
    ++sol_man->nj;
    
  }
}

/* this function interprets what the given type means for argument of jacobian equations.
// this function was developed to deal with cases that we wanna treat jacobians as the components
// of a tensor to expand them.
// ->return value: the type which is undestandable by make_jacobian_* functions. */
static char *interpret_type(const char *const type)
{
 char *jtype = calloc(MAX_STR_LEN,1);
 Uint i,len;
 
 IsNull(jtype);
 
 /*  possible given types are *_D?D?... or _U?U?... or U?D?... or D?U? ... 
 or simply df*_df */
 
 if (regex_search("^df[xyz]+_df$",type))
  sprintf(jtype,"%s",type);
 else if (regex_search("_([DU][[:digit:]])+$",type))
 {
  char *match = regex_find("_([DU][[:digit:]])+$",type);/* e.g _U2D0 */
  sprintf(jtype,"%s","df");
  len = (Uint)strlen(match);
  for (i = 2; i < len; ++i)/* starting from number */
  {
   if (match[i] == '0')      strcat(jtype,"x");
   else if (match[i] == '1') strcat(jtype,"y");
   else if (match[i] == '2') strcat(jtype,"z");
   else if (match[i] == 'U' || match[i] == 'D') continue;
   else Errors("This type '%s' has not been defined.\n",type);
  }
  
  strcat(jtype,"_df");
  
  free(match);
 }
 else
  Errors("The type '%s' cannot be realized!\n",type);
 
  return jtype;
}

/* making elements of Jacobian for equations at the inner mesh.
// types are pointers to string determining the type of jacobian
// e.g. *types[3] = {"dfxx_df","dfy_df",0}.Note: the number of
// types is found by null.
*/
void make_Js_jacobian_eq(Grid_T *const grid, const char * const* types)
{
  Js_Jacobian_eq_F *Jacobian = 0;
  Matrix_T *J = 0;
  JType_E jt_e = T_UNDEF;
  char *jtype = 0;
  Uint i,p,nn;
  
  /* selecting Jacobian method for making of jacobian equation */
  if (strcmp_i(Pgets("Making_Jacobian_Eq_Method"),"spectral"))
    Jacobian = make_jacobian_spectral_method;
  else if (strcmp_i(Pgets("Making_Jacobian_Eq_Method"),"direct"))
    Jacobian = make_jacobian_direct_method;
  else
    Error0(INCOMPLETE_FUNC);
  
  i = 0;
  while (types[i] != 0)
  {
    jtype = interpret_type(types[i]);
    jt_e = str2JType_E(jtype);
    
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      nn = total_nodes_patch(patch);
      J = alloc_matrix(REG_SF,nn,nn);
      Jacobian(J->reg->A,patch,jt_e);
      printf("This function is not ready yet!\n%s,%d\n",__FILE__,__LINE__);
      free_matrix(J);
    }
    free(jtype);
    i++;
  }
}

/* testing make_jacobian_eq function:
// it checks the consistensy between direct and spectral methods 
*/
void test_make_Js_jacobian_eq(Grid_T *const grid, const char * const* types)
{
  enum Method_E {Spectral_e = 0,Direct_e,N_Method_E};
  Js_Jacobian_eq_F *Jacobian[N_Method_E] = 
      {make_jacobian_spectral_method,make_jacobian_direct_method};
  double **cmp[N_Method_E];
  Matrix_T *J = 0;
  const char *path_par = Pgets("top_directory");
  char *path = make_directory(path_par,"Test_Jacobian_Eq");
  char file_name[MAX_STR_LEN];
  char line[MAX_STR_LEN]={'\0'};
  FILE *file = 0;
  double Err = 0;
  char *jtype = 0;
  JType_E jt_e;
  Uint i,p,nn,r,c;
  enum Method_E e;
  Flag_T flg = NONE;
  
  i = 0;
  while (types[i] != 0)
  {
    jtype = interpret_type(types[i]);
    jt_e = str2JType_E(jtype);
    
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      nn = total_nodes_patch(patch);
      Err = CONST/nn;
      
      for (e = Spectral_e; e < N_Method_E; ++e)
      {
        J = alloc_matrix(REG_SF,nn,nn);
        Jacobian[e](J->reg->A,patch,jt_e);
        cmp[e] = J->reg->A;
        J->reg->A = 0;
        free_matrix(J);
      }
      
      printf("Testing Jacobian for Equations: patch=%s, type:%5s\t",patch->name,jtype);
      
      sprintf(file_name,"%s/%s_SepctalDirect.patch%u",path,jtype,patch->pn);
      file = Fopen(file_name,"w");
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
      Fclose(file);
      
      for (e = Spectral_e; e < N_Method_E; ++e)
      {
        free_2d_mem(cmp[e],nn);
      }
      
      flg = NO;
      /* check if the second line is empty so both approach are equal */
      file = Fopen(file_name,"r");
      fgets(line,sizeof(line),file);
      if(fgets(line,sizeof(line),file) == 0)
        flg = YES;
        
      if (flg == YES) printf("[+].\n");
      else	      printf("[-].\n");

    }
    i++;
    
    free(jtype);
  }
  
  free(path);
}

/* translating string to enum JType_E */
static JType_E str2JType_E(const char *const str)
{
  JType_E jt_e = T_UNDEF;
  
  if (strcmp_i(str,"dfx_df"))
    jt_e = T_x;
    
  else if (strcmp_i(str,"dfxx_df"))
    jt_e = T_xx;
    
  else if (strcmp_i(str,"dfy_df"))
    jt_e = T_y;
    
  else if (strcmp_i(str,"dfyy_df"))
    jt_e = T_yy;
    
  else if (strcmp_i(str,"dfz_df"))
    jt_e = T_z;
    
  else if (strcmp_i(str,"dfzz_df"))
    jt_e = T_zz;
    
  else if (strcmp_i(str,"dfxy_df") || strcmp_i(str,"dfyx_df"))
    jt_e = T_xy;

  else if (strcmp_i(str,"dfxz_df") || strcmp_i(str,"dfzx_df"))
    jt_e = T_xz;
    
  else if (strcmp_i(str,"dfyz_df") || strcmp_i(str,"dfzy_df"))
    jt_e = T_yz;
    
  else
    Error0(INCOMPLETE_FUNC);
  
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
    case T_xy:
      fill_jacobian_direct_method_2ndOrder(J,patch,T_xy);
      break;
    case T_xz:
      fill_jacobian_direct_method_2ndOrder(J,patch,T_xz);
      break;
      
    case T_y:
      fill_jacobian_direct_method_1stOrder(J,patch,T_y);
      break;
    case T_yy:
      fill_jacobian_direct_method_2ndOrder(J,patch,T_yy);
      break;
    case T_yz:
      fill_jacobian_direct_method_2ndOrder(J,patch,T_yz);
      break;
      
    case T_z:
      fill_jacobian_direct_method_1stOrder(J,patch,T_z);
      break;
    case T_zz:
      fill_jacobian_direct_method_2ndOrder(J,patch,T_zz);
      break;
    default:
      Error0("No such type for Jacobian defined!\n");
  }
}

/* making Jacobian using direct method in direction $
// d(df(i,j,k)/d$)/df(l,m,n) = (d(f+df)/d$-df/d$)/df
*/
static void fill_jacobian_direct_method_1stOrder(double **const J, Patch_T *const patch,const JType_E jt_e)
{
  Field_T *j_1st_deriv_field = 0;
  Patch_T temp_patch;
  const Uint nn = patch->nn;
  const double EPS = CONST/nn;
  char deriv_str[MAX_STR_LEN] ;
  Uint lmn;
  
  JType_E2str(jt_e,deriv_str);
  
  temp_patch = make_temp_patch(patch);
  j_1st_deriv_field = add_field("j_1st_deriv_field","(3dim)",&temp_patch,YES);
  
  for (lmn = 0; lmn < nn; ++lmn)
    j_1st_deriv_field->v[lmn] = CONST;
    
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Field_T *Jf = j_1st_deriv_field;
    double *J_deriv = 0;
    Uint ijk;
    
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
  
  Field_T *f = temp_patch.fields[LookUpField("j_1st_deriv_field",&temp_patch)];
  remove_field(f);
  free_temp_patch(&temp_patch);
}

/* making Jacobian using direct method in direction $&
// d(d^2f(i,j,k)/d$d&)/df(l,m,n) = (d^2(f+df)/d$d&-d^2f/d$d&)/df
*/
static void fill_jacobian_direct_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir)
{
  const Uint nn = patch->nn;
  const double EPS = CONST/nn;
  Field_T *j;
  Patch_T temp_patch;
  Uint lmn;
  char deriv_str[MAX_STR_LEN];
  
  JType_E2str(deriv_dir,deriv_str);

  temp_patch = make_temp_patch(patch);
  j = add_field("j","(3dim)",&temp_patch,YES);
  
  for (lmn = 0; lmn < nn; ++lmn)
    j->v[lmn] = CONST;
        
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Field_T *Jf = j;
    double *J_2nd = 0;
    Uint ijk;
    
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
  Field_T *f = temp_patch.fields[LookUpField("j",&temp_patch)];
  remove_field(f);
  free_temp_patch(&temp_patch);
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
    case T_xy:
      fill_jacobian_spectral_method_2ndOrder(J,patch,T_xy);
      break;
    case T_xz:
      fill_jacobian_spectral_method_2ndOrder(J,patch,T_xz);
      break;
      
    case T_y:
      fill_jacobian_spectral_method_1stOrder(J,patch,T_y);
      break;
    case T_yy:
      fill_jacobian_spectral_method_2ndOrder(J,patch,T_yy);
      break;
    case T_yz:
      fill_jacobian_spectral_method_2ndOrder(J,patch,T_yz);
      break;
      
    case T_z:
      fill_jacobian_spectral_method_1stOrder(J,patch,T_z);
      break;
    case T_zz:
      fill_jacobian_spectral_method_2ndOrder(J,patch,T_zz);
      break;
      
    default:
      Error0("No such type for Jacobian defined!\n");
  }
}

/* making Jacobian using spectral method in direction $
// d(df(i)/d$)/df(l) = j(N_i,$) *(2 \sum_{ip=1}^{n-2} dc(ip)/df(l)*dT(ip)/dN_i + dc(n-1)/df(l)*dT(n-1)/dN_i)
*/
void obsolete_fill_jacobian_spectral_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e)
{
  const Uint nn = patch->nn;
  const Uint *const N = patch->n;
  Dd_T q_dir = UNDEFINED_DIR;
  Uint ijk;
  
  JType_E2Dd_T(jt_e,&q_dir);
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double cj0 = dq2_dq1(patch,_N0_,q_dir,ijk);/* coordinate jacobian */
    double cj1 = dq2_dq1(patch,_N1_,q_dir,ijk);/* coordinate jacobian */
    double cj2 = dq2_dq1(patch,_N2_,q_dir,ijk);/* coordinate jacobian */
    double x,y,z;
    Uint lmn;
    Uint i,j,k;
    
    ijk_to_i_j_k(ijk,N,&i,&j,&k);
    x = ChebExtrema_1point(N[0],i);
    y = ChebExtrema_1point(N[1],j);
    z = ChebExtrema_1point(N[2],k);
    
    for (lmn = 0; lmn < nn; ++lmn)
    {
      double j0,j1,j2;
      Uint l,m,n,ip,jp,kp;
      
      ijk_to_i_j_k(lmn,N,&l,&m,&n);
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
          j0 *= 2;
          j0 += dc_df(N[0],ip,l)*dT_dx((int)ip,x);
          j0 *= cj0;
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
          j1 *= 2;
          j1 += dc_df(N[1],jp,m)*dT_dx((int)jp,y);
          j1 *= cj1;
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
          j2 *= 2;
          j2 += dc_df(N[2],kp,n)*dT_dx((int)kp,z);
          j2 *= cj2;
        }
      }
      
      J[ijk][lmn] = j0+j1+j2;  
    }/* end of for (lmn = 0; lmn < nn; ++lmn) */
    
  }/* end of for (ijk = 0; ijk < nn; ++ijk) */
  
}

/* making Jacobian using spectral method in direction $
// d(df(i)/d$)/df(l) = j(N_i,$) *(2 \sum_{ip=1}^{n-2} dc(ip)/df(l)*dT(ip)/dN_i + dc(n-1)/df(l)*dT(n-1)/dN_i)
*/
static void fill_jacobian_spectral_method_1stOrder(double **const J,Patch_T *const patch,const JType_E jt_e)
{
  const Uint nn = patch->nn;
  const Uint *const N = patch->n;
  Dd_T q_dir = UNDEFINED_DIR;
  Uint ijk;
  
  JType_E2Dd_T(jt_e,&q_dir);
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double dN0_dq = dq2_dq1(patch,_N0_,q_dir,ijk);/* coordinate jacobian */
    double dN1_dq = dq2_dq1(patch,_N1_,q_dir,ijk);/* coordinate jacobian */
    double dN2_dq = dq2_dq1(patch,_N2_,q_dir,ijk);/* coordinate jacobian */
    double x,y,z;
    Uint i,j,k,l,m,n;
    ijk_to_i_j_k(ijk,N,&i,&j,&k);
    
    x = ChebExtrema_1point(N[0],i);
    y = ChebExtrema_1point(N[1],j);
    z = ChebExtrema_1point(N[2],k);
    
    if (!EQL(dN0_dq,0))
    {
      for (l = 0; l < N[0]; ++l)
        J[ijk][i_j_k_to_ijk(N,l,j,k)] += dN0_dq*sum_0_N_dCi_dfj_by_dTi_dq(N[0],l,x);
    }
    if (!EQL(dN1_dq,0))
    {
      for (m = 0; m < N[1]; ++m)
        J[ijk][i_j_k_to_ijk(N,i,m,k)] += dN1_dq*sum_0_N_dCi_dfj_by_dTi_dq(N[1],m,y);
    }
    if (!EQL(dN2_dq,0))
    {
      for (n = 0; n < N[2]; ++n)
        J[ijk][i_j_k_to_ijk(N,i,j,n)] += dN2_dq*sum_0_N_dCi_dfj_by_dTi_dq(N[2],n,z);
    }
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
      Error0(INCOMPLETE_FUNC);
  }
}

/* making Jacobian using spectral method in direction @&
// d(d^2 f(i,j,k)/d@d&)/df(l,m,n) = d (d(df(i,j,k)/d@)/df(l,m,n))/d&
*/
static void fill_jacobian_spectral_method_2ndOrder(double **const J, Patch_T *const patch,const JType_E deriv_dir)
{
  const Uint nn = patch->nn;
  Field_T *j_1st_deriv_field = 0;
  Patch_T temp_patch;
  JType_E deriv_1st = T_UNDEF,deriv_2nd = T_UNDEF;
  char deriv_2nd_s[MAX_STR_LEN];
  Uint lmn;
  
  read_1st_and_2nd_deriv(deriv_dir,&deriv_1st,&deriv_2nd);
  JType_E2str(deriv_2nd,deriv_2nd_s);
  
  temp_patch = make_temp_patch(patch);
  j_1st_deriv_field = add_field("j_1st_deriv_field","(3dim)",&temp_patch,YES);
  
  fill_jacobian_spectral_method_1stOrder(J,patch,deriv_1st);/* -> J = d(df/d@)/df */
  
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Uint ijk;
    double *j_2nd_deriv_value = 0;
   
    for (ijk = 0; ijk < nn; ++ijk)
      j_1st_deriv_field->v[ijk] = J[ijk][lmn];
      
    j_2nd_deriv_value = Partial_Derivative(j_1st_deriv_field,deriv_2nd_s);
    /* since it was added v2 and info in j_1st_deriv_field we clean them 
    // to avoid using them again for new data in next iteration
    */
    free_coeffs(j_1st_deriv_field);
    
    for (ijk = 0; ijk < nn; ++ijk)
      J[ijk][lmn] = j_2nd_deriv_value[ijk];/* -> j_2nd_deriv_value */
    
    free(j_2nd_deriv_value);
   
  }/* end of for (lmn = 0; lmn < nn; ++lmn) */
  
  /* removing field and freeing memories */
  Field_T *f = temp_patch.fields[LookUpField("j_1st_deriv_field",&temp_patch)];
  remove_field(f);
  free_temp_patch(&temp_patch);
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
    case T_xy:
      *deriv_1st = T_x;
      *deriv_2nd = T_y;
      break;
    case T_xz:
      *deriv_1st = T_x;
      *deriv_2nd = T_z;
      break;
    case T_yz:
      *deriv_1st = T_y;
      *deriv_2nd = T_z;
      break;
    default:
      Error0(INCOMPLETE_FUNC);
  }
}

/* dc/df where c is coefficients of expansion in a direction with n nodes
// ->return value: dc(i)/df(l)
*/
static double dc_df(const Uint n,const Uint i,const Uint l)
{
  double dcdf = 0;
  const double SIGN[2] = {1.0,-1.0};
  
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
static double ChebExtrema_1point(const Uint n, const Uint p)
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
    case T_xy:
      sprintf(str,"x,y");
      break;
    case T_xz:
      sprintf(str,"x,z");
      break;  
    case T_yz:
      sprintf(str,"y,z");
      break;
    default:
      Error0(INCOMPLETE_FUNC);
  }
}

/* getting the appropriate reader for J matrices according
// to their format.
// ->return value: matrix reader for specific format
*/
fJs_T *get_j_reader(const Matrix_T *const m)
{
  fJs_T *reader = 0;
  
  if (m->ccs_f)
    reader = read_matrix_entry_ccs;
  else
    Error0(INCOMPLETE_FUNC);
  
  return reader;  
}

/* getting j_* matrix according to its type */
Matrix_T *get_j_matrix(const Patch_T *const patch,const char *type)
{
  Solving_Man_T *const sol_man = patch->solving_man;
  Matrix_T *j = 0;
  char *jtype = interpret_type(type);
  Uint i;
  
  if (!sol_man)
    return 0;
  if (!sol_man->jacobian)
    return 0;
  
  for (i = 0; i < sol_man->nj; ++i)
  {
    if (strcmp_i(sol_man->jacobian[i]->type,jtype))
    {
      j = sol_man->jacobian[i]->J;
      break;
    }
  }
  
  free(jtype);
  
  return j;
}

/* given matrix, row and column of a CCS format matrix,
// it returns the corresponing enteries of matrix.
// ->return value: m[i][j] in which m is in CCS format
*/
double read_matrix_entry_ccs(Matrix_T *const m, const long r,const long c)
{
  const int *const Ap    = m->ccs->Ap;
  const int *const Ai    = m->ccs->Ai;
  const double *const Ax = m->ccs->Ax;
  
  /* moving along none zero entries of the matrix at column c.
  // Note: it should not pass the given row.  */
  for (int i = Ap[c]; Ai[i] <= r && i < Ap[c+1]; ++i)
    if (Ai[i] == r) return Ax[i];
    
  return 0.;
}

/* calculating the give ccs matirx size.
// ->return value: size of matrix in Mb
*/
static double J_sizeMb_ccs(const Matrix_T *const m)
{
  long Uint n1,n2,n3;
  double sum = 0;
  
  if (m->ccs_f)
  {
    n1 = (long Uint)m->col;
    n2 = (long Uint)m->ccs->Ap[n1];/* number of none zero entries */
    n3 = (n1+1)*sizeof(int)/* size of Ap */ + 
         n2*sizeof(int)/* size of Ai */+
         n2*sizeof(double)/* size of aij */;
    sum = (double)n3/1E6;
  }
  else if (m->ccs_l_f)
  {
    n1 = (long Uint)m->col;
    n2 = (long Uint)m->ccs_long->Ap[n1];/* number of none zero entries */
    n3 =  (n1+1)*sizeof(int)/* size of Ap */ + 
          n2*sizeof(int)/* size of Ai */+
          n2*sizeof(double)/* size of aij */;
    sum = (double)n3/1E6;
  }
  else
  {
    Error0("This given matrix is not in CCS format.\n");
  }
  
  return sum;
}

/* supposed to write J in ccs format in disk. No completed yet! */
static void write_J_in_disk_ccs(void)
{
  Error0("Jacobian Exceeded max size and no function implemented yet."
          "To avoid this error increase \"Maximum_Size_of_J_Kept_in_Mb\" parameter.\n");
}

/* getting patch, subface and directive, it decides which d(interpolation)/df 
// be chosen and then returns the appropriate function.
// directives:
//	x derivative # calculate d(interp(f_x))/df
//	y derivative # calculate d(interp(f_y))/df
//	z derivative # calculate d(interp(f_z))/df
// 	none	     # calculate d(interp(f))/df
// ->return value: d(interp(?))/df */
fdInterp_dfs_T *get_dInterp_df(const Patch_T *const patch,const SubFace_T *const sf,const char *const dir)
{
  fdInterp_dfs_T *Func = 0;
  char type[_MAX_STR_] = {'\0'};
  Uint i;

  for (i = 0; i < 3; ++i)
  {
    if (patch->basis[i]       == Chebyshev_Tn_BASIS &&
        patch->collocation[i] == Chebyshev_Extrema)
        strcat(type,"Tn_Extrema,");
  }
  
  if (!strcmp(type,"Tn_Extrema,Tn_Extrema,Tn_Extrema,"))
  {
    if (!sf)/* if no subface is involved */
    {
      if      (!strcmp("x derivative",dir)) Func = dInterp_x_df_XYZ_Tn_Ex;
      else if (!strcmp("y derivative",dir)) Func = dInterp_y_df_XYZ_Tn_Ex;
      else if (!strcmp("z derivative",dir)) Func = dInterp_z_df_XYZ_Tn_Ex;
      else if (!strcmp("none",dir))	    Func = dInterp_df_XYZ_Tn_Ex;
      else
        Error0("No such directive defined for this function.\n");
    }
    else if (sf->sameX)
    {
      if      (!strcmp("x derivative",dir)) Func = dInterp_x_df_YZ_Tn_Ex;
      else if (!strcmp("y derivative",dir)) Func = dInterp_y_df_YZ_Tn_Ex;
      else if (!strcmp("z derivative",dir)) Func = dInterp_z_df_YZ_Tn_Ex;
      else if (!strcmp("none",dir))	    Func = dInterp_df_YZ_Tn_Ex;
      else
        Error0("No such directive defined for this function.\n");
    }
    else if (sf->sameY)
    {
      if      (!strcmp("x derivative",dir)) Func = dInterp_x_df_XZ_Tn_Ex;
      else if (!strcmp("y derivative",dir)) Func = dInterp_y_df_XZ_Tn_Ex;
      else if (!strcmp("z derivative",dir)) Func = dInterp_z_df_XZ_Tn_Ex;
      else if (!strcmp("none",dir))	    Func = dInterp_df_XZ_Tn_Ex;
      else
        Error0("No such directive defined for this function.\n");
    }
    else if (sf->sameZ)
    {
      if      (!strcmp("x derivative",dir)) Func = dInterp_x_df_XY_Tn_Ex;
      else if (!strcmp("y derivative",dir)) Func = dInterp_y_df_XY_Tn_Ex;
      else if (!strcmp("z derivative",dir)) Func = dInterp_z_df_XY_Tn_Ex;
      else if (!strcmp("none",dir))	    Func = dInterp_df_XY_Tn_Ex;
      else
        Error0("No such directive defined for this function.\n");
    }
    else if (!sf->sameX && !sf->sameY && !sf->sameZ)
    {
      if      (!strcmp("x derivative",dir)) Func = dInterp_x_df_XYZ_Tn_Ex;
      else if (!strcmp("y derivative",dir)) Func = dInterp_y_df_XYZ_Tn_Ex;
      else if (!strcmp("z derivative",dir)) Func = dInterp_z_df_XYZ_Tn_Ex;
      else if (!strcmp("none",dir))	    Func = dInterp_df_XYZ_Tn_Ex;
      else
        Error0("No such directive defined for this function.\n");
    }
    else 
    {
      Error0("Such surface's flag has not been defined for this function.\n");
    }
  }/* end of if (!strcmp(type,"Tn_Extrema,Tn_Extrema,Tn_Extrema,")) */
  else
  {
    Error0("No such option has been defined for this function.\n");
  }
  
  return Func;
}

/* d(interp(f_x))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in Y and Z direction using Cheb Tn bases
// with Extrema points.
// ->return value: d(interp(f_x))/df */
static double dInterp_x_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  const Uint X0 = plane;/* const. plane */
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  double J;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  J = dq2_dq1(patch,_N0_,_x_,i_j_k_to_ijk(n,X0,b,c));
  
  if (!EQL(J,0))
  {
    sum += J
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,q[0])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  }
  
  /* only changing of field at the plane X=X[0] must contribute */
  if (a == X0)
  {
    s = 0;
    for (r = 0; r < n[1]; ++r)
    {
      l = i_j_k_to_ijk(n,a,r,c);
      qr = General2ChebyshevExtrema(node[l]->X[1],1,patch);

      s += dq2_dq1(patch,_N1_,_x_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],r,q[1]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
    sum += s;

    s = 0;
    for (r = 0; r < n[2]; ++r)
    {
      l = i_j_k_to_ijk(n,a,b,r);
      qr = General2ChebyshevExtrema(node[l]->X[2],2,patch);

      s += dq2_dq1(patch,_N2_,_x_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],r,q[2]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
    sum += s;

  }

  return sum;
}

/* d(interp(f_y))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in Y and Z direction using Cheb Tn bases
// with Extrema points.
// ->return value: d(interp(f_y))/df */
static double dInterp_y_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  const Uint X0 = plane;;/* const. plane */
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  double J;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  J = dq2_dq1(patch,_N0_,_y_,i_j_k_to_ijk(n,X0,b,c));
  
  if (!EQL(J,0))
  {
    sum += J
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,q[0])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  }
  
  /* only changing of field at the plane X=X[0] must contribute */
  if (a == X0)
  {
    s = 0;
    for (r = 0; r < n[1]; ++r)
    {
      l = i_j_k_to_ijk(n,a,r,c);
      qr = General2ChebyshevExtrema(node[l]->X[1],1,patch);

      s += dq2_dq1(patch,_N1_,_y_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],r,q[1]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
    sum += s;

    s = 0;
    for (r = 0; r < n[2]; ++r)
    {
      l = i_j_k_to_ijk(n,a,b,r);
      qr = General2ChebyshevExtrema(node[l]->X[2],2,patch);

      s += dq2_dq1(patch,_N2_,_y_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],r,q[2]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
    sum += s;

  }

  return sum;
}

/* d(interp(f_z))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in Y and Z direction using Cheb Tn bases
// with Extrema points.
// ->return value: d(interp(f_z))/df */
static double dInterp_z_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  const Uint X0 = plane;/* const. plane */
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  double J;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  J = dq2_dq1(patch,_N0_,_z_,i_j_k_to_ijk(n,X0,b,c));
  
  if (!EQL(J,0))
  {
    sum += J
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,q[0])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  }
  
  /* only changing of field at the plane X=X[0] must contribute */
  if (a == X0)
  {
    s = 0;
    for (r = 0; r < n[1]; ++r)
    {
      l = i_j_k_to_ijk(n,a,r,c);
      qr = General2ChebyshevExtrema(node[l]->X[1],1,patch);

      s += dq2_dq1(patch,_N1_,_z_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],r,q[1]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
    sum += s;

    s = 0;
    for (r = 0; r < n[2]; ++r)
    {
      l = i_j_k_to_ijk(n,a,b,r);
      qr = General2ChebyshevExtrema(node[l]->X[2],2,patch);

      s += dq2_dq1(patch,_N2_,_z_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],r,q[2]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
    sum += s;

  }

  return sum;
}

/* d(interp(f))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in Y and Z direction using Cheb Tn bases
// with Extrema points.
// plane argument is not used; it can be put to any number!
// ->return value: d(interp(f))/df */
static double dInterp_df_YZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  const double *point = patch->node[df]->X;
  double q[3];/* normalized coords */
  Uint i,j,k;
  
  ijk_to_i_j_k(df,n,&i,&j,&k);
  q[1] = General2ChebyshevExtrema(X[1],1,patch);
  q[2] = General2ChebyshevExtrema(X[2],2,patch);
  
  /* only changing of field at the plane X=X[0] must contribute */
  if (!EQL(point[0],X[0]))
    return 0;
  
  UNUSED(plane);
  return sum_0_N_dCi_dfj_by_Ti_q(n[1],j,q[1])
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[2],k,q[2]);
}

/* d(interp(f_x))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Z direction using Cheb Tn bases
// with Extrema points.
// ->return value: d(interp(f_x))/df */
static double dInterp_x_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  const Uint Y0 = plane;/* const. plane */
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  double J;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  J = dq2_dq1(patch,_N1_,_x_,i_j_k_to_ijk(n,a,Y0,c));
  
  if (!EQL(J,0))
  {
    sum += J
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,q[1])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  }
  
  /* only changing of field at the plane X=X[1] must contribute */
  if (b == Y0)
  {
    s = 0;
    for (r = 0; r < n[0]; ++r)
    {
      l = i_j_k_to_ijk(n,r,b,c);
      qr = General2ChebyshevExtrema(node[l]->X[0],0,patch);

      s += dq2_dq1(patch,_N0_,_x_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],r,q[0]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
    sum += s;

    s = 0;
    for (r = 0; r < n[2]; ++r)
    {
      l = i_j_k_to_ijk(n,a,b,r);
      qr = General2ChebyshevExtrema(node[l]->X[2],2,patch);

      s += dq2_dq1(patch,_N2_,_x_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],r,q[2]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0]);
    sum += s;

  }

  return sum;
}

/* d(interp(f_y))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Z direction using Cheb. Tn bases
// with Extrema points.
// ->return value: d(interp(f_y))/df */
static double dInterp_y_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  const Uint Y0 = plane;/* const. plane */
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  double J;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  J = dq2_dq1(patch,_N1_,_y_,i_j_k_to_ijk(n,a,Y0,c));
  
  if (!EQL(J,0))
  {
    sum += J
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,q[1])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  }
  
  /* only changing of field at the plane X=X[1] must contribute */
  if (b == Y0)
  {
    s = 0;
    for (r = 0; r < n[0]; ++r)
    {
      l = i_j_k_to_ijk(n,r,b,c);
      qr = General2ChebyshevExtrema(node[l]->X[0],0,patch);

      s += dq2_dq1(patch,_N0_,_y_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],r,q[0]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
    sum += s;

    s = 0;
    for (r = 0; r < n[2]; ++r)
    {
      l = i_j_k_to_ijk(n,a,b,r);
      qr = General2ChebyshevExtrema(node[l]->X[2],2,patch);

      s += dq2_dq1(patch,_N2_,_y_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],r,q[2]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0]);
    sum += s;

  }

  return sum;
}

/* d(interp(f_z))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Z direction using Cheb Tn bases
// with Extrema points.
// ->return value: d(interp(f_z))/df */
static double dInterp_z_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  const Uint Y0 = plane;/* const. plane */
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  double J;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  J = dq2_dq1(patch,_N1_,_z_,i_j_k_to_ijk(n,a,Y0,c));
  
  if (!EQL(J,0))
  {
    sum += J
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,q[1])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  }
  
  /* only changing of field at the plane X=X[1] must contribute */
  if (b == Y0)
  {
    s = 0;
    for (r = 0; r < n[0]; ++r)
    {
      l = i_j_k_to_ijk(n,r,b,c);
      qr = General2ChebyshevExtrema(node[l]->X[0],0,patch);

      s += dq2_dq1(patch,_N0_,_z_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],r,q[0]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
    sum += s;

    s = 0;
    for (r = 0; r < n[2]; ++r)
    {
      l = i_j_k_to_ijk(n,a,b,r);
      qr = General2ChebyshevExtrema(node[l]->X[2],2,patch);

      s += dq2_dq1(patch,_N2_,_z_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[2],r,q[2]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0]);
    sum += s;

  }

  return sum;
}
/* d(interp(f))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Z direction using Cheb Tn bases
// with Extrema points.
// plane argument is not used; it can be put to any number!
// ->return value: d(interp(f))/df */
static double dInterp_df_XZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  const double *point = patch->node[df]->X;
  double q[3];/* normalized coords */
  Uint i,j,k;
  
  ijk_to_i_j_k(df,n,&i,&j,&k);
  q[0] = General2ChebyshevExtrema(X[0],0,patch);
  q[2] = General2ChebyshevExtrema(X[2],2,patch);
  
  
  /* only changing of field at the plane X=X[1] must contribute */
  if (!EQL(point[1],X[1]))
    return 0;
  
  UNUSED(plane);
  return sum_0_N_dCi_dfj_by_Ti_q(n[0],i,q[0])
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[2],k,q[2]);
}

/* d(interp(f_x))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Y direction using Cheb Tn bases
// with Extrema points.
// ->return value: d(interp(f_x))/df */
static double dInterp_x_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  const Uint Z0 = plane;/* const. plane */
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  double J;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  J = dq2_dq1(patch,_N2_,_x_,i_j_k_to_ijk(n,a,b,Z0));
  
  if (!EQL(J,0))
  {
    sum += J
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,q[2])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
  }
  
  /* only changing of field at the plane X=X[2] must contribute */
  if (c == Z0)
  {
    s = 0;
    for (r = 0; r < n[0]; ++r)
    {
      l = i_j_k_to_ijk(n,r,b,c);
      qr = General2ChebyshevExtrema(node[l]->X[0],0,patch);

      s += dq2_dq1(patch,_N0_,_x_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],r,q[0]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
    sum += s;

    s = 0;
    for (r = 0; r < n[1]; ++r)
    {
      l = i_j_k_to_ijk(n,a,r,c);
      qr = General2ChebyshevExtrema(node[l]->X[1],1,patch);

      s += dq2_dq1(patch,_N1_,_x_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],r,q[1]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0]);
    sum += s;

  }

  return sum;
}

/* d(interp(f_y))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Y direction using Cheb Tn bases
// with Extrema points.
// ->return value: d(interp(f_y))/df */
static double dInterp_y_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  const Uint Z0 = plane;/* const. plane */
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  double J;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  J = dq2_dq1(patch,_N2_,_y_,i_j_k_to_ijk(n,a,b,Z0));
  
  if (!EQL(J,0))
  {
    sum += J
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,q[2])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
  }
  
  /* only changing of field at the plane X=X[2] must contribute */
  if (c == Z0)
  {
    s = 0;
    for (r = 0; r < n[0]; ++r)
    {
      l = i_j_k_to_ijk(n,r,b,c);
      qr = General2ChebyshevExtrema(node[l]->X[0],0,patch);

      s += dq2_dq1(patch,_N0_,_y_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],r,q[0]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
    sum += s;

    s = 0;
    for (r = 0; r < n[1]; ++r)
    {
      l = i_j_k_to_ijk(n,a,r,c);
      qr = General2ChebyshevExtrema(node[l]->X[1],1,patch);

      s += dq2_dq1(patch,_N1_,_y_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],r,q[1]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0]);
    sum += s;

  }
  
  return sum;
}

/* d(interp(f_z))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Y direction using Cheb Tn bases
// with Extrema points.
// ->return value: d(interp(f_z))/df */
static double dInterp_z_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  const Uint Z0 = plane;/* const. plane */
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  double J;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  J = dq2_dq1(patch,_N2_,_z_,i_j_k_to_ijk(n,a,b,Z0));
  
  if (!EQL(J,0))
  {
    sum += J
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,q[2])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
  }
  
  /* only changing of field at the plane X=X[2] must contribute */
  if (c == Z0)
  {
    s = 0;
    for (r = 0; r < n[0]; ++r)
    {
      l = i_j_k_to_ijk(n,r,b,c);
      qr = General2ChebyshevExtrema(node[l]->X[0],0,patch);

      s += dq2_dq1(patch,_N0_,_z_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[0],r,q[0]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
    sum += s;

    s = 0;
    for (r = 0; r < n[1]; ++r)
    {
      l = i_j_k_to_ijk(n,a,r,c);
      qr = General2ChebyshevExtrema(node[l]->X[1],1,patch);

      s += dq2_dq1(patch,_N1_,_z_,l)
           *
           sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,qr)
           *
           sum_0_N_dCi_dfj_by_Ti_q(n[1],r,q[1]);
    }
    s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0]);
    sum += s;

  }
  
  return sum;
}

/* d(interp(f))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Y direction using Cheb Tn bases
// with Extrema points.
// plane argument is not used; it can be put to any number!
// ->return value: d(interp(f))/df */
static double dInterp_df_XY_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  const double *point = patch->node[df]->X;
  double q[3];/* normalized coords */
  Uint i,j,k;
  
  ijk_to_i_j_k(df,n,&i,&j,&k);
  q[0] = General2ChebyshevExtrema(X[0],0,patch);
  q[1] = General2ChebyshevExtrema(X[1],1,patch);
  
  /* only changing of field at the plane X=X[2] must contribute */
  if (!EQL(point[2],X[2]))
    return 0;
  
  UNUSED(plane);
  return sum_0_N_dCi_dfj_by_Ti_q(n[0],i,q[0])
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[1],j,q[1]);
}


/* d(interp(f_x))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Y and Z direction using Cheb Tn bases
// with Extrema points.
// plane argument is not used; it can be put to any number!
// ->return value: d(interp(f_x))/df */
static double dInterp_x_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  
  s = 0;
  for (r = 0; r < n[0]; ++r)
  {
    l  = i_j_k_to_ijk(n,r,b,c);
    qr = General2ChebyshevExtrema(node[l]->X[0],0,patch);
    s += dq2_dq1(patch,_N0_,_x_,l)
         *
         sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,qr)
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[0],r,q[0]);
  }
  s *= sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1])
       *
       sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  sum += s;
  
  s = 0;
  for (r = 0; r < n[1]; ++r)
  {
    l  = i_j_k_to_ijk(n,a,r,c);
    qr = General2ChebyshevExtrema(node[l]->X[1],1,patch);
    s += dq2_dq1(patch,_N1_,_x_,l)
         *
         sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,qr)
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[1],r,q[1]);
  }
  s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
       *
       sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  sum += s;
  
  s = 0;
  for (r = 0; r < n[2]; ++r)
  {
    l  = i_j_k_to_ijk(n,a,b,r);
    qr = General2ChebyshevExtrema(node[l]->X[2],2,patch);
    s += dq2_dq1(patch,_N2_,_x_,l)
         *
         sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,qr)
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[2],r,q[2]);
  }
  s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
       *
       sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
  sum += s;
  
  UNUSED(plane);
  return sum;
}

/* d(interp(f_y))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Y and Z direction using Cheb Tn bases
// with Extrema points.
// plane argument is not used; it can be put to any number!
// ->return value: d(interp(f_y))/df */
static double dInterp_y_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  
  s = 0;
  for (r = 0; r < n[0]; ++r)
  {
    l  = i_j_k_to_ijk(n,r,b,c);
    qr = General2ChebyshevExtrema(node[l]->X[0],0,patch);
    s += dq2_dq1(patch,_N0_,_y_,l)
         *
         sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,qr)
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[0],r,q[0]);
  }
  s *= sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1])
       *
       sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  sum += s;
  
  s = 0;
  for (r = 0; r < n[1]; ++r)
  {
    l  = i_j_k_to_ijk(n,a,r,c);
    qr = General2ChebyshevExtrema(node[l]->X[1],1,patch);
    s += dq2_dq1(patch,_N1_,_y_,l)
         *
         sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,qr)
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[1],r,q[1]);
  }
  s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
       *
       sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  sum += s;
  
  s = 0;
  for (r = 0; r < n[2]; ++r)
  {
    l  = i_j_k_to_ijk(n,a,b,r);
    qr = General2ChebyshevExtrema(node[l]->X[2],2,patch);
    s += dq2_dq1(patch,_N2_,_y_,l)
         *
         sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,qr)
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[2],r,q[2]);
  }
  s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
       *
       sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
  sum += s;
  
  UNUSED(plane);
  return sum;
}

/* d(interp(f_z))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Y and Z direction using Cheb Tn bases
// with Extrema points.
// plane argument is not used; it can be put to any number!
// ->return value: d(interp(f_z))/df */
static double dInterp_z_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  Node_T **const node = patch->node;
  double q[3];/* normalized coords. it is the same as N0, N1 and N2 */
  double sum = 0,s,qr;
  Uint a,b,c,r,l;
  
  ijk_to_i_j_k(df,n,&a,&b,&c);
  q[0]   = General2ChebyshevExtrema(X[0],0,patch);
  q[1]   = General2ChebyshevExtrema(X[1],1,patch);
  q[2]   = General2ChebyshevExtrema(X[2],2,patch);
  
  s = 0;
  for (r = 0; r < n[0]; ++r)
  {
    l  = i_j_k_to_ijk(n,r,b,c);
    qr = General2ChebyshevExtrema(node[l]->X[0],0,patch);
    s += dq2_dq1(patch,_N0_,_z_,l)
         *
         sum_0_N_dCi_dfj_by_dTi_dq(n[0],a,qr)
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[0],r,q[0]);
  }
  s *= sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1])
       *
       sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  sum += s;
  
  s = 0;
  for (r = 0; r < n[1]; ++r)
  {
    l  = i_j_k_to_ijk(n,a,r,c);
    qr = General2ChebyshevExtrema(node[l]->X[1],1,patch);
    s += dq2_dq1(patch,_N1_,_z_,l)
         *
         sum_0_N_dCi_dfj_by_dTi_dq(n[1],b,qr)
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[1],r,q[1]);
  }
  s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
       *
       sum_0_N_dCi_dfj_by_Ti_q(n[2],c,q[2]);
  sum += s;
  
  s = 0;
  for (r = 0; r < n[2]; ++r)
  {
    l  = i_j_k_to_ijk(n,a,b,r);
    qr = General2ChebyshevExtrema(node[l]->X[2],2,patch);
    s += dq2_dq1(patch,_N2_,_z_,l)
         *
         sum_0_N_dCi_dfj_by_dTi_dq(n[2],c,qr)
         *
         sum_0_N_dCi_dfj_by_Ti_q(n[2],r,q[2]);
  }
  s *= sum_0_N_dCi_dfj_by_Ti_q(n[0],a,q[0])
       *
       sum_0_N_dCi_dfj_by_Ti_q(n[1],b,q[1]);
  sum += s;
  
  UNUSED(plane);
  return sum;
}

/* d(interp(f))/df at point X in general coordinates.
// df is node number which we are varying field at that point.
// interpolation takes place in X and Y and Z direction using Cheb Tn bases
// with Extrema points.
// plane argument is not used; it can be put to any number!
// ->return value: d(interp(f))/df */
static double dInterp_df_XYZ_Tn_Ex(Patch_T *const patch,const double *const X,const Uint df,const Uint plane)
{
  const Uint *const n = patch->n;
  double q[3];/* normalized coords */
  Uint i,j,k;
  
  ijk_to_i_j_k(df,n,&i,&j,&k);
  q[0] = General2ChebyshevExtrema(X[0],0,patch);
  q[1] = General2ChebyshevExtrema(X[1],1,patch);
  q[2] = General2ChebyshevExtrema(X[2],2,patch);
  
  UNUSED(plane);
  return sum_0_N_dCi_dfj_by_Ti_q(n[0],i,q[0])*
         sum_0_N_dCi_dfj_by_Ti_q(n[1],j,q[1])*
         sum_0_N_dCi_dfj_by_Ti_q(n[2],k,q[2]);
  
}

/* free thoroughly patch->interface */
void free_patch_SolMan_jacobian(Patch_T *const patch)
{
  Solving_Man_T *const SolMan = patch->solving_man;
  Uint i;
  
  if (!SolMan)
    return;
  
  for (i = 0; i < SolMan->nj; ++i)
  {
    free_matrix(SolMan->jacobian[i]->J);
    Free(SolMan->jacobian[i]);
  }
  Free(SolMan->jacobian);
  
  SolMan->jacobian = 0;
  SolMan->nj       = 0;
}


