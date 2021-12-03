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
// NOTE: to make this faster first call second order Jacobian
// e.g. dfxy_df and then dfx_df. this is b/c second orders save 
// first order too.  */
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
  
  /* patch->solving_man must not be empty */
  assert(sol_man);
  
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
  /* IMPORTANT NOTE: ALWAYS use sol_man->nj for loop and count, 
  // since it's varied in sub call function. */
  while (types[i] != 0)
  {
    jtype = interpret_type(types[i]);
    /* check if this type has been already made then skip this */
    Flag_T flg = NONE;
    Uint c;
    for (c = 0; c < sol_man->nj; ++c)
    {
      if (strcmp_i(sol_man->jacobian[c]->type,jtype))
      {
        /* if regular cast to ccs, this happens when 
        // second order called first */
        if (sol_man->jacobian[c]->J->reg_f)
        {
          Matrix_T *J_reg         = sol_man->jacobian[c]->J;
          sol_man->jacobian[c]->J = cast_matrix_ccs(J_reg);
          free_matrix(J_reg);
           
          /* to optimize ccs reader if required */
          #ifdef CCS_READER_OPTIMIZE
          
            int Nslice = PgetiEZ("matrix_ccs_reader_split");
            Nslice = (Nslice != INT_MAX ? Nslice : 1);
            coarse_grain_Ap_ccs_matrix(sol_man->jacobian[c]->J,Nslice);
          
          #endif
        }
        flg = FOUND;
        break;
      }
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
    
    /* saving jacobian elements in ccs but first check duplication. */
    for (c = 0; c < sol_man->nj; ++c)
     if (strcmp_i(sol_man->jacobian[c]->type,jtype))
      Errors("duplicated type for '%s'!",jtype);
    
    c = sol_man->nj;
    sol_man->jacobian = 
      realloc(sol_man->jacobian,(c+1)*sizeof(*sol_man->jacobian));
    IsNull(sol_man->jacobian);
    sol_man->jacobian[c] = calloc(1,sizeof(*sol_man->jacobian[c]));
    IsNull(sol_man->jacobian[c]);
    sprintf(sol_man->jacobian[c]->type,jtype);
    sol_man->jacobian[c]->J = cast_matrix_ccs(J);
    free_matrix(J);
    
    /* to optimize ccs reader if required */
    #ifdef CCS_READER_OPTIMIZE
    
      int Nslice = PgetiEZ("matrix_ccs_reader_split");
      Nslice = (Nslice != INT_MAX ? Nslice : 1);
      coarse_grain_Ap_ccs_matrix(sol_man->jacobian[c]->J,Nslice);
    
    #endif
    
    if (GRT(J_sizeMb_ccs(sol_man->jacobian[c]->J),max_j_size))
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
          if (GRT(ABSd(cmp[Spectral_e][r][c]-cmp[Direct_e][r][c]),Err))
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
            j0 += dc_df(N[0],ip,l)*dCheb_Tn_dx((int)ip,x);
          }
          j0 *= 2;
          j0 += dc_df(N[0],ip,l)*dCheb_Tn_dx((int)ip,x);
          j0 *= cj0;
        }
      }
      
      if (l == i && n == k)
      {
        if (!EQL(cj1,0))
        {
          for (jp = 1; jp < N[1]-1; ++jp)
          {
            j1 += dc_df(N[1],jp,m)*dCheb_Tn_dx((int)jp,y);
          }
          j1 *= 2;
          j1 += dc_df(N[1],jp,m)*dCheb_Tn_dx((int)jp,y);
          j1 *= cj1;
        }
      }
      
      if (l == i && m == j)
      {
        if (!EQL(cj2,0))
        {
          for (kp = 1; kp < N[2]-1; ++kp)
          {
            j2 += dc_df(N[2],kp,n)*dCheb_Tn_dx((int)kp,z);
          }
          j2 *= 2;
          j2 += dc_df(N[2],kp,n)*dCheb_Tn_dx((int)kp,z);
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
  Solving_Man_T *const sol_man = patch->solving_man;
  double **J_1st = 0;
  Field_T *j_1st_deriv_field = 0;
  Patch_T temp_patch;
  JType_E deriv_1st = T_UNDEF,deriv_2nd = T_UNDEF;
  char deriv_2nd_s[MAX_STR_LEN/2];
  char deriv_1st_s[MAX_STR_LEN/2];
  char aux[MAX_STR_LEN] = {'\0'};
  char *jtype_1st = 0;
  Uint lmn;
  Flag_T flg = NONE;
  Uint c;
  
  read_1st_and_2nd_deriv(deriv_dir,&deriv_1st,&deriv_2nd);
  JType_E2str(deriv_2nd,deriv_2nd_s);
  JType_E2str(deriv_1st,deriv_1st_s);
  sprintf(aux,"df%s_df",deriv_1st_s);
  jtype_1st = interpret_type(aux);
  
  /* see if Jacobian exists in reg_f already so use this */
  flg = NONE;
  for (c = 0; c < sol_man->nj; ++c)
  {
   if (strcmp_i(sol_man->jacobian[c]->type,jtype_1st))
   {
    if (sol_man->jacobian[c]->J->reg_f)
    {
     flg   = FOUND;
     J_1st = sol_man->jacobian[c]->J->reg->A;
    }
    else/* it has other format like ccs */
    {
     flg = INUSE;
     J_1st = 0;
    }
    break;
   }
  }
  /* if nothing, calculate the Jacobian and save it */
  if (flg == NONE)
  {
   Matrix_T *Jm_1st = alloc_matrix(REG_SF,nn,nn);
   J_1st            = Jm_1st->reg->A;
   /* -> J = d(df/d@)/df */
   fill_jacobian_spectral_method_1stOrder(J_1st,patch,deriv_1st);
   
   /* saving jacobian elements in reg_f but first check duplication. */
   for (c = 0; c < sol_man->nj; ++c)
     if (strcmp_i(sol_man->jacobian[c]->type,jtype_1st))
      Errors("duplicated type for '%s'!",jtype_1st);
    
   c = sol_man->nj;
   /* add to struct in reg_f to be used potentially later */
   sol_man->jacobian = 
     realloc(sol_man->jacobian,(c+1)*sizeof(*sol_man->jacobian));
   IsNull(sol_man->jacobian);
   sol_man->jacobian[c] = calloc(1,sizeof(*sol_man->jacobian[c]));
   IsNull(sol_man->jacobian[c]);
   sprintf(sol_man->jacobian[c]->type,jtype_1st);
   sol_man->jacobian[c]->J = Jm_1st;
   Jm_1st = 0;
   ++sol_man->nj;
  }
  /* we need to construct again */
  else if (flg == INUSE)
  {
   J_1st = J;/* using the same give memory */
   /* -> J = d(df/d@)/df */
   fill_jacobian_spectral_method_1stOrder(J_1st,patch,deriv_1st);
  }
  /* do nothing */
  else if (flg == FOUND)
  {
   ;
  }
  else
   Error0(NO_OPTION);
  
  temp_patch = make_temp_patch(patch);
  j_1st_deriv_field = add_field("j_1st_deriv_field","(3dim)",&temp_patch,YES);
  for (lmn = 0; lmn < nn; ++lmn)
  {
    Uint ijk;
    double *j_2nd_deriv_value = 0;
   
    for (ijk = 0; ijk < nn; ++ijk)
      j_1st_deriv_field->v[ijk] = J_1st[ijk][lmn];
      
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
  Free(jtype_1st);
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
  
  if (l == 0)
    dcdf = 1;
  else if (l == n-1)
    dcdf = Jsign(i);
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

#ifdef CCS_READER_OPTIMIZE
/* given matrix, row and column of a CCS format matrix,
// it returns the corresponing enteries of matrix.
// NOTE: this function is heavily used and must be
// super optimized.
// ->return value: m[i][j] in which m is in CCS format. */
double read_matrix_entry_ccs(Matrix_T *const m, const long r,const long c)
{
  const int *const Ai    = m->ccs->Ai;
  const double *const Ax = m->ccs->Ax;
  const int *const Ap_cg = m->ccs->Ap_cg;
  const int *const i_cg  = m->ccs->i_cg;
  const int i_i_max      = Ap_cg[c+1]-1;
  int i_i                = Ap_cg[c];/* i_{i} */
   
  /* find the interval(slice) where given row resides */
  while (i_i < i_i_max)
  {
    // if (Ai[i_cg[i_i]] <= r && r <= Ai[i_cg[i_i+1]-1])
    if (r <= Ai[i_cg[i_i+1]-1])
      break;
    
    ++i_i;
  }
  
  //for (int i = i_cg[i_i]; Ai[i] <= r && i < Ap[c+1]; ++i)
  const int i_max = i_cg[i_i+1];
  for (int i = i_cg[i_i]; i < i_max; ++i)
    if (Ai[i] == r) return Ax[i];
   
  return 0.;
}

#else
/* given matrix, row and column of a CCS format matrix,
// it returns the corresponing enteries of matrix.
// NOTE: this function is heavily used and must be
// super optimized.
// ->return value: m[i][j] in which m is in CCS format. */
double read_matrix_entry_ccs(Matrix_T *const m, const long r,const long c)
{
  const int *const Ai    = m->ccs->Ai;
  const double *const Ax = m->ccs->Ax;
  const int *const Ap    = m->ccs->Ap;
  const int i_max        = Ap[c+1];
  /* moving along none zero entries of the matrix at column c.
  // Note: it should not pass the given row.  */
  for (int i = Ap[c]; Ai[i] <= r && i < i_max; ++i)
    if (Ai[i] == r) return Ax[i];
  
  return 0.;
}
#endif


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
  
  for (i = 0; i < 3; i++)
  {
    Free(SolMan->jacobian_workspace->dT_dx[i]);
    Free(SolMan->jacobian_workspace->d2T_dx2[i]);
  }
}

/* slice Ap of ccs matrix to optimiza ccs reader.
// see explanations at Matrix_T->ccs. */
static void coarse_grain_Ap_ccs_matrix(Matrix_T *const m,const int Nslice)
{
  assert(Nslice);
  
  /* alloc and init */
  m->ccs->Nslice = Nslice;
  m->ccs->Ap_cg   = calloc((Uint)m->col+1,sizeof(*m->ccs->Ap_cg));
  IsNull(m->ccs->Ap_cg);
  m->ccs->i_cg   = calloc((Uint)(m->col*Nslice+1),sizeof(*m->ccs->i_cg));
  IsNull(m->ccs->i_cg);
  
  const int *const Ap = m->ccs->Ap;
  int *const Ap_cg = m->ccs->Ap_cg;
  int *const i_cg = m->ccs->i_cg;
  int i_i;
  long c;
    
  i_i = 0;
  for (c = 0; c < m->col; ++c)
  {
    if (Ap[c+1]-Ap[c] < Nslice)
      Error0("Parameter 'matrix_ccs_reader_split' is too large!\n");

    int quotient = (Ap[c+1]-Ap[c])/(Nslice);
    
    Ap_cg[c] = i_i;
    for (int s = 0; s < Nslice; ++s)
    {
      i_cg[i_i] = Ap[c] + s*quotient;
      i_i++;
    }
  }
  i_cg[i_i] = Ap[c];
  Ap_cg[c]   = i_i;
}


/* -> d/dX 2*sum_0^N (Tn(Xj) Tn(X))| X = Xi.
// X = cos(th). */
static double
d_dXi_2xsum_0_N_Tnj_Tni(double thi/* X_i = cos(theta_i) */,
                        double thj/* X_j = cos(theta_j) */,
                        Uint N/* the sum upper limit ( patch->n-1 ) */)
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
// X = cos(th), (optimized). */
static double
d2_dXi2_2xsum_0_N_Tnj_Tni_opt(const double thi/* X_i = cos(theta_i) */,
                              const double thj/* X_i = cos(theta_i) */,
                              Patch_T *const patch,
                              Uint X_axis)

{
  double sum;
  double lambda;
  double sin_half_lambda;
  double cos_half_lambda;
  double csc_half_lambda;
  double cot_half_lambda;
  
  if (EQL(thi,0.))
  {
    double cos_lambda;
    lambda = thj;
    sin_half_lambda = sin(0.5*lambda);
    cos_half_lambda = cos(0.5*lambda);
    csc_half_lambda = 1./sin_half_lambda;
    cot_half_lambda = cos_half_lambda/sin_half_lambda;
    cos_lambda      = cos(lambda);

    sum = 
      Jd4_dlambda4_sum_0_N_cos_nlambda_opt(X_axis,JW->N0[X_axis],lambda) +
      Jd2_dlambda2_sum_0_N_cos_nlambda_opt(X_axis,JW->N0[X_axis],lambda);
    sum *= 2./3.;
  }
  else if (EQL(thi,M_PI))
  {
    double cos_lambda;
    lambda = thj+M_PI;
    sin_half_lambda = sin(0.5*lambda);
    cos_half_lambda = cos(0.5*lambda);
    csc_half_lambda = 1./sin_half_lambda;
    cot_half_lambda = cos_half_lambda/sin_half_lambda;
    cos_lambda      = cos(lambda);
    
    sum = 
      Jd4_dlambda4_sum_0_N_cos_nlambda_opt(X_axis,JW->N0[X_axis],lambda) +
      Jd2_dlambda2_sum_0_N_cos_nlambda_opt(X_axis,JW->N0[X_axis],lambda);
    
    lambda = thj-M_PI;
    sin_half_lambda = sin(0.5*lambda);
    cos_half_lambda = cos(0.5*lambda);
    csc_half_lambda = 1./sin_half_lambda;
    cot_half_lambda = cos_half_lambda/sin_half_lambda;
    cos_lambda      = cos(lambda);
    
    sum += 
      Jd4_dlambda4_sum_0_N_cos_nlambda_opt(X_axis,JW->N0[X_axis],lambda) +
      Jd2_dlambda2_sum_0_N_cos_nlambda_opt(X_axis,JW->N0[X_axis],lambda);
      
    sum /= 3.;
  }
  else
  {
    double sin_thi   = JW->sin_thi[X_axis];
    double d2thi_dX2 = -JW->cos_thi[X_axis]/(Pow3(sin_thi));
    double dthi_dX   = -1./sin_thi;
    lambda    = thi+thj;
    sin_half_lambda = sin(0.5*lambda);
    cos_half_lambda = cos(0.5*lambda);
    csc_half_lambda = 1./sin_half_lambda;
    cot_half_lambda = cos_half_lambda/sin_half_lambda;
    
    
    sum = d2thi_dX2*Jd_dlambda_sum_0_N_cos_nlambda_opt(X_axis,JW->N0[X_axis],lambda) +
          Pow2(dthi_dX)*Jd2_dlambda2_sum_0_N_cos_nlambda_opt(X_axis,JW->N0[X_axis],lambda);
    
    lambda = thi-thj;
    sin_half_lambda = sin(0.5*lambda);
    cos_half_lambda = cos(0.5*lambda);
    csc_half_lambda = 1./sin_half_lambda;
    cot_half_lambda = cos_half_lambda/sin_half_lambda;
    
    sum += d2thi_dX2*Jd_dlambda_sum_0_N_cos_nlambda_opt(X_axis,JW->N0[X_axis],lambda) +
           Pow2(dthi_dX)*Jd2_dlambda2_sum_0_N_cos_nlambda_opt(X_axis,JW->N0[X_axis],lambda);
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


/* ->: compute d(df/du)/dx, in which x is a Cartesian coords. */
double
  d2f_dxdu_spectral_Jacobian_analytic(Patch_T *const patch,
                                      const Uint dx_axis, 
                                      const Uint ijk,const Uint lmn)
{
  Uint i,j,k;
  Uint l,m,n;
  
  ijk_to_i_j_k(ijk,patch->n,&i,&j,&k);
  ijk_to_i_j_k(lmn,patch->n,&l,&m,&n);
    
  return
    Jd2f_dudx(patch,dx_axis,0,ijk,lmn,i,l)*JKD(j,m)*JKD(k,n)+
    Jd2f_dudx(patch,dx_axis,1,ijk,lmn,j,m)*JKD(i,l)*JKD(k,n)+
    Jd2f_dudx(patch,dx_axis,2,ijk,lmn,k,n)*JKD(i,l)*JKD(j,m);
}

/* ->: compute d(df/du)/dx, in which x is a Cartesian coords (optimized).
// in this optimized version, many of the quantites are saved 
// in patch->solving_man->jacobian_workspace */
double
  d2f_dxdu_optimized_spectral_Jacobian_analytic(Patch_T *const patch,
                                                const Uint dx_axis) 
{
  double sum;
  
  switch(JW->kd)
  {
    case JKD_zero:
    case JKD_il:
    case JKD_jm:
    case JKD_kn:
      sum = 0;
    break;
    
    case JKD_iljm:
      sum = Jd2f_dudx_opt(patch,dx_axis,2,JW->ijk,JW->k,JW->n);
    break;
    
    case JKD_ilkn:
      sum = Jd2f_dudx_opt(patch,dx_axis,1,JW->ijk,JW->j,JW->m);
    break;
    
    case JKD_jmkn:
      sum = Jd2f_dudx_opt(patch,dx_axis,0,JW->ijk,JW->i,JW->l);
    break;
    
    default:/* case JKD_iljmkn */
      sum = 
        Jd2f_dudx_opt(patch,dx_axis,0,JW->ijk,JW->i,JW->l)+
        Jd2f_dudx_opt(patch,dx_axis,1,JW->ijk,JW->j,JW->m)+
        Jd2f_dudx_opt(patch,dx_axis,2,JW->ijk,JW->k,JW->n);
  }
  
  return sum;
}

/* ->: compute d^2(df/du)/dxdy, in which x and y are Cartesian coords.(optimized) */
double
  d3f_dxdydu_optimized_spectral_Jacobian_analytic(Patch_T *const patch,
                                                  const Uint dxdy_axis)
{
  double sum;
  Uint dx_axis, dy_axis;
  
  /* set dx_axis & dy_axis.
  // convention for n:
  // 0 = (x,x), 1=(x,y), 2=(x,z), 3=(y,y), 4=(y,z), 5=(z,z) */
  if (dxdy_axis == 5)
  {
    dx_axis = dy_axis = 2;
  }
  else
  {
    dx_axis = dxdy_axis/3;
    dy_axis = dxdy_axis/3 + dxdy_axis%3;
  }
  
  switch(JW->kd)
  {
    case JKD_zero:
      sum = 0;
    break;
    
    case JKD_il:
      sum = Jd2f_dudx_opt(patch,dx_axis,1,JW->ijk,JW->j,JW->m)*
            Jd2f_dudx_opt(patch,dy_axis,2,JW->ijk,JW->k,JW->n)
            +
            Jd2f_dudx_opt(patch,dx_axis,2,JW->ijk,JW->k,JW->n)*
            Jd2f_dudx_opt(patch,dy_axis,1,JW->ijk,JW->j,JW->m);
    break;
    
    case JKD_jm:
      sum = Jd2f_dudx_opt(patch,dx_axis,0,JW->ijk,JW->i,JW->l)*
              (
                Jd2f_dudx_opt(patch,dy_axis,2,JW->ijk,JW->k,JW->n)
              ) +
            Jd2f_dudx_opt(patch,dx_axis,2,JW->ijk,JW->k,JW->n)*
              (
                Jd2f_dudx_opt(patch,dy_axis,0,JW->ijk,JW->i,JW->l)
              );
    break;
    
    case JKD_kn:
        sum = Jd2f_dudx_opt(patch,dx_axis,0,JW->ijk,JW->i,JW->l)*
                (
                  Jd2f_dudx_opt(patch,dy_axis,1,JW->ijk,JW->j,JW->m)
                ) +
              Jd2f_dudx_opt(patch,dx_axis,1,JW->ijk,JW->j,JW->m)*
                (
                  Jd2f_dudx_opt(patch,dy_axis,0,JW->ijk,JW->i,JW->l)
                );
    break;
    
    case JKD_iljm:
      sum = Jd2f_dudx_opt(patch,dx_axis,0,JW->ijk,JW->i,JW->l)*
              (
                Jd2f_dudx_opt(patch,dy_axis,2,JW->ijk,JW->k,JW->n)
              ) +
            Jd2f_dudx_opt(patch,dx_axis,1,JW->ijk,JW->j,JW->m)*
              (
                Jd2f_dudx_opt(patch,dy_axis,2,JW->ijk,JW->k,JW->n)
              ) +
            Jd3f_dudxdy_opt(patch,dx_axis,dy_axis,dxdy_axis,2,JW->ijk,JW->k,JW->n) +
            Jd2f_dudx_opt(patch,dx_axis,2,JW->ijk,JW->k,JW->n)*
              (
                Jd2f_dudx_opt(patch,dy_axis,1,JW->ijk,JW->j,JW->m) +
                Jd2f_dudx_opt(patch,dy_axis,0,JW->ijk,JW->i,JW->l)
              );
    break;
    
    case JKD_ilkn:
      sum = Jd2f_dudx_opt(patch,dx_axis,0,JW->ijk,JW->i,JW->l)*
              (
                Jd2f_dudx_opt(patch,dy_axis,1,JW->ijk,JW->j,JW->m)
              ) +
            Jd3f_dudxdy_opt(patch,dx_axis,dy_axis,dxdy_axis,1,JW->ijk,JW->j,JW->m) +
            Jd2f_dudx_opt(patch,dx_axis,1,JW->ijk,JW->j,JW->m)*
              (
                Jd2f_dudx_opt(patch,dy_axis,0,JW->ijk,JW->i,JW->l) +
                Jd2f_dudx_opt(patch,dy_axis,2,JW->ijk,JW->k,JW->n)
              ) +
            Jd2f_dudx_opt(patch,dx_axis,2,JW->ijk,JW->k,JW->n)*
              (
                Jd2f_dudx_opt(patch,dy_axis,1,JW->ijk,JW->j,JW->m)
              );
    break;
    
    case JKD_jmkn:
      sum = Jd3f_dudxdy_opt(patch,dx_axis,dy_axis,dxdy_axis,0,JW->ijk,JW->i,JW->l) +
            Jd2f_dudx_opt(patch,dx_axis,0,JW->ijk,JW->i,JW->l)*
              (
                Jd2f_dudx_opt(patch,dy_axis,1,JW->ijk,JW->j,JW->m) +
                Jd2f_dudx_opt(patch,dy_axis,2,JW->ijk,JW->k,JW->n)
              ) +
            Jd2f_dudx_opt(patch,dx_axis,1,JW->ijk,JW->j,JW->m)*
              (
                Jd2f_dudx_opt(patch,dy_axis,0,JW->ijk,JW->i,JW->l)
              ) +
            Jd2f_dudx_opt(patch,dx_axis,2,JW->ijk,JW->k,JW->n)*
              (
                Jd2f_dudx_opt(patch,dy_axis,0,JW->ijk,JW->i,JW->l)
              );
    break;
    
    default:/* case JKD_iljmkn */
      sum = Jd3f_dudxdy_opt(patch,dx_axis,dy_axis,dxdy_axis,0,JW->ijk,JW->i,JW->l)+
            Jd2f_dudx_opt(patch,dx_axis,0,JW->ijk,JW->i,JW->l)*
              (
                Jd2f_dudx_opt(patch,dy_axis,1,JW->ijk,JW->j,JW->m) +
                Jd2f_dudx_opt(patch,dy_axis,2,JW->ijk,JW->k,JW->n)
              ) +
            Jd3f_dudxdy_opt(patch,dx_axis,dy_axis,dxdy_axis,1,JW->ijk,JW->j,JW->m) +
            Jd2f_dudx_opt(patch,dx_axis,1,JW->ijk,JW->j,JW->m)*
              (
                Jd2f_dudx_opt(patch,dy_axis,0,JW->ijk,JW->i,JW->l) +
                Jd2f_dudx_opt(patch,dy_axis,2,JW->ijk,JW->k,JW->n)
              ) +
            Jd3f_dudxdy_opt(patch,dx_axis,dy_axis,dxdy_axis,2,JW->ijk,JW->k,JW->n) +
            Jd2f_dudx_opt(patch,dx_axis,2,JW->ijk,JW->k,JW->n)*
              (
                Jd2f_dudx_opt(patch,dy_axis,1,JW->ijk,JW->j,JW->m) +
                Jd2f_dudx_opt(patch,dy_axis,0,JW->ijk,JW->i,JW->l)
              );
  }
  
  
  return sum;
  
  /*  for reference:
    sum = 
    Jd3f_dudxdy_opt(patch,dx_axis,dy_axis,dxdy_axis,0,JW->ijk,JW->i,JW->l)*JKD(JW->j,JW->m)*JKD(JW->k,JW->n) +
    Jd2f_dudx_opt(patch,dx_axis,0,JW->ijk,JW->i,JW->l)*
      (
        JKD(JW->k,JW->n)*Jd2f_dudx_opt(patch,dy_axis,1,JW->ijk,JW->j,JW->m) +
        JKD(JW->j,JW->m)*Jd2f_dudx_opt(patch,dy_axis,2,JW->ijk,JW->k,JW->n)
      ) +
    
    Jd3f_dudxdy_opt(patch,dx_axis,dy_axis,dxdy_axis,1,JW->ijk,JW->j,JW->m)*JKD(JW->i,JW->l)*JKD(JW->k,JW->n) +
    Jd2f_dudx_opt(patch,dx_axis,1,JW->ijk,JW->j,JW->m)*
      (
        JKD(JW->k,JW->n)*Jd2f_dudx_opt(patch,dy_axis,0,JW->ijk,JW->i,JW->l) +
        JKD(JW->i,JW->l)*Jd2f_dudx_opt(patch,dy_axis,2,JW->ijk,JW->k,JW->n)
      ) +
      
    Jd3f_dudxdy_opt(patch,dx_axis,dy_axis,dxdy_axis,2,JW->ijk,JW->k,JW->n)*JKD(JW->j,JW->m)*JKD(JW->i,JW->l) +
    Jd2f_dudx_opt(patch,dx_axis,2,JW->ijk,JW->k,JW->n)*
      (
        JKD(JW->i,JW->l)*Jd2f_dudx_opt(patch,dy_axis,1,JW->ijk,JW->j,JW->m) +
        JKD(JW->j,JW->m)*Jd2f_dudx_opt(patch,dy_axis,0,JW->ijk,JW->i,JW->l)
      ); */
}

/* ->: compute d^2(df/du)/dxdy, in which x and y are Cartesian coords. */
double
  d3f_dxdydu_spectral_Jacobian_analytic(Patch_T *const patch,
                                        const int dxdy_axis,
                                        const Uint ijk,const Uint lmn)
{
  int dx_axis, dy_axis;
  Uint i,j,k;
  Uint l,m,n;
  
  /* set dx_axis & dy_axis.
  // convention for n:
  // 0 = (x,x), 1=(x,y), 2=(x,z), 3=(y,y), 4=(y,z), 5=(z,z) */
  if (dxdy_axis == 5)
  {
    dx_axis = dy_axis = 2;
  }
  else
  {
    dx_axis = dxdy_axis/3;
    dy_axis = dxdy_axis/3 + dxdy_axis%3;
  }
  
  //printf("dxdy = %d, dx = %d, dy = %d\n",dxdy_axis,dx_axis,dy_axis);
  
  ijk_to_i_j_k(ijk,patch->n,&i,&j,&k);
  ijk_to_i_j_k(lmn,patch->n,&l,&m,&n);
  
  return
    Jd3f_dudxdy(patch,dx_axis,dy_axis,dxdy_axis,0,ijk,lmn,i,l)*JKD(j,m)*JKD(k,n) +
    Jd2f_dudx(patch,dx_axis,0,ijk,lmn,i,l)*
      (
        JKD(k,n)*Jd2f_dudx(patch,dy_axis,1,ijk,lmn,j,m) +
        JKD(j,m)*Jd2f_dudx(patch,dy_axis,2,ijk,lmn,k,n)
      ) +
    
    Jd3f_dudxdy(patch,dx_axis,dy_axis,dxdy_axis,1,ijk,lmn,j,m)*JKD(i,l)*JKD(k,n) +
    Jd2f_dudx(patch,dx_axis,1,ijk,lmn,j,m)*
      (
        JKD(k,n)*Jd2f_dudx(patch,dy_axis,0,ijk,lmn,i,l) +
        JKD(i,l)*Jd2f_dudx(patch,dy_axis,2,ijk,lmn,k,n)
      ) +
      
    Jd3f_dudxdy(patch,dx_axis,dy_axis,dxdy_axis,2,ijk,lmn,k,n)*JKD(j,m)*JKD(i,l) +
    Jd2f_dudx(patch,dx_axis,2,ijk,lmn,k,n)*
      (
        JKD(i,l)*Jd2f_dudx(patch,dy_axis,1,ijk,lmn,j,m) +
        JKD(j,m)*Jd2f_dudx(patch,dy_axis,0,ijk,lmn,i,l)
      );
}


/* set Solving_Man_T->jacobian_workspace.
// NOTE: if the flag "set" is set, it won't populate the jacobian_workspace again. */
void set_Solving_Man_jacobian_workspace(Patch_T *const patch)
{
  Solving_Man_T *const solving_man = patch->solving_man;
  
  /* some checks */
  if (!solving_man) return;
  if (solving_man->jacobian_workspace->set) return;
  
  Uint *const nm1 = solving_man->jacobian_workspace->nm1;
  double *const pi_o_nm1 = solving_man->jacobian_workspace->pi_o_nm1;
  double *const norm = solving_man->jacobian_workspace->norm;
  double *const N0   = solving_man->jacobian_workspace->N0;
  double *const c1_d2 = solving_man->jacobian_workspace->c1_d2;
  double *const c2_d2 = solving_man->jacobian_workspace->c2_d2;
  double *const c1_d4 = solving_man->jacobian_workspace->c1_d4;
  double *const c2_d4 = solving_man->jacobian_workspace->c2_d4;
  double *const c3_d4 = solving_man->jacobian_workspace->c3_d4;
  double *const c4_d4 = solving_man->jacobian_workspace->c4_d4;
  double *const c5_d4 = solving_man->jacobian_workspace->c5_d4;
  Uint i;
  
  for (i = 0; i < 3; ++i)
  {
    nm1[i]      = patch->n[i]-1;
    pi_o_nm1[i] = M_PI/nm1[i];
    norm[i] = 0.5/nm1[i];
    N0[i]   = 0.5 + nm1[i];
    c1_d2[i] = -(Pow3(nm1[i])/3.+Pow2(nm1[i])/2.+nm1[i]/6.);
    c2_d2[i] = -1. - 4.*Pow2(N0[i]);
    c1_d4[i] = (Pow4(nm1[i])*(nm1[i]/5.+0.5)+Pow3(nm1[i])/3.-nm1[i]/30.);
    c2_d4[i] = 4.*Pow2(N0[i]);
    c3_d4[i] = 115. - 120.*Pow2(N0[i]) + 48.*Pow4(N0[i]);
    c4_d4[i] = (76. + 96.*Pow2(N0[i]) - 64.*Pow4(N0[i]));
    c5_d4[i] = (1. + 24.*Pow2(N0[i]) + 16.*Pow4(N0[i]));
  }
  
  /* set d^n Cheb/dx^n */
  double *dT_dx[3]   = {0};
  double *d2T_dx2[3] = {0};
  
  dT_dx[0] = alloc_double(patch->nn);
  dT_dx[1] = alloc_double(patch->nn);
  dT_dx[2] = alloc_double(patch->nn);
  
  d2T_dx2[0] = alloc_double(patch->nn);
  d2T_dx2[1] = alloc_double(patch->nn);
  d2T_dx2[2] = alloc_double(patch->nn);
  
  /* set */
  for (Uint ijk = 0; ijk < patch->nn; ++ijk)
  {
    Uint ip,jp,kp;
    double x[3];
    
    ijk_to_i_j_k(ijk,patch->n,&ip,&jp,&kp);
    x[0] =  cos(ip*pi_o_nm1[0]);
    x[1] =  cos(jp*pi_o_nm1[1]);
    x[2] =  cos(kp*pi_o_nm1[2]);
    
    dT_dx[0][ijk] = dCheb_Tn_dx((int)(nm1[0]),x[0]);
    dT_dx[1][ijk] = dCheb_Tn_dx((int)(nm1[1]),x[1]);
    dT_dx[2][ijk] = dCheb_Tn_dx((int)(nm1[2]),x[2]);
    
    d2T_dx2[0][ijk] = d2Cheb_Tn_dx2((int)(nm1[0]),x[0]);
    d2T_dx2[1][ijk] = d2Cheb_Tn_dx2((int)(nm1[1]),x[1]);
    d2T_dx2[2][ijk] = d2Cheb_Tn_dx2((int)(nm1[2]),x[2]);
  }
  
  /* save */
  for (i = 0; i < 3; ++i)
  {
    solving_man->jacobian_workspace->dT_dx[i]   = dT_dx[i];
    solving_man->jacobian_workspace->d2T_dx2[i] = d2T_dx2[i];
  }
  
  /* fully set */
  solving_man->jacobian_workspace->set = 1;
}
