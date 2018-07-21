/*
// Alireza Rashti
// July 2018
*/

#include "tests.h"
#define ArgM(a) a,#a/*used for being more accurate in naming and fast */
#define MAXSTR 400

/* testing:
// ========
//
// given a set of analytic functions, it expands them in these basis and
// take its various derivatives ans finally 
// compares the analytic results with numeric results.
// note: the results will be printed accordingly in 
// "DerivativeTest" folder.
// note: only those patches that use basis will be compared.
// ->return value: EXIT_SUCCESS;
*/
int DerivativeTest(Grid_T *const grid)
{
  sFunc_Patch2Pdouble_T **DataBase_func;
  
  char *path;
  char der_s[MAXSTR];
  unsigned fi;
  Flag_T flg;
  
  path = get_parameter_value_S("output_directory_path",&flg);
  parameterEr(flg);

  path = dup_s(path);
  make_directory(&path,"DerivativeTest",YES);

  
  init_func_Patch2Pdouble(&DataBase_func);
  /* below is all the available analytic functions with their derivatives. 
  // one can add more function, with the same style below.
  // note: derivative must be shown by underline '_' and direction of
  // derivative as it has been shown. it is "recommended" to use the same
  // notation for naming of functions.
  // note: functions are defined in Analytic folder in Maths.
  */
/*  add_func_Patch2Pdouble(&DataBase_func,ArgM(c_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(x_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(y_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(z_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f));*/
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly3_f));
  //add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f));
/*  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(c_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(c_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(c_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(x_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(x_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(x_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(x_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(y_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(y_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(y_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(y_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(z_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(z_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(z_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(z_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinx_f_xyz));*/
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly3_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly3_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly3_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly3_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly3_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly3_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly3_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly3_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly3_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly3_f_xyz));
/*  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(r_f_xyz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xyz));
*/  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(sinxyz_f_xyz));
 /* add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(mix2_f_xyz));*/
  
  FOR_ALL(fi,DataBase_func)
  {
    /* avoid counting some of funcs which has already been counted */
    if (DataBase_func[fi]->flg == 1) continue;
    
    sFunc_Patch2Pdouble_T *F[N_FUNC];
    
    double *anac[N_FUNC];/* analytic */
    double *numc[N_FUNC];/* numeric */
    enum FUNC_E e;
    unsigned p;
    
    /* initializing, make them to point to 0 */
    for (e = FUNC; e < N_FUNC; ++e)
    {
      anac[e] = 0;
      numc[e] = 0;
    }
    
    /* read F from data base */
    flg = read_F(F,DataBase_func,fi);
    
    /* if the is no derivative for this function it continues */
    if (flg == NO) continue;
    
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      Field_T *df_num = add_field("Numerica_derivative","(3dim)",patch,NO);
      
      /* compute anac if any */
      for (e = FUNC; e < N_FUNC; ++e)
        if (F[e])  anac[e] = F[e]->func(patch);
      
      for (e = FUNC_x; e < N_FUNC; ++e)
      {
        /* if anac is define and so not null, compare with numeric */
        if (anac[e])
        {
          printf("Testing Derivative for %15s\t",F[e]->name);
          df_num->v = anac[FUNC];
          free_v2(df_num);
          enum2str(e,der_s);
          numc[e] = Df(df_num,der_s);
          compare_derivative(F[e]->name,numc[e],anac[e],patch,path);
          free(anac[e]);
          free(numc[e]);
          printf("[X] Ready.\n");
        }
      }/* end of for (e = FUNC_x; e < N_FUNC; ++e) */
      free(anac[FUNC]);
      df_num->v = 0;
      remove_field(df_num);
  
    }/* end of FOR_ALL_PATCHES(pa,grid) */
    
  }/* end of FOR_ALL(fi,DataBase_func) */
  
  free_func_Patch2Pdouble(DataBase_func);
  free(path);
  
  return EXIT_SUCCESS;
}

/* get a function and based on its name, find all of its 
// derivative in data base.
// ->return value: if related derivative found YES, NO otherwise. 
*/
static Flag_T read_F(sFunc_Patch2Pdouble_T **const F,sFunc_Patch2Pdouble_T **const DataBase_func,const enum FUNC_E fn)
{
  const char *fname = DataBase_func[fn]->name;/* function name */
  char fname_derivative[MAXSTR];/* e.g. poly3_f_xy */
  enum FUNC_E e;
  Flag_T flg = NO;
  
  F[FUNC] = get_func_Patch2Pdouble(fname,DataBase_func);
  
  if (!F[FUNC])
    abortEr_s("There is no function %s .\n",fname);
  else
    F[FUNC]->flg = 1;
    
  for (e = FUNC_x; e < N_FUNC; ++e)
  {
    sprintf(fname_derivative,"%s",fname);
    enum2strcat(e,fname_derivative);
    F[e] = get_func_Patch2Pdouble(fname_derivative,DataBase_func);
    
    /* if there is no such a derivative continue */
    if(!F[e]) continue;
    
    F[e]->flg = 1;
    flg = YES;
  }
  
  return flg;
}

/* appending a string based on the given enum e to the fname_derivatives */
static void enum2strcat(enum FUNC_E e,char *const fname_derivative)
{
  switch (e)
  {
    case FUNC_x:
      strcat(fname_derivative,"_x");
      break;
    case FUNC_y:
      strcat(fname_derivative,"_y");
      break;
    case FUNC_z:
      strcat(fname_derivative,"_z");
      break;
    case FUNC_xx:
      strcat(fname_derivative,"_xx");
      break;
    case FUNC_yy:
      strcat(fname_derivative,"_yy");
      break;
    case FUNC_zz:
      strcat(fname_derivative,"_zz");
      break;
    case FUNC_xy:
      strcat(fname_derivative,"_xy");
      break;
    case FUNC_xz:
      strcat(fname_derivative,"_xz");
      break;
    case FUNC_yz:
      strcat(fname_derivative,"_yz");
      break;
    case FUNC_xyz:
      strcat(fname_derivative,"_xyz");
      break;
    default:
      abortEr("There is no such derivative defined.\n"
      "If you added more kind of derivative please add" 
        "to enum FUNC_E and consequently other locations.\n");
  }
}

/* filling a string str based on the given enum e  */
static void enum2str(enum FUNC_E e,char *const str)
{
  switch (e)
  {
    case FUNC_x:
      sprintf(str,"x");
      break;
    case FUNC_y:
      sprintf(str,"y");
      break;
    case FUNC_z:
      sprintf(str,"z");
      break;
    case FUNC_xx:
      sprintf(str,"x,x");
      break;
    case FUNC_yy:
      sprintf(str,"y,y");
      break;
    case FUNC_zz:
      sprintf(str,"z,z");
      break;
    case FUNC_xy:
      sprintf(str,"x,y");
      break;
    case FUNC_xz:
      sprintf(str,"x,z");
      break;
    case FUNC_yz:
      sprintf(str,"y,z");
      break;
    case FUNC_xyz:
      sprintf(str,"x,y,z");
      break;
    default:
      abortEr("There is no such derivative defined.\n"
      "If you added more kind of derivative please add" 
        "to enum FUNC_E and consequently other locations.\n");
  }
}

/* comparing the values obtained from numeric and with analytic one */
static void compare_derivative(const char *const name,const double *const numc,const double *const anac,const Patch_T *const patch,const char *const path)
{
  char prefix[MAXSTR];
  
  sprintf(prefix,"%s/%s.DiffByNode",path,name);
  pr_derivatives_DiffByNode(numc,anac,patch,prefix);
}
