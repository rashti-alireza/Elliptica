/*
// Alireza Rashti
// July 2018
*/

#include "tests.h"
#define ArgM(a) a,#a/*used for being more accurate in naming and fast */
#define MAXSTR 200

/* testing:
// ========
//
// given a set of analytic functions, it expands them in these basis and
// take its various derivatives ans finally 
// compares the analytic results with numeric results.
// note: the results will be printed accordingly in 
// "ExpantionTests_Derivatives" folder.
// note: only those patches that use basis will be compared.
// ->return value: EXIT_SUCCESS;
*/
int DerivativeTest(Grid_T *const grid)
{
  sFunc_Grid2Pdouble_T **DataBase_func;
  unsigned fi;
  
  init_func_Grid2Pdouble(&DataBase_func);
  /* below is all the available analytic functions with their derivatives. 
  // one can add more function, with the same style below.
  // note: derivative must be shown by underline '_' and direction of
  // derivative as it has been shown. it is "recommended" to use the same
  // notation for naming of functions.
  // note: functions are defined in Analytic folder in Maths.
  */
  add_func_Grid2Pdouble(&DataBase_func,ArgM(c_f));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(x_f));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(y_f));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(z_f));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(poly3_f));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(r_f));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(sinxyz_f));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(mix2_f));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(c_f_x));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(c_f_y));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(c_f_z));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(x_f_x));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(x_f_y));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(x_f_z));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(x_f_xx));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(y_f_x));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(y_f_y));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(y_f_z));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(y_f_yy));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(z_f_x));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(z_f_y));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(z_f_z));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(z_f_zz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(poly3_f_x));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(poly3_f_y));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(poly3_f_z));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(poly3_f_xx));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(poly3_f_yy));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(poly3_f_zz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(poly3_f_xy));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(poly3_f_xz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(poly3_f_yz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(poly3_f_xyz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(r_f_x));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(r_f_y));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(r_f_z));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(r_f_xx));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(r_f_yy));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(r_f_zz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(r_f_xy));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(r_f_xz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(r_f_yz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(r_f_xyz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(inv_rP1_f));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(inv_rP1_f_x));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(inv_rP1_f_y));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(inv_rP1_f_z));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xx));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(inv_rP1_f_yy));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(inv_rP1_f_zz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xy));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(inv_rP1_f_yz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(inv_rP1_f_xyz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(sinxyz_f_x));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(sinxyz_f_y));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(sinxyz_f_z));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(sinxyz_f_xx));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(sinxyz_f_yy));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(sinxyz_f_zz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(sinxyz_f_xy));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(sinxyz_f_xz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(sinxyz_f_yz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(sinxyz_f_xyz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(mix2_f_x));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(mix2_f_y));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(mix2_f_z));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(mix2_f_xx));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(mix2_f_yy));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(mix2_f_zz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(mix2_f_xy));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(mix2_f_xz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(mix2_f_yz));
  add_func_Grid2Pdouble(&DataBase_func,ArgM(mix2_f_xyz));
  
  FOR_ALL(fi,DataBase_func)
  {
    /* avoid counting some of funcs which has already been counted */
    if (DataBase_func[fi]->flg == 1) continue;
    
    sFunc_Grid2Pdouble_T *F[N_FUNC];
    double *anal[N_FUNC];/* analytic */
    double *numc[N_FUNC];/* numeric */
    enum FUNC_E e;
    Flag_T flg;
    
    /* initializing, make them to point to 0 */
    for (e = FUNC; e < N_FUNC; ++e)
    {
      anal[e] = 0;
      numc[e] = 0;
    }
    
    /* read F from data base */
    flg = read_F(F,DataBase_func,fi);
    
    /* if the is no derivative for this function it continues */
    if (flg == NO) continue;
    
    /* compute anal if any */
    for (e = FUNC; e < N_FUNC; ++e)
      if (F[e])  anal[e] = F[e]->func(grid);
    
    for (e = FUNC_x; e < N_FUNC; ++e)
    {
      /* if anal is define and so not null, compare with numeric */
      if (anal[e])
      {
        //test
        printf("%s\n",F[e]->name);
        //end
        //numc[e] = derivative(anal[FUNC],e);
        //analyze(numc[e],anal[e],grid);
      }
    }
    
    /* freeing */
    for (e = FUNC; e < N_FUNC; ++e)
    {
      if (anal[e]) free(anal[e]);
      //if (numc[e]) free(numc[e]);
    }
  }/* end of FOR_ALL(fi,DataBase_func) */
  
  free_func_Grid2Pdouble(DataBase_func);
  
  return EXIT_SUCCESS;
}

/* get a function and based on its name, find all of its 
// derivative in data base.
// ->return value: if related derivative found YES, NO otherwise. 
*/
static Flag_T read_F(sFunc_Grid2Pdouble_T **const F,sFunc_Grid2Pdouble_T **const DataBase_func,const enum FUNC_E fn)
{
  const char *fname = DataBase_func[fn]->name;/* function name */
  char fname_derivative[MAXSTR];/* e.g. poly3_f_xy */
  enum FUNC_E e;
  Flag_T flg = NO;
  
  F[FUNC] = get_func_Grid2Pdouble(fname,DataBase_func);
  
  if (!F[FUNC])
    abortEr_s("There is no function %s .\n",fname);
  else
    F[FUNC]->flg = 1;
    
  for (e = FUNC_x; e < N_FUNC; ++e)
  {
    sprintf(fname_derivative,"%s",fname);
    enum2strcat(e,fname_derivative);
    F[e] = get_func_Grid2Pdouble(fname_derivative,DataBase_func);
    
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

/* testing: Chebyshev first kind basis. */
static void ChebyshevFirstKindBasis_DerivativeTest(const Patch_T *const patch)
{
}