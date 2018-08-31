/*
// Alireza Rashti
// July 2018
*/

#include "approximation_tests.h"
#define ArgM(a) a,#a/*used for being more accurate in naming and fast */
#define MAXSTR 400

/* testing interpolation functions.
// takes bunch of random points intside each patch and compares
// the value of interpolation and exact value for a field.
// ->return value: EXIT_SUCCESS
*/
int interpolation_tests(Grid_T *const grid)
{
  unsigned p;
  double *x,*y,*z;
  Field_T *field;
  int status;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    const unsigned *n = patch->n;
    const double *min = patch->min;
    const double *max = patch->max;
    
    /* making a field which have analytic properties */
    field = add_field("interpolant","(3dim)",patch,NO);
    field->v = poly5_f(patch);/* a polynomial */
    
    /* making 2n random number in each direction */
    x = make_random_number(2*n[0],min[0],max[0]);
    y = make_random_number(2*n[1],min[1],max[1]);
    z = make_random_number(2*n[2],min[2],max[2]);
    
    if (DO)
    {
      printf("Interpolation test:      x direction, patch %10s =>",patch->name);
      status = interpolation_tests_x(field,x,2*n[0]);
      check_test_result(status);
    }
    
    if (DO)
    {
      printf("Interpolation test:      y direction, patch %10s =>",patch->name);
      status = interpolation_tests_y(field,y,2*n[1]);
      check_test_result(status);
    }
    if (DO)
    {
      printf("Interpolation test:      z direction, patch %10s =>",patch->name);
      status = interpolation_tests_z(field,z,2*n[2]);
      check_test_result(status);
    }
    if (DO)
    {
      printf("Interpolation test: x & y directions, patch %10s =>",patch->name);
      status = interpolation_tests_xy(field,x,y,2*n[0],2*n[1]);
      check_test_result(status);
    }
    if (DO)
    {
      printf("Interpolation test: x & z directions, patch %10s =>",patch->name);
      status = interpolation_tests_xz(field,x,z,2*n[0],2*n[2]);
      check_test_result(status);
    } 
    if (DO)
    {
      printf("Interpolation test: y & z directions, patch %10s =>",patch->name);
      status = interpolation_tests_yz(field,y,z,2*n[1],2*n[2]);
      check_test_result(status);
    }
    if (DO)
    {
      printf("Interpolation test:              3-D, patch %10s =>",patch->name);
      status = interpolation_tests_xyz(field,x,y,z,2*n[0],2*n[1],2*n[2]);
      check_test_result(status);
    }
    
    /* freeing */
    free(x);
    free(y);
    free(z);
    remove_field(field);
  }
  
  return EXIT_SUCCESS;
}

/* testing interpolation in x direction and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_x(Field_T *const field,const double *const x,const unsigned N)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const unsigned *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  double tol = 0;
  double y,z;
  const char *const par = GetParameterS("Test_Interpolation");
  char *tok,*save;
  Flag_T flg;
  unsigned i,j;
  
  if (patch->coordsys != Cartesian)
    abortEr("In order to use this function in different coordinates,\n"
      "one needs to consider that analytic function calculates only in \n"
      "Cartesian coordinate system, so the transformation between\n"
      "Curvilinear and Cartesian for random point must be done.\n");
  
  /* find the tolerance defined in input file if any */
  assert(par);
  tok = dup_s(par);
  tok_s(tok,',',&save);
  tol = fabs(atof(save));
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->X_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (j = 0; j < N; ++j)/* -> choose different slices for testing */
  {
    interp_s->J = (unsigned) floor(random_double(0,n[1],j));
    interp_s->K = (unsigned) floor(random_double(0,n[2],1));
    
    y = node[L(n,0,interp_s->J,0)]->x[1];
    z = node[L(n,0,0,interp_s->K)]->x[2];
    
    for (i = 0; i < N; ++i)
    {
      interp_s->X = x[i];
      diff = poly5_f_point(x[i],y,z)-execute_interpolation(interp_s);
      if (GRT(fabs(diff),tol))
      {
        flg = FOUND;
        break;
      }
    }
    
    if (flg == FOUND)
      break;
  }
  
  free_interpolation(interp_s);
  free(tok);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in y direction and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_y(Field_T *const field,const double *const y,const unsigned N)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const unsigned *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  double tol = 0;
  double x,z;
  const char *const par = GetParameterS("Test_Interpolation");
  char *tok,*save;
  Flag_T flg;
  unsigned i,j;
  
  if (patch->coordsys != Cartesian)
    abortEr("In order to use this function in different coordinates,\n"
      "one needs to consider that analytic function calculates only in \n"
      "Cartesian coordinate system, so the transformation between\n"
      "Curvilinear and Cartesian for random point must be done.\n");
  
  /* find the tolerance defined in input file if any */
  assert(par);
  tok = dup_s(par);
  tok_s(tok,',',&save);
  tol = fabs(atof(save));
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->Y_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (j = 0; j < N; ++j)/* -> choose different slices for testing */
  {
    interp_s->I = (unsigned) floor(random_double(0,n[0],j));
    interp_s->K = (unsigned) floor(random_double(0,n[2],1));
    
    x = node[L(n,interp_s->I,0,0)]->x[0];
    z = node[L(n,0,0,interp_s->K)]->x[2];
    
    for (i = 0; i < N; ++i)
    {
      interp_s->Y = y[i];
      diff = poly5_f_point(x,y[i],z)-execute_interpolation(interp_s);
      if (GRT(fabs(diff),tol))
      {
        flg = FOUND;
        break;
      }
    }
    
    if (flg == FOUND)
      break;
  }
  
  free_interpolation(interp_s);
  free(tok);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in z direction and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_z(Field_T *const field,const double *const z,const unsigned N)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const unsigned *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  double tol = 0;
  double y,x;
  const char *const par = GetParameterS("Test_Interpolation");
  char *tok,*save;
  Flag_T flg;
  unsigned i,j;
  
  if (patch->coordsys != Cartesian)
    abortEr("In order to use this function in different coordinates,\n"
      "one needs to consider that analytic function calculates only in \n"
      "Cartesian coordinate system, so the transformation between\n"
      "Curvilinear and Cartesian for random point must be done.\n");
  
  /* find the tolerance defined in input file if any */
  assert(par);
  tok = dup_s(par);
  tok_s(tok,',',&save);
  tol = fabs(atof(save));
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->Z_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (j = 0; j < N; ++j)/* -> choose different slices for testing */
  {
    interp_s->J = (unsigned) floor(random_double(0,n[1],j));
    interp_s->I = (unsigned) floor(random_double(0,n[0],1));
    
    y = node[L(n,0,interp_s->J,0)]->x[1];
    x = node[L(n,interp_s->I,0,0)]->x[0];
    
    for (i = 0; i < N; ++i)
    {
      interp_s->Z = z[i];
      diff = poly5_f_point(x,y,z[i])-execute_interpolation(interp_s);
      if (GRT(fabs(diff),tol))
      {
        flg = FOUND;
        break;
      }
    }
    
    if (flg == FOUND)
      break;
  }
  
  free_interpolation(interp_s);
  free(tok);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in x&y directions and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_xy(Field_T *const field,const double *const x,const double *const y,const unsigned Nx,const unsigned Ny)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const unsigned *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  double tol = 0;
  double z;
  const char *const par = GetParameterS("Test_Interpolation");
  char *tok,*save;
  Flag_T flg;
  unsigned a,b,c;
  
  if (patch->coordsys != Cartesian)
    abortEr("In order to use this function in different coordinates,\n"
      "one needs to consider that analytic function calculates only in \n"
      "Cartesian coordinate system, so the transformation between\n"
      "Curvilinear and Cartesian for random point must be done.\n");
  
  /* find the tolerance defined in input file if any */
  assert(par);
  tok = dup_s(par);
  tok_s(tok,',',&save);
  tol = fabs(atof(save));
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->XY_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (a = 0; a < Nx; ++a)/* -> choose different slices for testing */
  {
    interp_s->K = (unsigned) floor(random_double(0,n[2],a));
    
    z = node[L(n,0,0,interp_s->K)]->x[2];
    
    for (b = 0; b < Nx; ++b)
    {
      interp_s->X = x[b];
      for (c = 0; c < Ny; ++c)
      {
        interp_s->Y = y[c];
        diff = poly5_f_point(x[b],y[c],z)-execute_interpolation(interp_s);
        if (GRT(fabs(diff),tol))
        {
          flg = FOUND;
          break;
        }
      }
    
      if (flg == FOUND)
        break;
    }
    if (flg == FOUND)
        break;
  }
  
  free_interpolation(interp_s);
  free(tok);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in x&z directions and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_xz(Field_T *const field,const double *const x,const double *const z,const unsigned Nx,const unsigned Nz)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const unsigned *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  double tol = 0;
  double y;
  const char *const par = GetParameterS("Test_Interpolation");
  char *tok,*save;
  Flag_T flg;
  unsigned a,b,c;
  
  if (patch->coordsys != Cartesian)
    abortEr("In order to use this function in different coordinates,\n"
      "one needs to consider that analytic function calculates only in \n"
      "Cartesian coordinate system, so the transformation between\n"
      "Curvilinear and Cartesian for random point must be done.\n");
  
  /* find the tolerance defined in input file if any */
  assert(par);
  tok = dup_s(par);
  tok_s(tok,',',&save);
  tol = fabs(atof(save));
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->XZ_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (a = 0; a < Nx; ++a)/* -> choose different slices for testing */
  {
    interp_s->J = (unsigned) floor(random_double(0,n[1],a));
    
    y = node[L(n,0,interp_s->J,0)]->x[1];
    
    for (b = 0; b < Nx; ++b)
    {
      interp_s->X = x[b];
      for (c = 0; c < Nz; ++c)
      {
        interp_s->Z = z[c];
        diff = poly5_f_point(x[b],y,z[c])-execute_interpolation(interp_s);
        if (GRT(fabs(diff),tol))
        {
          flg = FOUND;
          break;
        }
      }
    
      if (flg == FOUND)
        break;
    }
    if (flg == FOUND)
        break;
  }
  
  free_interpolation(interp_s);
  free(tok);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in y&z directions and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_yz(Field_T *const field,const double *const y,const double *const z,const unsigned Ny,const unsigned Nz)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  const unsigned *const n = patch->n;
  Node_T **const node = patch->node;
  double diff = 0;
  double tol = 0;
  double x;
  const char *const par = GetParameterS("Test_Interpolation");
  char *tok,*save;
  Flag_T flg;
  unsigned a,b,c;
  
  if (patch->coordsys != Cartesian)
    abortEr("In order to use this function in different coordinates,\n"
      "one needs to consider that analytic function calculates only in \n"
      "Cartesian coordinate system, so the transformation between\n"
      "Curvilinear and Cartesian for random point must be done.\n");
  
  /* find the tolerance defined in input file if any */
  assert(par);
  tok = dup_s(par);
  tok_s(tok,',',&save);
  tol = fabs(atof(save));
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->YZ_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (a = 0; a < Ny; ++a)/* -> choose different slices for testing */
  {
    interp_s->I = (unsigned) floor(random_double(0,n[0],a));
    
    x = node[L(n,interp_s->I,0,0)]->x[0];
    
    for (b = 0; b < Ny; ++b)
    {
      interp_s->Y = y[b];
      for (c = 0; c < Nz; ++c)
      {
        interp_s->Z = z[c];
        diff = poly5_f_point(x,y[b],z[c])-execute_interpolation(interp_s);
        if (GRT(fabs(diff),tol))
        {
          flg = FOUND;
          break;
        }
      }
    
      if (flg == FOUND)
        break;
    }
    if (flg == FOUND)
        break;
  }
  
  free_interpolation(interp_s);
  free(tok);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}

/* testing interpolation in x&y&z directions and comparing
// to analytical value and returning the result.
// ->return value: result of test.
*/
static int interpolation_tests_xyz(Field_T *const field,const double *const x,const double *const y,const double *const z,const unsigned Nx,const unsigned Ny,const unsigned Nz)
{
  Interpolation_T *interp_s = init_interpolation();
  Patch_T *const patch = field->patch;
  double diff = 0;
  double tol = 0;
  const char *const par = GetParameterS("Test_Interpolation");
  char *tok,*save;
  Flag_T flg;
  unsigned a,b,c;
  
  if (patch->coordsys != Cartesian)
    abortEr("In order to use this function in different coordinates,\n"
      "one needs to consider that analytic function calculates only in \n"
      "Cartesian coordinate system, so the transformation between\n"
      "Curvilinear and Cartesian for random point must be done.\n");
  
  /* find the tolerance defined in input file if any */
  assert(par);
  tok = dup_s(par);
  tok_s(tok,',',&save);
  tol = fabs(atof(save));
  
  /* setting up interpolation */
  interp_s->field = field;
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  
  flg = NONE;
  for (a = 0; a < Nx; ++a)/* -> choose different slices for testing */
  {
    interp_s->X = x[a];
    for (b = 0; b < Ny; ++b)
    {
      interp_s->Y = y[b];
      for (c = 0; c < Nz; ++c)
      {
        interp_s->Z = z[c];
        diff = poly5_f_point(x[a],y[b],z[c])-execute_interpolation(interp_s);
        if (GRT(fabs(diff),tol))
        {
          flg = FOUND;
          break;
        }
      }
    
      if (flg == FOUND)
        break;
    }
    if (flg == FOUND)
        break;
  }
  
  free_interpolation(interp_s);
  free(tok);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}

/* testing:
// ========
//
// given a set of analytic functions, it expands them in these basis and
// take its various derivatives ans finally 
// compares the analytic results with numeric results.
// note: the results will be printed accordingly in 
// "Derivative_Tests" folder.
// note: only those patches that use basis will be compared.
// ->return value: EXIT_SUCCESS;
*/
int derivative_tests(Grid_T *const grid)
{
  sFunc_Patch2Pdouble_T **DataBase_func;
  
  const char *path_par;
  char *path;
  char der_s[MAXSTR];
  unsigned fi;
  Flag_T flg;
  
  path_par = GetParameterS_E("output_directory_path");
  path = make_directory(path_par,"Derivative_Tests");

  
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
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f));
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
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_x));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_y));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_z));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_xx));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_yy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_zz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_xy));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_xz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_yz));
  add_func_Patch2Pdouble(&DataBase_func,ArgM(poly5_f_xyz));
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
      Field_T *df_num = add_field("$Numerica_derivative_TEST_FUNCTION$","(3dim)",patch,NO);
      
      /* compute anac if any */
      for (e = FUNC; e < N_FUNC; ++e)
        if (F[e])  anac[e] = F[e]->func(patch);
      
      for (e = FUNC_x; e < N_FUNC; ++e)
      {
        /* if anac is define and so not null, compare with numeric */
        if (anac[e])
        {
          printf("Testing Derivative: patch=%s, function:%15s\t",patch->name,F[e]->name);
          df_num->v = anac[FUNC];
          free_v2(df_num);
          enum2str(e,der_s);
          numc[e] = Partial_Derivative(df_num,der_s);
          flg = compare_derivative(F[e]->name,numc[e],anac[e],patch,path);
          free(anac[e]);
          free(numc[e]);
          
          if (flg == YES)
            printf("[+].\n");
          else
            printf("[-].\n");
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
  char fname_derivative[MAXSTR];/* e.g. poly5_f_xy */
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
static Flag_T compare_derivative(const char *const name,const double *const numc,const double *const anac,const Patch_T *const patch,const char *const path)
{
  char prefix[MAXSTR];
  Flag_T flg;
  
  sprintf(prefix,"%s/%s.DiffByNode",path,name);
  flg = pr_derivatives_DiffByNode(numc,anac,patch,prefix);
  
  return flg;
}

/* making N random number in the within min and max.
// ->return value: random number(s)
*/
static double *make_random_number(const unsigned N,const double min,const double max)
{
  double *x = calloc(N,sizeof(*x));
  unsigned i;
  
  for (i = 0; i < N; ++i)
  {
    x[i] = random_double(min,max,i);
  }
  
  return x;
}
