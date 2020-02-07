/*
// Alireza Rashti
// Feb 2018
*/

#include "solvings_tests.h"

/* testing root finders */
void test_root_finders(Grid_T *const grid)
{
  if (DO)
  {
    printf("\nRoot Finder test: Steepest Descent method:\n\n");
    root_finder_SteepestDescent(grid);
  }

}

/* testing steepest Descent root finder
// ->return value: EXIT_SUCCESS */
static int root_finder_SteepestDescent(Grid_T *const grid)
{
  double (*f0)(void *params,const double *const x) = root_finder_f0_eq;
  double (*f1)(void *params,const double *const x) = root_finder_f1_eq;
  double (*f2)(void *params,const double *const x) = root_finder_f2_eq;
  double (*df0_dx)(void *params,const double *const x,const unsigned dir) = root_finder_df0_dx_eq;
  double (*df1_dx)(void *params,const double *const x,const unsigned dir) = root_finder_df1_dx_eq;
  double (*df2_dx)(void *params,const double *const x,const unsigned dir) = root_finder_df2_dx_eq;
  Root_Finder_T *root;
  const double x_analytic[3] = {0.,0.1,1.};
  double *x_sol;
  const char *const par = Pgets_E("Test_RootFinders");
  unsigned MaxIter = 20;
  
  if (regex_search("[[:digit:]]+",par))
  {
    char *s = regex_find("[[:digit:]]+",par);
    MaxIter = (unsigned)atoi(s);
    _free(s);
  }
  
  /* testing with the derivatives are given: */
  root = init_root_finder(3);
  root->verbose       = 1;
  root->type          = "Steepest_Descent";
  plan_root_finder(root);
  root->description = "solving f = 0, knowing df_dx";
  root->tolerance   = 10E-15;
  root->MaxIter     = MaxIter;
  root->f[0]        = f0;
  root->f[1]        = f1;
  root->f[2]        = f2;
  root->df_dx[0]    = df0_dx;
  root->df_dx[1]    = df1_dx;
  root->df_dx[2]    = df2_dx;
  x_sol             = execute_root_finder(root);
  print_root_finder_exit_status(root);
  
  /* printing solutions: */
  printf("->Steepest Descent Method root finder found (derivatives were given):\n");
  printf("x0  = %+e, x1  = %+e, x2  = %+e\n",x_sol[0],x_sol[1],x_sol[2]);
  printf("dx0 = %+e, dx1 = %+e, dx2 = %+e\n",
         fabs(x_sol[0]-x_analytic[0]),fabs(x_sol[1]-x_analytic[1]),fabs(x_sol[2]-x_analytic[2]));
  printf("Residual = %e\n\n",root->residual);
  free_root_finder(root);
  free(x_sol);
  
  /* testing without the derivatives are given: */
  root = init_root_finder(3);
  root->verbose       = 1;
  root->type          = "Steepest_Descent";
  plan_root_finder(root);
  root->description = "solving f = 0, df_dx is calculated by finite difference";
  root->tolerance   = 10E-15;
  root->MaxIter     = MaxIter;
  root->f[0]        = f0;
  root->f[1]        = f1;
  root->f[2]        = f2;
  x_sol             = execute_root_finder(root);
  /* printing solutions: */
  printf("->Steepest Descent Method root finder found (No derivatives are given):\n");
  printf("x0  = %+e, x1  = %+e, x2  = %+e\n",x_sol[0],x_sol[1],x_sol[2]);
  printf("dx0 = %+e, dx1 = %+e, dx2 = %+e\n",
         fabs(x_sol[0]-x_analytic[0]),fabs(x_sol[1]-x_analytic[1]),fabs(x_sol[2]-x_analytic[2]));
  printf("Residual = %e\n\n",root->residual);
  free_root_finder(root);
  free(x_sol);
  
  return EXIT_SUCCESS;
  UNUSED(grid);
}

/* ->return value: d(root_finder_f0_eq)/dx^{dir} */
static double root_finder_df0_dx_eq(void *params,const double *const x,const unsigned dir)
{
  double df_dx = 0;
  
  if (dir == 0)
  {
    df_dx = 1 - x[1]*x[2]*sin(x[0]*x[1]*x[2]);
  }
  else if (dir == 1)
  {
    df_dx = - x[0]*x[2]*sin(x[0]*x[1]*x[2]);
  }
  else if (dir == 2)
  {
    df_dx = - x[0]*x[1]*sin(x[0]*x[1]*x[2]);
  }
  else
    abortEr("Bad argument for derivative.\n");
  
  return df_dx;
  UNUSED(params);
}

/* ->return value: d(root_finder_f1_eq)/dx^{dir} */
static double root_finder_df1_dx_eq(void *params,const double *const x,const unsigned dir)
{
  double df_dx = 0;
  
  if (dir == 0)
  {
    df_dx = -0.25*pow(1-x[0],-0.75);
  }
  else if (dir == 1)
  {
    df_dx = 1;
  }
  else if (dir == 2)
  {
    df_dx = 0.1*x[2]-0.15;
  }
  else
    abortEr("Bad argument for derivative.\n");
  
  return df_dx;
  UNUSED(params);
}

/* ->return value: d(root_finder_f2_eq)/dx^{dir} */
static double root_finder_df2_dx_eq(void *params,const double *const x,const unsigned dir)
{
  double df_dx = 0;
  
  if (dir == 0)
  {
    df_dx = -2*x[0];
  }
  else if (dir == 1)
  {
    df_dx = -0.2 *x[1]+0.01;
  }
  else if (dir == 2)
  {
    df_dx = 1;
  }
  else
    abortEr("Bad argument for derivative.\n");
  
  return df_dx;
  UNUSED(params);
}

/* ->return value: x0+cos(x0 x1 x2) -1 */
static double root_finder_f0_eq(void *params,const double *const x)
{
  return x[0]+cos(x[0]*x[1]*x[2])-1;
  UNUSED(params);
}

/* ->return value: (1-x0)^(1/4) + x1 + 0.05*x2^2-0.15*x2-1 */
static double root_finder_f1_eq(void *params,const double *const x)
{
  return pow(1-x[0],0.25)+x[1]+0.05*SQR(x[2])-0.15*x[2]-1;
  UNUSED(params);
}

/* ->return value: -x0^2-0.1 x1^2 + 0.01 x1 + x2 -1 */
static double root_finder_f2_eq(void *params,const double *const x)
{
  return -SQR(x[0])-0.1*SQR(x[1])+0.01*x[1]+x[2]-1;
  UNUSED(params);
}

/* testing if the value of dfs_df are correct */
void test_dfs_df_values(Grid_T *const grid)
{
  const char *const types[] = {"dfx_df","dfy_df","dfz_df",
                               "dfxx_df","dfxy_df","dfxz_df",
                               "dfyy_df","dfyz_df","dfzz_df",0};
  const double start = get_time_sec();
  test_make_Js_jacobian_eq(grid,types);
  pr_spent_time(start,"Making Jacobian");
}

/* testing various d(Interpolation)/df */
void test_dInterp_a_df(Grid_T *const grid)
{
  unsigned p;
  const unsigned Num_Tests = 3;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    add_field("phi_field",0,patch,YES);
    Field_T *phi_field = patch->pool[Ind("phi_field")];
    double *phi = phi_field->v;
    unsigned i;
    
    /* initialize phi */
    FOR_ALL_POINTS(i,patch)
      phi[i] = sin(y_(i))*SQR(x_(i))+cos(z_(i))*SQR(y_(i))+SQR(z_(i))*SQR(y_(i));
    
    /* to test differet random points the tests placed in a loop */
    for (i = 0; i < Num_Tests; ++i)
    {
      if (DO)
        test_dInterp_x_df_YZ_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_y_df_YZ_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_z_df_YZ_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_df_YZ_Tn_Ex(phi_field);
        
      if (DO)
        test_dInterp_x_df_XZ_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_y_df_XZ_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_z_df_XZ_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_df_XZ_Tn_Ex(phi_field);
      
      if (DO)
        test_dInterp_x_df_XY_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_y_df_XY_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_z_df_XY_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_df_XY_Tn_Ex(phi_field);
        
      if (DO)
        test_dInterp_x_df_XYZ_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_y_df_XYZ_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_z_df_XYZ_Tn_Ex(phi_field);
      if (DO)
        test_dInterp_df_XYZ_Tn_Ex(phi_field);
    }/* end of for (i = 0; i < Num_Tests; ++i) */
    remove_field(phi_field);
  }
}

/* test function: dInterp_x_df_YZ_Tn_Ex */
static void test_dInterp_x_df_YZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_x",0,patch,NO);
  Field_T *phi_field_x = patch->pool[Ind("phi_field_x")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  unsigned a,b,plane,l;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned df;

  phi_field_x->v = Partial_Derivative(phi_field,"x");
  l    = (unsigned)floor(random_double(0,nn-1,0));
  IJK(l,n,&plane,&a,&b);
  X[0] = patch->node[l]->X[0];
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field_x;  
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);

  sf->sameX = 1;
  sf->sameY = 0;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"x derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"x");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,plane);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }

  free_coeffs(phi_field);
  remove_field(phi_field_x);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_x_df_YZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_x_df_YZ_Tn_Ex function: [+].\n",patch->name);
}
  
/* test function: dInterp_y_df_YZ_Tn_Ex */
static void test_dInterp_y_df_YZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_y",0,patch,NO);
  Field_T *phi_field_y = patch->pool[Ind("phi_field_y")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  unsigned a,b,plane,l;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  phi_field_y->v = Partial_Derivative(phi_field,"y");
  l    = (unsigned)floor(random_double(0,nn-1,0));
  IJK(l,n,&plane,&a,&b);
  X[0] = patch->node[l]->X[0];
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field_y;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 1;
  sf->sameY = 0;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"y derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"y");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,plane);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  remove_field(phi_field_y);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_y_df_YZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_y_df_YZ_Tn_Ex function: [+].\n",patch->name);
}

/* test function: dInterp_z_df_YZ_Tn_Ex */
static void test_dInterp_z_df_YZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_z",0,patch,NO);
  Field_T *phi_field_z = patch->pool[Ind("phi_field_z")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  unsigned a,b,plane,l;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  phi_field_z->v = Partial_Derivative(phi_field,"z");
  l    = (unsigned)floor(random_double(0,nn-1,0));
  IJK(l,n,&plane,&a,&b);
  X[0] = patch->node[l]->X[0];
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field_z;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 1;
  sf->sameY = 0;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"z derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"z");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,plane);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  remove_field(phi_field_z);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_z_df_YZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_z_df_YZ_Tn_Ex function: [+].\n",patch->name);
}

/* test function: dInterp_df_YZ_Tn_Ex */
static void test_dInterp_df_YZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned df,l;
    
  l    = (unsigned)floor(random_double(0,nn-1,0));
  X[0] = patch->node[l]->X[0];
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 1;
  sf->sameY = 0;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"none");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,0);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_df_YZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_df_YZ_Tn_Ex function: [+].\n",patch->name);
}


/* test function: dInterp_x_df_XZ_Tn_Ex */
static void test_dInterp_x_df_XZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_x",0,patch,NO);
  Field_T *phi_field_x = patch->pool[Ind("phi_field_x")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned a,b,plane,l;
  unsigned df;

  phi_field_x->v = Partial_Derivative(phi_field,"x");
  l    = (unsigned)floor(random_double(0,nn-1,0));
  IJK(l,n,&a,&plane,&b);
  X[1] = patch->node[l]->X[1];
  X[0] = random_double(patch->min[0],patch->max[0],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field_x;  
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);

  sf->sameX = 0;
  sf->sameY = 1;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"x derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"x");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,plane);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }

  free_coeffs(phi_field);
  remove_field(phi_field_x);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_x_df_XZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_x_df_XZ_Tn_Ex function: [+].\n",patch->name);
}
  
/* test function: dInterp_y_df_XZ_Tn_Ex */
static void test_dInterp_y_df_XZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_y",0,patch,NO);
  Field_T *phi_field_y = patch->pool[Ind("phi_field_y")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned a,b,plane,l;
  unsigned df;
    
  phi_field_y->v = Partial_Derivative(phi_field,"y");
  l    = (unsigned)floor(random_double(0,nn-1,0));
  IJK(l,n,&a,&plane,&b);
  X[1] = patch->node[l]->X[1];
  X[0] = random_double(patch->min[0],patch->max[0],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field_y;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 0;
  sf->sameY = 1;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"y derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"y");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,plane);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  remove_field(phi_field_y);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_y_df_XZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_y_df_XZ_Tn_Ex function: [+].\n",patch->name);
}

/* test function: dInterp_z_df_XZ_Tn_Ex */
static void test_dInterp_z_df_XZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_z",0,patch,NO);
  Field_T *phi_field_z = patch->pool[Ind("phi_field_z")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned a,b,plane,l;
  unsigned df;
    
  phi_field_z->v = Partial_Derivative(phi_field,"z");
  l    = (unsigned)floor(random_double(0,nn-1,0));
  IJK(l,n,&a,&plane,&b);
  X[1] = patch->node[l]->X[1];
  X[0] = random_double(patch->min[0],patch->max[0],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field_z;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 0;
  sf->sameY = 1;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"z derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"z");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,plane);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  remove_field(phi_field_z);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_z_df_XZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_z_df_XZ_Tn_Ex function: [+].\n",patch->name);
}

/* test function: dInterp_df_XZ_Tn_Ex */
static void test_dInterp_df_XZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned df,l;
  
  l    = (unsigned)floor(random_double(0,nn-1,0));
  X[1] = patch->node[l]->X[1];
  X[0] = random_double(patch->min[0],patch->max[0],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 0;
  sf->sameY = 1;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"none");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,0);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_df_XZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_df_XZ_Tn_Ex function: [+].\n",patch->name);
}


/* test function: dInterp_x_df_XY_Tn_Ex */
static void test_dInterp_x_df_XY_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_x",0,patch,NO);
  Field_T *phi_field_x = patch->pool[Ind("phi_field_x")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned a,b,plane,l;
  unsigned df;

  phi_field_x->v = Partial_Derivative(phi_field,"x");
  l    = (unsigned)floor(random_double(0,nn-1,0));
  IJK(l,n,&a,&b,&plane);
  X[2] = patch->node[l]->X[2];
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[0] = random_double(patch->min[0],patch->max[0],1);
  interp_s->field = phi_field_x;  
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);

  sf->sameX = 0;
  sf->sameY = 0;
  sf->sameZ = 1;
  dInterp_spec = get_dInterp_df(patch,sf,"x derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"x");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,plane);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }

  free_coeffs(phi_field);
  remove_field(phi_field_x);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_x_df_XY_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_x_df_XY_Tn_Ex function: [+].\n",patch->name);
}
  
/* test function: dInterp_y_df_XY_Tn_Ex */
static void test_dInterp_y_df_XY_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_y",0,patch,NO);
  Field_T *phi_field_y = patch->pool[Ind("phi_field_y")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned a,b,plane,l;
  unsigned df;
    
  phi_field_y->v = Partial_Derivative(phi_field,"y");
  l    = (unsigned)floor(random_double(0,nn-1,0));
  IJK(l,n,&a,&b,&plane);
  X[2] = patch->node[l]->X[2];
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[0] = random_double(patch->min[0],patch->max[0],1);
  interp_s->field = phi_field_y;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 0;
  sf->sameY = 0;
  sf->sameZ = 1;
  dInterp_spec = get_dInterp_df(patch,sf,"y derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"y");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,plane);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  remove_field(phi_field_y);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_y_df_XY_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_y_df_XY_Tn_Ex function: [+].\n",patch->name);
}

/* test function: dInterp_z_df_XY_Tn_Ex */
static void test_dInterp_z_df_XY_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_z",0,patch,NO);
  Field_T *phi_field_z = patch->pool[Ind("phi_field_z")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned a,b,plane,l;
  unsigned df;
    
  phi_field_z->v = Partial_Derivative(phi_field,"z");
  l    = (unsigned)floor(random_double(0,nn-1,0));
  IJK(l,n,&a,&b,&plane);
  X[2] = patch->node[l]->X[2];
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[0] = random_double(patch->min[0],patch->max[0],1);
  interp_s->field = phi_field_z;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 0;
  sf->sameY = 0;
  sf->sameZ = 1;
  dInterp_spec = get_dInterp_df(patch,sf,"z derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"z");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,plane);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  remove_field(phi_field_z);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_z_df_XY_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_z_df_XY_Tn_Ex function: [+].\n",patch->name);
}

/* test function: dInterp_df_XY_Tn_Ex */
static void test_dInterp_df_XY_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  X[2] = patch->node[(unsigned)floor(random_double(0,nn-1,0))]->X[2];
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[0] = random_double(patch->min[0],patch->max[0],1);
  interp_s->field = phi_field;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 0;
  sf->sameY = 0;
  sf->sameZ = 1;
  dInterp_spec = get_dInterp_df(patch,sf,"none");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,0);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_df_XY_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_df_XY_Tn_Ex function: [+].\n",patch->name);
}

/* test function: dInterp_x_df_XYZ_Tn_Ex */
static void test_dInterp_x_df_XYZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_x",0,patch,NO);
  Field_T *phi_field_x = patch->pool[Ind("phi_field_x")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned df;

  phi_field_x->v = Partial_Derivative(phi_field,"x");
  X[0] = random_double(patch->min[0],patch->max[0],0);
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field_x;  
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);

  sf->sameX = 0;
  sf->sameY = 0;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"x derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"x");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,0);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }

  free_coeffs(phi_field);
  remove_field(phi_field_x);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_x_df_XYZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_x_df_XYZ_Tn_Ex function: [+].\n",patch->name);
}
  
/* test function: dInterp_y_df_XYZ_Tn_Ex */
static void test_dInterp_y_df_XYZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_y",0,patch,NO);
  Field_T *phi_field_y = patch->pool[Ind("phi_field_y")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  phi_field_y->v = Partial_Derivative(phi_field,"y");
  X[0] = random_double(patch->min[0],patch->max[0],0);
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field_y;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 0;
  sf->sameY = 0;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"y derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"y");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,0);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  remove_field(phi_field_y);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_y_df_XYZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_y_df_XYZ_Tn_Ex function: [+].\n",patch->name);
}

/* test function: dInterp_z_df_XYZ_Tn_Ex */
static void test_dInterp_z_df_XYZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  add_field("phi_field_z",0,patch,NO);
  Field_T *phi_field_z = patch->pool[Ind("phi_field_z")];
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  phi_field_z->v = Partial_Derivative(phi_field,"z");
  X[0] = random_double(patch->min[0],patch->max[0],0);
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field_z;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 0;
  sf->sameY = 0;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"z derivative");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    free(interp_s->field->v);
    interp_s->field->v = Partial_Derivative(phi_field,"z");
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,0);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  remove_field(phi_field_z);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_z_df_XYZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_z_df_XYZ_Tn_Ex function: [+].\n",patch->name);
}

/* test function: dInterp_df_XYZ_Tn_Ex */
static void test_dInterp_df_XYZ_Tn_Ex(Field_T *const phi_field)
{
  Patch_T *const patch = phi_field->patch;
  Flag_T flg = NONE;
  SubFace_T sf[1] = {0};
  fdInterp_dfs_T *dInterp_spec = 0;
  Interpolation_T *interp_s = init_interpolation();
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const double E_i = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);/* interpolation error */
  const double EPS = E_i > 1.0 ? E_i*nn : 1.0;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  X[0] = random_double(patch->min[0],patch->max[0],0);
  X[1] = random_double(patch->min[1],patch->max[1],1);
  X[2] = random_double(patch->min[2],patch->max[2],1);
  interp_s->field = phi_field;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  f1 = execute_interpolation(interp_s);
  
  sf->sameX = 0;
  sf->sameY = 0;
  sf->sameZ = 0;
  dInterp_spec = get_dInterp_df(patch,sf,"none");
  for (df = 0; df < nn; ++df)
  {
    phi_field->v[df] += EPS;
    free_coeffs(phi_field);
    free_coeffs(interp_s->field);
    f2 = execute_interpolation(interp_s);
    phi_field->v[df] -= EPS;
    
    spec_cal = dInterp_spec(patch,X,df,0);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),E_i))
    {
      printf("spec=%0.15f,direct=%0.15f\n",spec_cal,(f2-f1)/EPS);
      flg = FOUND;
      break;
    }
  }
  
  free_coeffs(phi_field);
  free_interpolation(interp_s);
  
  if (flg == FOUND)
    printf("%s: Testing dInterp_df_XYZ_Tn_Ex function: [-].\n",patch->name);
  else
    printf("%s: Testing dInterp_df_XYZ_Tn_Ex function: [+].\n",patch->name);
}
