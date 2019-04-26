/*
// Alireza Rashti
// Feb 2018
*/

#include "solvings_tests.h"

/* testing if the value of dfs_df are correct */
void test_dfs_df_values(Grid_T *const grid)
{
  const char *const types[] = {"dfx_df","dfxx_df","dfy_df","dfyy_df","dfz_df","dfzz_df",0};
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
      if (NOT_DO)
        test_dInterp_y_df_YZ_Tn_Ex(phi_field);
      if (NOT_DO)
        test_dInterp_z_df_YZ_Tn_Ex(phi_field);
      if (NOT_DO)
        test_dInterp_df_YZ_Tn_Ex(phi_field);
        
      if (NOT_DO)
        test_dInterp_x_df_XZ_Tn_Ex(phi_field);
      if (NOT_DO)
        test_dInterp_y_df_XZ_Tn_Ex(phi_field);
      if (NOT_DO)
        test_dInterp_z_df_XZ_Tn_Ex(phi_field);
      if (NOT_DO)
        test_dInterp_df_XZ_Tn_Ex(phi_field);
      
      if (NOT_DO)
        test_dInterp_x_df_XY_Tn_Ex(phi_field);
      if (NOT_DO)
        test_dInterp_y_df_XY_Tn_Ex(phi_field);
      if (NOT_DO)
        test_dInterp_z_df_XY_Tn_Ex(phi_field);
      if (NOT_DO)
        test_dInterp_df_XY_Tn_Ex(phi_field);
        
      if (NOT_DO)
        test_dInterp_x_df_XYZ_Tn_Ex(phi_field);
      if (NOT_DO)
        test_dInterp_y_df_XYZ_Tn_Ex(phi_field);
      if (NOT_DO)
        test_dInterp_z_df_XYZ_Tn_Ex(phi_field);
      if (NOT_DO)
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  const double interp_error = n[0]*n[1]*n[2]*spectral_derivative_max_error(phi_field,1);
  double X[3],f1,f2,spec_cal;
  unsigned df;

  phi_field_x->v = Partial_Derivative(phi_field,"x");
  X[0] = patch->node[(unsigned)floor(random_double(0,n[0]-1,0))]->X[0];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (GRT(fabs(spec_cal-(f2-f1)/EPS),interp_error))
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  phi_field_y->v = Partial_Derivative(phi_field,"y");
  X[0] = patch->node[(unsigned)floor(random_double(0,n[0]-1,0))]->X[0];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  phi_field_z->v = Partial_Derivative(phi_field,"z");
  X[0] = patch->node[(unsigned)floor(random_double(0,n[0]-1,0))]->X[0];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  X[0] = patch->node[(unsigned)floor(random_double(0,n[0]-1,0))]->X[0];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  double X[3],f1,f2,spec_cal;
  unsigned df;

  phi_field_x->v = Partial_Derivative(phi_field,"x");
  X[1] = patch->node[(unsigned)floor(random_double(0,n[1]-1,0))]->X[1];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  phi_field_y->v = Partial_Derivative(phi_field,"y");
  X[1] = patch->node[(unsigned)floor(random_double(0,n[1]-1,0))]->X[1];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  phi_field_z->v = Partial_Derivative(phi_field,"z");
  X[1] = patch->node[(unsigned)floor(random_double(0,n[1]-1,0))]->X[1];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  X[1] = patch->node[(unsigned)floor(random_double(0,n[1]-1,0))]->X[1];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  double X[3],f1,f2,spec_cal;
  unsigned df;

  phi_field_x->v = Partial_Derivative(phi_field,"x");
  X[2] = patch->node[(unsigned)floor(random_double(0,n[2]-1,0))]->X[2];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  phi_field_y->v = Partial_Derivative(phi_field,"y");
  X[2] = patch->node[(unsigned)floor(random_double(0,n[2]-1,0))]->X[2];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  phi_field_z->v = Partial_Derivative(phi_field,"z");
  X[2] = patch->node[(unsigned)floor(random_double(0,n[2]-1,0))]->X[2];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const unsigned *n = patch->n;
  const double CONST = 1.;
  const double EPS = CONST/nn;
  double X[3],f1,f2,spec_cal;
  unsigned df;
    
  X[2] = patch->node[(unsigned)floor(random_double(0,n[2]-1,0))]->X[2];
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const double CONST = 1.;
  const double EPS = CONST/nn;
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const double CONST = 1.;
  const double EPS = CONST/nn;
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const double CONST = 1.;
  const double EPS = CONST/nn;
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
  const double CONST = 1.;
  const double EPS = CONST/nn;
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
    
    spec_cal = dInterp_spec(patch,X,df);
    if (!EQL(spec_cal,(f2-f1)/EPS))
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
