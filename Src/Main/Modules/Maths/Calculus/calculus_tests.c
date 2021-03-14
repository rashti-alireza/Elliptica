/*
// Alireza Rashti
// June 2019
*/

#include "calculus_tests.h"

/* test integral */
int integration_tests(Grid_T *const grid)
{
  int status;
  
  if (DO)
  {
    printf("\nIntegration test: Composite Simpson's Rule 1D => \n");
    status = csr_1d(grid);
    check_test_result(status);
  }
  if (DO)
  {
    printf("\nIntegration test: Gaussian Quadrature Chebyshev Extrema: \n");
    GQ_ChebExtrema(grid);
  }
  if (DO)
  {
    printf("\nIntegration test: Gaussian Quadrature Lobatto method: \n");
    GQ_Lobatto(grid);
  }
  if (DO)
  {
    printf("\nIntegration test: Gaussian Quadrature Legendre method: \n");
    GQ_Legendre(grid);
  }
  if (DO)
  {
    printf("\nIntegration test: Integral{f(x)dV}, Spectral method: \n\n");
    fdV_spectral(grid);
  }
  if (DO)
  {
    printf("\nIntegration test: Integral{f(x)dS}, Spectral method: \n\n");
    fdS_spectral(grid);
  }
  
  return EXIT_SUCCESS;
}

/* testing Integral{fdS}, Spectral method.
// ->return value: TEST_SUCCESSFUL */
static int fdS_spectral(Grid_T *const grid)
{
  Integration_T *I;
  Field_T *f;
  Patch_T *patch;
  Uint nn,ijk,p;
  const Uint *n;
  double analytic = 0,numeric = 0;
  
  /* to test this function we use f = 1, so it means we calculate the volume */
  if (grid->kind == Grid_Box)
  {
    analytic = 0;
    numeric  = 0;
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      I       = init_integration();
      I->type = "Integral{f(x)dS},Spectral";
      nn      = patch->nn;
      n       = patch->n;
      f       = add_field("int f",0,patch,YES);
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      double x1 = patch->min[0];
      double x2 = patch->max[0];
      double y1 = patch->min[1];
      double y2 = patch->max[1];
      double z1 = patch->min[2];
      double z2 = patch->max[2];
      
      for (ijk = 0; ijk < nn; ++ijk)
      {
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = pow(y_(ijk), 2)*exp(x_(ijk))*pow(sin(z_(ijk)), 2);
      }
      
      /* analytic calculation */
      analytic += (1.0/3.0)*(pow(y1, 3)*(exp(x1) - exp(x2)) + pow(y2, 3)*(-exp(x1) + exp(x2)))*pow(sin(z1), 2);
      analytic += (1.0/3.0)*(pow(y1, 3)*(exp(x1) - exp(x2)) + pow(y2, 3)*(-exp(x1) + exp(x2)))*pow(sin(z2), 2);
      analytic += (1.0/2.0)*pow(y1, 2)*((-z1 + (1.0/2.0)*sin(2*z1))*exp(x2) + (z1 - 1.0/2.0*sin(2*z1))*exp(x1) + (-z2 + (1.0/2.0)*sin(2*z2))*exp(x1) + (z2 - 1.0/2.0*sin(2*z2))*exp(x2));
      analytic += (1.0/2.0)*pow(y2, 2)*((-z1 + (1.0/2.0)*sin(2*z1))*exp(x2) + (z1 - 1.0/2.0*sin(2*z1))*exp(x1) + (-z2 + (1.0/2.0)*sin(2*z2))*exp(x1) + (z2 - 1.0/2.0*sin(2*z2))*exp(x2));
      analytic += (1.0/6.0)*(pow(y1, 3)*(z1 - 1.0/2.0*sin(2*z1)) + pow(y1, 3)*(-z2 + (1.0/2.0)*sin(2*z2)) + pow(y2, 3)*(-z1 + (1.0/2.0)*sin(2*z1)) + pow(y2, 3)*(z2 - 1.0/2.0*sin(2*z2)))*exp(x1);
      analytic += (1.0/6.0)*(pow(y1, 3)*(z1 - 1.0/2.0*sin(2*z1)) + pow(y1, 3)*(-z2 + (1.0/2.0)*sin(2*z2)) + pow(y2, 3)*(-z1 + (1.0/2.0)*sin(2*z1)) + pow(y2, 3)*(z2 - 1.0/2.0*sin(2*z2)))*exp(x2);
      
      /* lets's do numeric calculation */  
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;
      
      /* for X = const */
      I->Spectral->X_surface = 1;
      I->Spectral->I         = 0; 
      plan_integration(I);
      numeric += execute_integration(I);
      
      I->Spectral->X_surface = 1;
      I->Spectral->I         = n[0]-1; 
      plan_integration(I);
      numeric += execute_integration(I);

      /* for Y = const */
      I->Spectral->X_surface = 0;
      I->Spectral->Y_surface = 1;
      I->Spectral->J         = 0; 
      plan_integration(I);
      numeric += execute_integration(I);
      
      I->Spectral->Y_surface = 1;
      I->Spectral->J         = n[1]-1; 
      plan_integration(I);
      numeric += execute_integration(I);

      /* for Z = const */
      I->Spectral->Y_surface = 0;
      I->Spectral->Z_surface = 1;
      I->Spectral->K         = 0; 
      plan_integration(I);
      numeric += execute_integration(I);
      
      I->Spectral->Z_surface = 1;
      I->Spectral->K         = n[2]-1; 
      plan_integration(I);
      numeric += execute_integration(I);
      
      free_integration(I);
      remove_field(f);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("Cartesian grid:\n");
    printf("numeric = %e, analytic = %e, difference = %e\n",
            numeric,analytic,numeric-analytic);
    
  }
  else if (grid->kind == Grid_SplitCubedSpherical_BHNS)
  {
    double r;
    printf("BHNS_SplitCubedSpherical_grid:\n");
    
    // testing NS sphere surface:
    printf("\n--> Integral{f(x)dS}|at NS surface section:\n");
    
    r = Pgetd("NS_radius");
    analytic = 4*M_PI*pow(r,2);
    numeric  = 0;
    // go over all patches
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      if (!IsItCovering(patch,"NS_OB"))
        continue;

      I  = init_integration();
      I->type = "Integral{f(x)dS},Spectral";
      f  = add_field("int f",0,patch,YES);
      n  = patch->n;
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      // for Z = const
      I->Spectral->Z_surface = 1;
      I->Spectral->K         = n[2]-1; 
      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
           
    // testing outermost0 sphere surface:
    printf(Pretty0"Integral{f(x)dS}|at sphere surface of outermost0 section:\n");
    
    r = Pgetd("grid_outermost_radius");
    analytic = 4*M_PI*pow(r,2);
    numeric  = 0;
    // go over all patches
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      if (!IsItCovering(patch,"outermost_OB"))
        continue;
        
      I  = init_integration();
      I->type = "Integral{f(x)dS},Spectral";
      f  = add_field("int f",0,patch,YES);
      n  = patch->n;
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      // for Z = const
      I->Spectral->Z_surface = 1;
      I->Spectral->K         = n[2]-1; 
      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
           
    // testing outermost0 cube surface:
    printf(Pretty0"Integral{f(x)dS}|at plane surface of outermost0 section:\n");
    
    r = Pgetd("BHNS_separation");
    analytic = 10.*Pow2(r);
    numeric  = 0;
    // go over all patches
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      if (!IsItCovering(patch,"outermost_IB"))
        continue;
      
      I  = init_integration();
      I->type = "Integral{f(x)dS},Spectral";
      f  = add_field("int f",0,patch,YES);
      n  = patch->n;
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      // for Z = const
      I->Spectral->Z_surface = 1;
      I->Spectral->K         = 0; 
      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
  }
  else if (grid->kind == Grid_CubedSpherical_BHNS)
  {
    double r;
    printf("BBN_CubedSpherical_grid:\n");
    
    // testing NS sphere surface:
    printf("\n--> Integral{f(x)dS}|at NS surface section:\n");
    
    r = Pgetd("NS_radius");
    analytic = 4*M_PI*pow(r,2);
    numeric  = 0;
    // go over all patches
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      if (!IsItCovering(patch,"NS"))
        continue;
      if (strstr(patch->name,"left_central_box"))
        continue;

      I  = init_integration();
      I->type = "Integral{f(x)dS},Spectral";
      f  = add_field("int f",0,patch,YES);
      n  = patch->n;
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      // for Z = const
      I->Spectral->Z_surface = 1;
      I->Spectral->K         = n[2]-1; 
      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
           
    // testing outermost0 sphere surface:
    printf("\n--> Integral{f(x)dS}|at sphere surface of outermost0 section:\n");
    
    r = Pgetd("Outermost0_radius");
    analytic = 4*M_PI*pow(r,2);
    numeric  = 0;
    // go over all patches
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      if (!strstr(patch->name,"outermost0"))
        continue;
      
      I  = init_integration();
      I->type = "Integral{f(x)dS},Spectral";
      f  = add_field("int f",0,patch,YES);
      n  = patch->n;
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      // for Z = const
      I->Spectral->Z_surface = 1;
      I->Spectral->K         = n[2]-1; 
      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
           
    // testing outermost0 cube surface:
    printf("\n--> Integral{f(x)dS}|at plane surface of outermost0 section:\n");
    
    r = 2*Pgetd("BHNS_separation");
    analytic = 6*Pow2(r);
    numeric  = 0;
    // go over all patches
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      if (!strstr(patch->name,"outermost0"))
        continue;
      
      I  = init_integration();
      I->type = "Integral{f(x)dS},Spectral";
      f  = add_field("int f",0,patch,YES);
      n  = patch->n;
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      // for Z = const
      I->Spectral->Z_surface = 1;
      I->Spectral->K         = 0; 
      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
  }
  else
    Error0(NO_OPTION);
    
  return EXIT_SUCCESS;
}

/* testing Integral{fdV}, Spectral method.
// ->return value: TEST_SUCCESSFUL */
static int fdV_spectral(Grid_T *const grid)
{
  Integration_T *I;
  Field_T *f;
  Patch_T *patch;
  Uint nn,ijk,p;
  double analytic = 0,numeric = 0;
  
  /* to test this function we use f = 1, so it means we calculate the volume */
  if (grid->kind == Grid_Box)
  {
    analytic = 0;
    numeric = 0;
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      I       = init_integration();
      I->type = "Integral{f(x)dV},Spectral";
      nn      = patch->nn;
      f       = add_field("int f",0,patch,YES);
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      double x1 = patch->min[0];
      double x2 = patch->max[0];
      double y1 = patch->min[1];
      double y2 = patch->max[1];
      double z1 = patch->min[2];
      double z2 = patch->max[2];
      
      /* analytic answer of the integral */
      analytic += (Power(E,y1 + y2 + z1 + z2)*(Power(E,2*x1)*(-2 + Cos(2*x1) + Sin(2*x1)) - 
                   Power(E,2*x2)*(-2 + Cos(2*x2) + Sin(2*x2)))*Sinh(y1 - y2)*Sinh(z1 - z2))/8.;
      
      for (ijk = 0; ijk < nn; ++ijk)
      {
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = Power(E,2*x_(ijk) + 2*y_(ijk) + 2*z_(ijk))*Power(Sin(x_(ijk)),2);
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;
      plan_integration(I);
     
      numeric += execute_integration(I);
     
      free_integration(I);
      remove_field(f);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("Cartesian grid::\n");
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
            numeric,analytic,numeric-analytic);
  }
  else if (grid->kind == Grid_SplitCubedSpherical_BHNS)
  {
    double r;

    printf("BHNS_SplitCubedSpherical_grid:\n");
    
    // testing outermost0: 
    printf(Pretty0"Integral{f(x)dV}|at outermost0 section:\n");
    
    r = Pgetd("grid_outermost_radius");
    analytic = 4./3.*M_PI*pow(r,3)-2.*pow(Pgetd("BHNS_separation"),3);
    numeric  = 0;
    
    // go over all patches 
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      
      if (!IsItCovering(patch,"outermost"))
        continue;
      
      I  = init_integration();
      I->type = "Integral{f(x)dV},Spectral";
      f  = add_field("int f",0,patch,YES);
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
    
    // testing NS:
    printf(Pretty0"Integral{f(x)dV}|at NS section:\n");
    
    r = Pgetd("NS_radius");
    analytic = 4./3.*M_PI*pow(r,3);
    numeric  = 0;
    
    // go over all patches 
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      if (!IsItCovering(patch,"NS"))
        continue;
      
      I  = init_integration();
      I->type = "Integral{f(x)dV},Spectral";
      f  = add_field("int f",0,patch,YES);
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
           
    // testing NS around: 
    printf(Pretty0"Integral{f(x)dV}|at NS around section:\n");
    
    r = Pgetd("NS_radius");
    analytic = pow(Pgetd("BHNS_separation"),3)-4./3.*M_PI*pow(r,3);
    numeric  = 0;
    
    // go over all patches 
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      if (!IsItCovering(patch,"NS_around"))
        continue;
      
      I  = init_integration();
      I->type = "Integral{f(x)dV},Spectral";
      f  = add_field("int f",0,patch,YES);
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
           
  }
  else if (grid->kind == Grid_CubedSpherical_BHNS)
  {
    double r;

    printf("BBN_CubedSpherical_grid:\n");
    
    // testing outermost0: 
    printf("\n--> Integral{f(x)dV}|at outermost0 section:\n");
    
    r = Pgetd("Outermost0_radius");
    analytic = 4./3.*M_PI*pow(r,3)-pow(2*Pgetd("BHNS_separation"),3);
    numeric  = 0;
    
    // go over all patches 
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      if (!strstr(patch->name,"outermost0"))
        continue;
      
      I  = init_integration();
      I->type = "Integral{f(x)dV},Spectral";
      f  = add_field("int f",0,patch,YES);
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
           
    // testing NS: 
    printf("\n--> Integral{f(x)dV}|at NS section:\n");
    
    r = Pgetd("NS_radius");
    analytic = 4./3.*M_PI*pow(r,3);
    numeric  = 0;
    
    // go over all patches 
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      if (!IsItCovering(patch,"NS"))
        continue;
      
      I  = init_integration();
      I->type = "Integral{f(x)dV},Spectral";
      f  = add_field("int f",0,patch,YES);
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
           
    // testing NS around: 
    printf("\n--> Integral{f(x)dV}|at NS around section:\n");
    
    r = Pgetd("NS_radius");
    analytic = pow(Pgetd("BHNS_separation"),3)-4./3.*M_PI*pow(r,3);
    numeric  = 0;
    
    // go over all patches 
    FOR_ALL_PATCHES(p,grid)
    {
      patch   = grid->patch[p];
      if (!strstr(patch->name,"left_NS_around"))
        continue;
      
      I  = init_integration();
      I->type = "Integral{f(x)dV},Spectral";
      f  = add_field("int f",0,patch,YES);
      nn = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);

      for (ijk = 0; ijk < nn; ++ijk)
      { 
        g00[ijk] = g11[ijk] = g22[ijk] = 1.;
        f->v[ijk] = 1.;
      }
        
      I->Spectral->f = f;
      I->g00 = g00;
      I->g01 = g01;
      I->g02 = g02;
      I->g11 = g11;
      I->g12 = g12;
      I->g22 = g22;

      plan_integration(I);
      
      double s0 = execute_integration(I);
      printf("... Integral{f(x)dS}|%s: %e\n",patch->name,s0);
      numeric += s0;
      
      remove_field(f);
      free_integration(I);
      free(g00);
      free(g01);
      free(g02);
      free(g11);
      free(g12);
      free(g22);
    }
    printf("=> numeric = %e, analytic = %e, diff. = %e\n",
           numeric,analytic,numeric-analytic);
           
  }
  else
    Error0(NO_OPTION);
    
  return EXIT_SUCCESS;
}

/* testing Gaussian Quadrature Chebyshev Extrema integration.
// ->return value: TEST_SUCCESSFUL */
static int GQ_ChebExtrema(Grid_T *const grid)
{
  Uint N = 0;
  Integration_T *I = init_integration();
  const char *const par = Pgets("Test_Integration");
  double *f;
  double sf,an;/* resultant */
  double t0,x;
  Uint i;
  
  if (regex_search("[[:digit:]]+",par))
  {
    char *s = regex_find("[[:digit:]]+",par);
    N = (Uint)atoi(s);
    Free(s);
  }
  else
    N = 20000;
 
  f = alloc_double(N);/* integrant */
  t0 = M_PI/(N-1);
  I->type = "Gaussian Quadrature Chebyshev Extrema";
  plan_integration(I);

  /* [-1,1] */
  for (i = 0; i < N; ++i)
  {
    x    = -cos(i*t0);
    f[i] = pow(x,2)+pow(x,4)+10*pow(x,6)+pow(x,3);/* \int f(x)/(1-x^2)dx */
  }
    
  an = 4*M_PI;/* analytic answer from -1 to 1 */

  I->GQ_ChebyshevExtrema->n = N;
  I->GQ_ChebyshevExtrema->f = f;
  sf = execute_integration(I);
  
  if (isfinite(I->err))
    printf("Max expected error for N = %u is %e\n",N,I->err);
  printf("Numeric = %e, Analytic = %e, diff = %e\n",sf,an,fabs(sf-an));

  free(f);
  free_integration(I);
  
  UNUSED(grid);
  return TEST_SUCCESSFUL;
}

/* testing Gaussian Quadrature Lobatto method integration.
// ->return value:  TEST_SUCCESSFUL */
static int GQ_Lobatto(Grid_T *const grid)
{
  Uint N = 0;
  Integration_T *I = init_integration();
  const char *const par = Pgets("Test_Integration");
  double *f;
  double sf,an;/* resultant */
  double x;
  Uint i;
  
  init_Lobatto_root_function();
  
  if (regex_search("[[:digit:]]+",par))
  {
    char *s = regex_find("[[:digit:]]+",par);
    N = (Uint)atoi(s);
    Free(s);
  }
  else
    N = 14;
 
  f = alloc_double(N);/* integrant */
  I->type = "Gaussian Quadrature Lobatto";
  plan_integration(I);

  /* [-1,1] */
  x = -1.; f[0]   = pow(x,2)+pow(x,4)+10*pow(x,6)+pow(x,3);
  x = 1. ; f[N-1] = pow(x,2)+pow(x,4)+10*pow(x,6)+pow(x,3);
  for (i = 1; i <= N-2; ++i)
  {
    x    = Lobatto_root_function(i-1,N-1);
    f[i] = pow(x,2)+pow(x,4)+10*pow(x,6)+pow(x,3);
  }
    
  an = 3.923809523809524;
  I->GQ_Lobatto->n = N;
  I->GQ_Lobatto->f = f;
  sf = execute_integration(I);
  
  printf("Max error for N = %u is %e\n",N,I->err);
  printf("Numeric = %e, Analytic = %e, diff = %e\n ",sf,an,fabs(sf-an));
  
  free(f);
  free_integration(I);
  
  UNUSED(grid);
  return TEST_SUCCESSFUL;
}

/* testing Gaussian Quadrature Legendre method integration.
// ->return value:  TEST_SUCCESSFUL */
static int GQ_Legendre(Grid_T *const grid)
{
  Uint N = 0;
  Integration_T *I = init_integration();
  const char *const par = Pgets("Test_Integration");
  double *f;
  double sf,an;/* resultant */
  double x;
  Uint i;
  
  init_Legendre_root_function();
  
  if (regex_search("[[:digit:]]+",par))
  {
    char *s = regex_find("[[:digit:]]+",par);
    N = (Uint)atoi(s);
    Free(s);
  }
  else
    N = 14;
 
  f = alloc_double(N);/* integrant */
  I->type = "Gaussian Quadrature Legendre";
  plan_integration(I);

  /* [-1,1] */
  for (i = 0; i < N; ++i)
  {
    x    = Legendre_root_function(i,N);
    f[i] = pow(x,2)+pow(x,4)+10*pow(x,6)+pow(x,3);
  }
    
  an = 3.923809523809524;
  I->GQ_Legendre->n = N;
  I->GQ_Legendre->f = f;
  sf = execute_integration(I);
  
  printf("Max error for N = %u is %e\n",N,I->err);
  printf("Numeric = %e, Analytic = %e, diff = %e\n ",sf,an,fabs(sf-an));
  
  free(f);
  free_integration(I);
  
  UNUSED(grid);
  return TEST_SUCCESSFUL;
}

/* testing composite simpson's Rule 1D.
// ->return value: if successful => TEST_SUCCESSFUL, otherwise TEST_UNSUCCESSFUL */
static int csr_1d(Grid_T *const grid)
{
  Uint N = 0;
  const double a = M_PI;
  const double b = 3./2.*M_PI;
  Integration_T *I = init_integration();
  const char *const par = Pgets("Test_Integration");
  double *f;/* integrant */
  double sf,an;/* resultant */
  double dx;
  double err;/* expected error from theory */
  Uint i;
  
  if (regex_search("[[:digit:]]+",par))
  {
    char *s = regex_find("[[:digit:]]+",par);
    N = (Uint)atoi(s);
    Free(s);
  }
  else
    N = 2019;
 
  f = alloc_double(N);
  UNUSED(grid);
  /* testing composite simpson's Rule 1D */
  dx = (b-a)/(N-1);
  for (i = 0; i < N; ++i)
    f[i] = sin(a+i*dx)+pow(a+i*dx,2);/* sin(x)+x^2 */
    
  an = 23.546635705237353;/* from Pi to 3/2*Pi */
  I->type = "Composite Simpson's Rule 1D";
  I->Composite_Simpson_1D->b = b;
  I->Composite_Simpson_1D->a = a;
  I->Composite_Simpson_1D->n = N;
  I->Composite_Simpson_1D->f = f;
  plan_integration(I);
  sf = execute_integration(I);
  err = fabs((b-a)/180*pow(dx,4));
  free(f);
  free_integration(I);
  
  printf("Numeric = %0.15f, Analytic = %0.15f, Expected error = %g\n=> ",sf,an,err);
  
  if (!EQL(fabs(sf-an),err))
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}
