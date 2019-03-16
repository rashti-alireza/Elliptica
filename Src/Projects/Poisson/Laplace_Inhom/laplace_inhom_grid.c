/*
// Alireza Rashti
// June 2018
*/

#include "laplace_inhom_grid.h"

/* making grid for Laplace_Inhom project.
// ->return value: made grid.
*/
Grid_T *Laplace_Inhom_make_grid(void)
{
  Grid_T *grid = alloc_grid();/* adding a new grid */
  grid_characteristics_Laplace_Inhome(grid);
  make_patches(grid);/* making patch(es) to cover the grid */
  realize_geometry(grid);/* realizing the geometry of whole grid
                     // including the way patches have been sewed,
                     // normal to the boundary, 
                     // outer-boundary, inner boundary and etc. */
  return grid;
}

/* calculating the main characteristics of grid needed for setting up
// the patches, like lengths, radii and etc */
static void grid_characteristics_Laplace_Inhome (Grid_T *const grid)
{
  const char *kind;
  
  /* finding the kind of grid */
  kind = GetParameterS_E("grid_kind");
  grid->kind = dup_s(kind);
  
  if (strcmp_i(grid->kind,"Cartesian_grid"))
    characteristics_Cartesian_grid(grid);
  
  else if (strcmp_i(grid->kind,"BNS_Projective_grid"))
    characteristics_BNS_Projective_grid(grid); 
    
  else
    abortEr_s("There is no such %s grid kind.\n",grid->kind);
  
}

/* calculating the main characteristic of grid for Cartesian grid */
static void characteristics_Cartesian_grid(Grid_T *const grid)
{
  /* this type of grid is so simple; nothing to calculate. */
  UNUSED(grid);
}

/* calculating the main characteristic of grid for BNS_Projective grid */
static void characteristics_BNS_Projective_grid(Grid_T *const grid)
{
  const unsigned gn   = grid->gn;
  const double CONST  = 5.0;
  const double C      = GetParameterD_E("Centers_Distance");
  const double R_NS_l = GetParameterD_E("NS_left_radius");/* assuming perfect sphere */
  const double R_NS_r = GetParameterD_E("NS_right_radius");/* assuming perfect sphere */
  const unsigned N_Outermost_Split = (unsigned)GetParameterI_E("Number_of_Outermost_Split"); 
  double O,O_l,O_r,
         R_Surr_l,R_Surr_r,
         R_inside_l,R_inside_r,
         box_size_l,box_size_r,
         *R_outmost_l = alloc_double(N_Outermost_Split),
         *R_outmost_r = alloc_double(N_Outermost_Split),
         *R0 = alloc_double(N_Outermost_Split);
  unsigned n_box_l, n_box_r,n_r,n_l;
  double m,M,s;
  char var[100] = {'\0'};
  char par[100] = {'\0'};
  char val[100] = {'\0'};
  unsigned i;
  
  /* making field of NS's radius */
  make_field_of_NS_radius(grid);
  
  n_l = (unsigned)GetParameterI("n_c");
  i   = (unsigned)GetParameterI("NS_left_n_c");
  if (i != INT_MAX) n_l = i;
  if (n_l == INT_MAX)
    abortEr("n_l could not be set.");
    
  n_r = (unsigned)GetParameterI("n_c");
  i   = (unsigned)GetParameterI("NS_right_n_c");
  if (i != INT_MAX) n_r = i;
  if (n_r == INT_MAX)
    abortEr("n_r could not be set.");
  
  O = 2*(C-R_NS_l-R_NS_r);
  O_l = O/2+R_NS_l;
  O_r = O/2+R_NS_r;
  
  R_inside_l = R_NS_l/n_l;
  R_inside_r = R_NS_r/n_r;
  
  box_size_l = 2*R_NS_l/(n_l-1);
  box_size_r = 2*R_NS_r/(n_r-1);
  
  m = O_l > O_r ? O_r : O_l;
  M = O_l > O_r ? O_l : O_r;
  s = SQR(m/M);
  R_Surr_l = sqrt(SQR(O_l)+s);
  R_Surr_r = sqrt(SQR(O_r)+s);
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R0[i] = GetParameterD_E(var);
    
    if (i > 0)
      if (LSS(R0[i],R0[i-1]))
        abortEr("The radius of outermost must be increasing.");
        
    R_outmost_l[i] = sqrt(SQR(O_l)+SQR(R0[i]));
    R_outmost_r[i] = sqrt(SQR(O_r)+SQR(R0[i]));
  }
  
  n_box_l = (unsigned)(3 + n_l/CONST);
  n_box_r = (unsigned)(3 + n_r/CONST);
  
  /* adding the results to the parameter data base */
  
  /* n_a, n_b, n_c */
  sprintf(par,"grid%u_centeral_box_left_n_abc",gn);
  sprintf(val,"%u",n_box_l);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_centeral_box_right_n_abc",gn);
  sprintf(val,"%u",n_box_r);
  add_parameter_string(par,val);
  
  /* size a,b,c */
  sprintf(par,"grid%u_centeral_box_left_size_a",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_centeral_box_left_size_b",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_centeral_box_left_size_c",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_centeral_box_right_size_a",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_centeral_box_right_size_b",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_centeral_box_right_size_c",gn);
  add_parameter_double(par,box_size_r);
  
  /* R inside */
  sprintf(par,"grid%u_NS_R_inside_left",gn);
  add_parameter_double(par,R_inside_l);
  
  sprintf(par,"grid%u_NS_R_inside_right",gn);
  add_parameter_double(par,R_inside_r);
  
  /* R2 surroundings */
  sprintf(par,"grid%u_NS_Surrounding_R2_left",gn);
  add_parameter_double(par,R_Surr_l);
  
  sprintf(par,"grid%u_NS_Surrounding_R2_right",gn);
  add_parameter_double(par,R_Surr_r);
  
  /* R1 and R2 outermost */
  for (i = 0; i < N_Outermost_Split; i++)
  {
    /* R1: */
    if (i == 0)
    {
      sprintf(par,"grid%u_Outermost%u_R1_left",gn,i);
      add_parameter_double(par,R_Surr_l);
      
      sprintf(par,"grid%u_Outermost%u_R1_right",gn,i);
      add_parameter_double(par,R_Surr_r);
    }
    else
    {
      sprintf(par,"grid%u_Outermost%u_R1_left",gn,i);
      add_parameter_double(par,R_outmost_l[i-1]);
      
      sprintf(par,"grid%u_Outermost%u_R1_right",gn,i);
      add_parameter_double(par,R_outmost_r[i-1]);
    }
    
    /* R2: */
    sprintf(par,"grid%u_Outermost%u_R2_left",gn,i);
    add_parameter_double(par,R_outmost_l[i]);
    
    sprintf(par,"grid%u_Outermost%u_R2_right",gn,i);
    add_parameter_double(par,R_outmost_r[i]);
  }
  
  /* assuming the center of left NS at (0,-O_l,0) */
  sprintf(par,"grid%u_NS_left_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_NS_left_center_b",gn);
  add_parameter_double(par,-O_l);
  
  sprintf(par,"grid%u_NS_left_center_c",gn);
  add_parameter_double(par,0.0);
  
  /* assuming the center of right NS at (0,O_r,0) */
  sprintf(par,"grid%u_NS_right_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_NS_right_center_b",gn);
  add_parameter_double(par,O_r);
  
  sprintf(par,"grid%u_NS_right_center_c",gn);
  add_parameter_double(par,0.0);
  
  free(R0);
  free(R_outmost_r);
  free(R_outmost_l);
}

/* making field of NS's radius */
static void make_field_of_NS_radius(Grid_T *const grid)
{
  const double R_NS_l = GetParameterD_E("NS_left_radius");/* assuming perfect sphere */
  const double R_NS_r = GetParameterD_E("NS_right_radius");/* assuming perfect sphere */
  double *R_field;
  char var[100] = {'\0'};
  unsigned N[3],n,i,j,N_total;
  
  /* left side: */
  
  /* filling n */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("NS_left_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("NS_left_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("NS_left_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  R_field = alloc_double(N_total);
  
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R_field[L(N,i,j,0)] = R_NS_l;
      
  sprintf(var,"grid%u_NS_left_radius",grid->gn);
  add_parameter_array(var,R_field,N_total);
  free(R_field);
  
  /* right side: */
  
  /* filling n */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("NS_right_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("NS_right_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("NS_right_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  R_field = alloc_double(N_total);
  
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R_field[L(N,i,j,0)] = R_NS_r;
  
  sprintf(var,"grid%u_NS_left_radius",grid->gn);
  add_parameter_array(var,R_field,N_total);
  free(R_field);
}

