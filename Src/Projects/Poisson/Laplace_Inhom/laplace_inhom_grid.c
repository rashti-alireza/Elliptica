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
  const double R_NS_l = GetParameterD_E("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = GetParameterD_E("right_NS_radius");/* assuming perfect sphere */
  double R_max_l,R_max_r;/* maximum distance from the center of each star */
  const unsigned N_Outermost_Split = (unsigned)GetParameterI_E("Number_of_Outermost_Split"); 
  double O,O_l,O_r,
         R_Surr_l,R_Surr_r,
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
  
  assert(GRT(C,0));
  assert(GRT(R_NS_l,0));
  assert(GRT(R_NS_r,0));
  
  n_l = (unsigned)GetParameterI("n_c");
  i   = (unsigned)GetParameterI("left_NS_n_c");
  if (i != INT_MAX) 	n_l = i;
  if (n_l == INT_MAX)   abortEr("n_l could not be set.");
  assert(n_l > 2);
  
  n_r = (unsigned)GetParameterI("n_c");
  i   = (unsigned)GetParameterI("right_NS_n_c");
  if (i != INT_MAX) 	n_r = i;
  if (n_r == INT_MAX)   abortEr("n_r could not be set.");
  assert(n_r > 2);
  
  /* making field of NS's radius and 
  // finding the max distance from the center of the star */
  make_field_of_NS_radii(grid,&R_max_l,&R_max_r);
  
  box_size_l = 2*R_max_l/(n_l-1);
  box_size_r = 2*R_max_r/(n_r-1);
  
  O = C-R_NS_l-R_NS_r;
  assert(GRT(O,0));
  O_l = O/2+R_NS_l;
  O_r = O/2+R_NS_r;
  
  M = GRT(O_l,O_r) ? O_l : O_r;
  m = LSS(R_max_l,R_max_r) ? R_max_l : R_max_r;
  s = SQR(m/M);
  R_Surr_l = sqrt(SQR(O_l)+s);
  R_Surr_r = sqrt(SQR(O_r)+s);
  assert(LSS(R_Surr_l-O_l,O/2));
  assert(LSS(R_Surr_r-O_r,O/2));
  
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
  sprintf(par,"grid%u_left_centeral_box_n_abc",gn);
  sprintf(val,"%u",n_box_l);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_right_centeral_box_n_abc",gn);
  sprintf(val,"%u",n_box_r);
  add_parameter_string(par,val);
  
  /* size a,b,c */
  sprintf(par,"grid%u_left_centeral_box_size_a",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_b",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_c",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_right_centeral_box_size_a",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_right_centeral_box_size_b",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_right_centeral_box_size_c",gn);
  add_parameter_double(par,box_size_r);
  
  /* R2 surroundings */
  sprintf(par,"grid%u_left_NS_surrounding_R2",gn);
  add_parameter_double(par,R_Surr_l);
  
  sprintf(par,"grid%u_right_NS_surrounding_R2",gn);
  add_parameter_double(par,R_Surr_r);
  
  /* R1 and R2 outermost */
  for (i = 0; i < N_Outermost_Split; i++)
  {
    /* R1: */
    if (i == 0)
    {
      sprintf(par,"grid%u_left_outermost%u_R1",gn,i);
      add_parameter_double(par,R_Surr_l);
      
      sprintf(par,"grid%u_right_outermost%u_R1",gn,i);
      add_parameter_double(par,R_Surr_r);
    }
    else
    {
      sprintf(par,"grid%u_left_outermost%u_R1",gn,i);
      add_parameter_double(par,R_outmost_l[i-1]);
      
      sprintf(par,"grid%u_right_outermost%u_R1",gn,i);
      add_parameter_double(par,R_outmost_r[i-1]);
    }
    
    /* R2: */
    sprintf(par,"grid%u_left_outermost%u_R2",gn,i);
    add_parameter_double(par,R_outmost_l[i]);
    
    sprintf(par,"grid%u_right_outermost%u_R2",gn,i);
    add_parameter_double(par,R_outmost_r[i]);
  }
  
  /* assuming the center of left NS at (0,-O_l,0) */
  sprintf(par,"grid%u_left_NS_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_left_NS_center_b",gn);
  add_parameter_double(par,-O_l);
  
  sprintf(par,"grid%u_left_NS_center_c",gn);
  add_parameter_double(par,0.0);
  
  /* assuming the center of right NS at (0,O_r,0) */
  sprintf(par,"grid%u_right_NS_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_right_NS_center_b",gn);
  add_parameter_double(par,O_r);
  
  sprintf(par,"grid%u_right_NS_center_c",gn);
  add_parameter_double(par,0.0);
  
  free(R0);
  free(R_outmost_r);
  free(R_outmost_l);
}

/* making field of NS's inside (R1) and surface (R2) radius and 
// finding the max distance from the center of the star */
static void make_field_of_NS_radii(Grid_T *const grid,double *const R_max_l,double *const R_max_r)
{
  const double R_NS_l = GetParameterD_E("left_NS_radius");/* assuming perfect sphere */
  const double R_NS_r = GetParameterD_E("right_NS_radius");/* assuming perfect sphere */
  double *R1_down,*R2_down;
  double *R1_up,*R2_up;
  char par[100] = {'\0'};
  unsigned N[3],n,i,j,N_total;
  
  /* left NS */
  *R_max_l = R_NS_l;/* since this star is perfect sphere */
  
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("left_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("left_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("left_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface: */
  
  /* up side */
  R2_up = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R2_up[L(N,i,j,0)] = R_NS_l;
      
  sprintf(par,"grid%u_left_NS_R2_up",grid->gn);
  add_parameter_array(par,R2_up,N_total);
  /* end of up side */
  /* down side */
  R2_down = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R2_down[L(N,i,j,0)] = R_NS_l;
      
  sprintf(par,"grid%u_left_NS_R2_down",grid->gn);
  add_parameter_array(par,R2_down,N_total);
  /* end of down side */
  
  /* inside: */
  
  /* up side */
  R1_up = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R1_up[L(N,i,j,0)] = R2_up[L(N,i,j,0)]/N[2];
      
  sprintf(par,"grid%u_left_NS_R1_up",grid->gn);
  add_parameter_array(par,R1_up,N_total);
  free(R1_up);
  free(R2_up);
  /* end of up side */
  /* down side */
  R1_down = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R1_down[L(N,i,j,0)] = R2_down[L(N,i,j,0)]/N[2];
      
  sprintf(par,"grid%u_left_NS_R1_down",grid->gn);
  add_parameter_array(par,R1_down,N_total);
  free(R1_down);
  free(R2_down);
  /* end of down side */
  
  /* right NS */
  
  *R_max_r = R_NS_r;/* since this star is perfect sphere */
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("right_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("right_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("right_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface: */
  
  /* up side */
  R2_up = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R2_up[L(N,i,j,0)] = R_NS_r;
      
  sprintf(par,"grid%u_right_NS_R2_up",grid->gn);
  add_parameter_array(par,R2_up,N_total);
  /* end of up side */
  /* down side */
  R2_down = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R2_down[L(N,i,j,0)] = R_NS_r;
      
  sprintf(par,"grid%u_right_NS_R2_down",grid->gn);
  add_parameter_array(par,R2_down,N_total);
  /* end of down side */
  
  /* inside: */
  
  /* up side */
  R1_up = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R1_up[L(N,i,j,0)] = R2_up[L(N,i,j,0)]/N[2];
      
  sprintf(par,"grid%u_right_NS_R1_up",grid->gn);
  add_parameter_array(par,R1_up,N_total);
  free(R1_up);
  free(R2_up);
  /* end of up side */
  /* down side */
  R1_down = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      R1_down[L(N,i,j,0)] = R2_down[L(N,i,j,0)]/N[2];
      
  sprintf(par,"grid%u_right_NS_R1_down",grid->gn);
  add_parameter_array(par,R1_down,N_total);
  free(R1_down);
  free(R2_down);
  /* end of down side */
  
}

