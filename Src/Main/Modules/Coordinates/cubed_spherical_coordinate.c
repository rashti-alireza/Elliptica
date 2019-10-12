/*
// Alireza Rashti
// April 2019
*/

#include "cubed_spherical_coordinate.h"

/* filling cubed spherical + box coordinate patches for SNS grid */
void fill_patches_SNS_CubedSpherical_Box_grid(Grid_T *const grid)
{
  const unsigned N_outermost_split = (unsigned) GetParameterI_E("Number_of_Outermost_Split");
  unsigned i,pn;
  
  pn = 0; /* patch number */
  populate_left_NS_central_box(grid,pn++);/* +1 */
  populate_left_NS(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_left_NS_surrounding(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_right_box_sns(grid,pn);
  pn += 1; /* +1 box */
  populate_filling_box_CubedSpherical(grid,pn++,UP);
  populate_filling_box_CubedSpherical(grid,pn++,DOWN);
  populate_filling_box_CubedSpherical(grid,pn++,BACK);
  populate_filling_box_CubedSpherical(grid,pn++,FRONT);
  
  for (i = 0; i < N_outermost_split; i++)
  {
    populate_outermost(grid,pn,i);
    pn += 6; /* +6 cubed sphere */
  }

}


/* filling cubed spherical coordinate patches for BNS grid */
void fill_patches_BNS_CubedSpherical_grid(Grid_T *const grid)
{
  const unsigned N_outermost_split = (unsigned) GetParameterI_E("Number_of_Outermost_Split");
  unsigned i,pn;
  
  pn = 0; /* patch number */
  populate_left_NS_central_box(grid,pn++);/* +1 */
  populate_left_NS(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_left_NS_surrounding(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_right_NS_central_box(grid,pn++);/*+1 */
  populate_right_NS(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_right_NS_surrounding(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_filling_box_CubedSpherical(grid,pn++,UP);
  populate_filling_box_CubedSpherical(grid,pn++,DOWN);
  populate_filling_box_CubedSpherical(grid,pn++,BACK);
  populate_filling_box_CubedSpherical(grid,pn++,FRONT);
  
  for (i = 0; i < N_outermost_split; i++)
  {
    populate_outermost(grid,pn,i);
    pn += 6; /* +6 cubed sphere */
  }

}

/* filling cubed spherical coordinate patches for BBN grid */
void fill_patches_BBN_CubedSpherical_grid(Grid_T *const grid)
{
  const unsigned N_outermost_split = (unsigned) GetParameterI_E("Number_of_Outermost_Split");
  unsigned i,pn;
  
  pn = 0; /* patch number */
  populate_left_NS_central_box(grid,pn++);/* +1 */
  populate_left_NS(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_left_NS_surrounding(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_right_BH_surrounding(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_filling_box_CubedSpherical(grid,pn++,UP);
  populate_filling_box_CubedSpherical(grid,pn++,DOWN);
  populate_filling_box_CubedSpherical(grid,pn++,BACK);
  populate_filling_box_CubedSpherical(grid,pn++,FRONT);
  
  for (i = 0; i < N_outermost_split; i++)
  {
    populate_outermost(grid,pn,i);
    pn += 6; /* +6 cubed sphere */
  }

}

/* making value of coords. it is a general function for cubed spherical type */
void make_nodes_CubedSpherical_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
  const Flag_T type = patch->CoordSysInfo->CubedSphericalCoord->type;
  double S; /* sign */
  unsigned a,b,c;/* permuted indices */
  const unsigned nn = patch->nn;
  const unsigned *const n = patch->n;
  const Field_T *const R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f,
                *const R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
  const double xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1,
               xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2,
                R1 = patch->CoordSysInfo->CubedSphericalCoord->R1,
                R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
               
  const double *const C = patch->c;/* center of origine translated */
  unsigned i,j,k,l;
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  initialize_collocation_struct(patch,&coll_s[2],2);

  /* setting up sign and indices for cubde sphere based on side*/  
  SignAndIndex_permutation_CubedSphere(side,&a,&b,&c,&S);
  
  switch (type)
  {
    case NS_T_CS:
      for (l = 0; l < nn; l++)
      {
        double *X = alloc_double(3);
        double *x = patch->node[l]->x;
        double x1,x2,d;
        
        IJK(l,n,&i,&j,&k);
        X[0] = point_value(i,&coll_s[0]);
        X[1] = point_value(j,&coll_s[1]);
        X[2] = point_value(k,&coll_s[2]);
        d = sqrt(1+SQR(X[0])+SQR(X[1]));
        patch->node[l]->X = X;
        x1 = xc1;
        x2 = S*R2_f->v[L(n,i,j,0)]/d;
        
        x[c] = x1+(x2-x1)*X[2];
        x[a] = S*x[c]*X[0];
        x[b] = S*x[c]*X[1];
        
        x[a]+= C[a];
        x[b]+= C[b];
        x[c]+= C[c];
      }
    break;
    case SR_T_CS:
      for (l = 0; l < nn; l++)
      {
        double *X = alloc_double(3);
        double *x = patch->node[l]->x;
        double x1,x2,d;
        
        IJK(l,n,&i,&j,&k);
        X[0] = point_value(i,&coll_s[0]);
        X[1] = point_value(j,&coll_s[1]);
        X[2] = point_value(k,&coll_s[2]);
        d = sqrt(1+SQR(X[0])+SQR(X[1]));
        patch->node[l]->X = X;
        x2 = xc2;
        x1 = S*R1_f->v[L(n,i,j,0)]/d;
        
        x[c] = x1+(x2-x1)*X[2];
        x[a] = S*x[c]*X[0];
        x[b] = S*x[c]*X[1];
        
        x[a]+= C[a];
        x[b]+= C[b];
        x[c]+= C[c];
      }
    break;
    case OT_T1_CS:
      for (l = 0; l < nn; l++)
      {
        double *X = alloc_double(3);
        double *x = patch->node[l]->x;
        double x1,x2,d;
        
        IJK(l,n,&i,&j,&k);
        X[0] = point_value(i,&coll_s[0]);
        X[1] = point_value(j,&coll_s[1]);
        X[2] = point_value(k,&coll_s[2]);
        d = sqrt(1+SQR(X[0])+SQR(X[1]));
        patch->node[l]->X = X;
        x1 = xc1;
        x2 = S*R2/d;
        
        x[c] = x1+(x2-x1)*X[2];
        x[a] = S*x[c]*X[0];
        x[b] = S*x[c]*X[1];
        
        x[a]+= C[a];
        x[b]+= C[b];
        x[c]+= C[c];
      }
    break;
    case OT_T2_CS:
      for (l = 0; l < nn; l++)
      {
        double *X = alloc_double(3);
        double *x = patch->node[l]->x;
        double x1,x2,d;
        
        IJK(l,n,&i,&j,&k);
        X[0] = point_value(i,&coll_s[0]);
        X[1] = point_value(j,&coll_s[1]);
        X[2] = point_value(k,&coll_s[2]);
        d = sqrt(1+SQR(X[0])+SQR(X[1]));
        patch->node[l]->X = X;
        x1 = S*R1/d;
        x2 = S*R2/d;
        
        x[c] = x1+(x2-x1)*X[2];
        x[a] = S*x[c]*X[0];
        x[b] = S*x[c]*X[1];
        
        x[a]+= C[a];
        x[b]+= C[b];
        x[c]+= C[c];
      }
    break;
    default:
      abortEr(NO_OPTION);
  }
  
}

/* making Jacobian transformation for cubed spherical  coord. */
void make_JacobianT_CubedSpherical_coord(Patch_T *const patch)
{
  const Flag_T type = patch->CoordSysInfo->CubedSphericalCoord->type;
  const Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
  
  if (type == NS_T_CS)
  {
    switch (side)
    {
      case UP:
        patch->JacobianT->j      = JT_NS_T_CS_up;
      break;
      case DOWN:
        patch->JacobianT->j      = JT_NS_T_CS_down;
      break;
      case LEFT:
        patch->JacobianT->j      = JT_NS_T_CS_left;
      break;
      case RIGHT:
        patch->JacobianT->j      = JT_NS_T_CS_right;
      break;
      case BACK:
        patch->JacobianT->j      = JT_NS_T_CS_back;
      break;
      case FRONT:
        patch->JacobianT->j      = JT_NS_T_CS_front;
      break;
      default:
        abortEr(NO_JOB);
    }/* end of switch */
    R2_derivative(patch);/* surface function derivative */
  }/* end of if (type == NS_T_CS) */
  
  else if (type == SR_T_CS)
  {
    switch (side)
    {
      case UP:
        patch->JacobianT->j      = JT_SR_T_CS_up;
      break;
      case DOWN:
        patch->JacobianT->j      = JT_SR_T_CS_down;
      break;
      case LEFT:
        patch->JacobianT->j      = JT_SR_T_CS_left;
      break;
      case RIGHT:
        patch->JacobianT->j      = JT_SR_T_CS_right;
      break;
      case BACK:
        patch->JacobianT->j      = JT_SR_T_CS_back;
      break;
      case FRONT:
        patch->JacobianT->j      = JT_SR_T_CS_front;
      break;
      default:
        abortEr(NO_JOB);
    }/* end of switch */
    R1_derivative(patch);/* surface function derivative */
  }/* end of else if (type == SR_T_CS) */
  
  else if (type == OT_T1_CS)
  {
    switch (side)
    {
      case UP:
        patch->JacobianT->j      = JT_OT_T1_CS_up;
      break;
      case DOWN:
        patch->JacobianT->j      = JT_OT_T1_CS_down;
      break;
      case LEFT:
        patch->JacobianT->j      = JT_OT_T1_CS_left;
      break;
      case RIGHT:
        patch->JacobianT->j      = JT_OT_T1_CS_right;
      break;
      case BACK:
        patch->JacobianT->j      = JT_OT_T1_CS_back;
      break;
      case FRONT:
        patch->JacobianT->j      = JT_OT_T1_CS_front;
      break;
      default:
        abortEr(NO_JOB);
    }/* end of switch */
  }/* end of else if (type == OT_T1_CS) */
  else if (type == OT_T2_CS)
  {
    switch (side)
    {
      case UP:
        patch->JacobianT->j      = JT_OT_T2_CS_up;
      break;
      case DOWN:
        patch->JacobianT->j      = JT_OT_T2_CS_down;
      break;
      case LEFT:
        patch->JacobianT->j      = JT_OT_T2_CS_left;
      break;
      case RIGHT:
        patch->JacobianT->j      = JT_OT_T2_CS_right;
      break;
      case BACK:
        patch->JacobianT->j      = JT_OT_T2_CS_back;
      break;
      case FRONT:
        patch->JacobianT->j      = JT_OT_T2_CS_front;
      break;
      default:
        abortEr(NO_JOB);
    }/* end of switch */
  }/* end of else if (type == OT_T2_CS) */
  
}

/* preparing R1 derivatives of Cubed Spherical coords.
// NOTE: one must remove dR_? after each updating of R. */
static void R1_derivative(Patch_T *const patch)
{
  Field_T *dR1_dX = add_field("dR1_dX",0,patch,NO),
          *dR1_dY = add_field("dR1_dY",0,patch,NO),
          *dR1_dx = add_field("dR1_dx",0,patch,YES),
          *dR1_dy = add_field("dR1_dy",0,patch,YES),
          *dR1_dz = add_field("dR1_dz",0,patch,YES);
  Field_T *const R1 = patch->pool[Ind("surface_function")];
  const unsigned nn = patch->nn;
  unsigned ijk;
          
  dR1_dX->v = Partial_Derivative(R1,"a");
  dR1_dY->v = Partial_Derivative(R1,"b");
    
  for (ijk = 0; ijk < nn; ++ijk)
  {
    dR1_dx->v[ijk] = dR1_dX->v[ijk]*dq2_dq1(patch,_a_,_x_,ijk)+
                     dR1_dY->v[ijk]*dq2_dq1(patch,_b_,_x_,ijk);
    dR1_dy->v[ijk] = dR1_dX->v[ijk]*dq2_dq1(patch,_a_,_y_,ijk)+
                     dR1_dY->v[ijk]*dq2_dq1(patch,_b_,_y_,ijk);
    dR1_dz->v[ijk] = dR1_dX->v[ijk]*dq2_dq1(patch,_a_,_z_,ijk)+
                     dR1_dY->v[ijk]*dq2_dq1(patch,_b_,_z_,ijk);
  }
                      
  remove_field(dR1_dX);
  remove_field(dR1_dY);
  
  patch->CoordSysInfo->CubedSphericalCoord->dR1_dx = dR1_dx;
  patch->CoordSysInfo->CubedSphericalCoord->dR1_dy = dR1_dy;
  patch->CoordSysInfo->CubedSphericalCoord->dR1_dz = dR1_dz;
}

/* preparing R2 derivatives of Cubed Spherical coords.
// NOTE: one must remove dR_? after each updating of R. */
static void R2_derivative(Patch_T *const patch)
{
  Field_T *dR2_dX = add_field("dR2_dX",0,patch,NO),
          *dR2_dY = add_field("dR2_dY",0,patch,NO),
          *dR2_dx = add_field("dR2_dx",0,patch,YES),
          *dR2_dy = add_field("dR2_dy",0,patch,YES),
          *dR2_dz = add_field("dR2_dz",0,patch,YES);
  Field_T *const R2 = patch->pool[Ind("surface_function")];
  const unsigned nn = patch->nn;
  unsigned ijk;
          
  dR2_dX->v = Partial_Derivative(R2,"a");
  dR2_dY->v = Partial_Derivative(R2,"b");
    
  for (ijk = 0; ijk < nn; ++ijk)
  {
    dR2_dx->v[ijk] = dR2_dX->v[ijk]*dq2_dq1(patch,_a_,_x_,ijk)+
                     dR2_dY->v[ijk]*dq2_dq1(patch,_b_,_x_,ijk);
    dR2_dy->v[ijk] = dR2_dX->v[ijk]*dq2_dq1(patch,_a_,_y_,ijk)+
                     dR2_dY->v[ijk]*dq2_dq1(patch,_b_,_y_,ijk);
    dR2_dz->v[ijk] = dR2_dX->v[ijk]*dq2_dq1(patch,_a_,_z_,ijk)+
                     dR2_dY->v[ijk]*dq2_dq1(patch,_b_,_z_,ijk);
  }
                      
  remove_field(dR2_dX);
  remove_field(dR2_dY);
  
  patch->CoordSysInfo->CubedSphericalCoord->dR2_dx = dR2_dx;
  patch->CoordSysInfo->CubedSphericalCoord->dR2_dy = dR2_dy;
  patch->CoordSysInfo->CubedSphericalCoord->dR2_dz = dR2_dz;
}

/* interpolation of 2d field of radius R (surface function)
// for Cubed Spherical coordinate.
// R is only function of X[0] and X[1].
// Note: R must be populated like a 3-d field f(X,Y,Z)
// but its value is equal on all slices of Z; thus, df/dZ = 0.
// X = the curvilinear coord in which we want R(X).
// ->return value: R(X) */
double R_interpolation_CS(Field_T *const R,const double *const X)
{
  double interp;
  Interpolation_T *interp_s = init_interpolation();
  
  interp_s->field = R;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->XY_dir_flag = 1;
  interp_s->K = 0;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return interp;
}

/* populating properties of patch for left NS */
static void populate_left_NS(Grid_T *const grid,const unsigned pn)
{
  unsigned p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R2 = add_field("surface_function",0,patch,NO);
    double *R2_array;
    Flag_T side = p-pn;
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    unsigned n,ijk;
    
    /* filling flags */
    patch->CoordSysInfo->CubedSphericalCoord->side = side;
    patch->CoordSysInfo->CubedSphericalCoord->type = NS_T_CS;
    
    /* filling grid */
    patch->grid = grid;
    
    /* filling patch number */
    patch->pn = p;
    
    /* filling inner boundary */
    patch->innerB = 0;
    
    /* filling n */
    patch->n[0] = (unsigned)GetParameterI("n_a");
    patch->n[1] = (unsigned)GetParameterI("n_b");
    patch->n[2] = (unsigned)GetParameterI("n_c");
    /* check for override */
    sprintf(var,"left_NS");
    make_keyword_parameter(&ret,var,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      abortEr("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      abortEr("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      abortEr("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,"grid%u_left_NS_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_up",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_left_centeral_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,"grid%u_left_NS_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_down",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_left_centeral_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,"grid%u_left_NS_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_left",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_left_centeral_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];

      break;
      case RIGHT:
        /* filling name */
        sprintf(name,"grid%u_left_NS_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_right",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_left_centeral_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,"grid%u_left_NS_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_back",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_left_centeral_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,"grid%u_left_NS_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_front",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_left_centeral_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      default:
        abortEr(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,"grid%u_left_NS_center_a",grid->gn);
    patch->c[0] = GetParameterDoubleF_E(var);
    sprintf(var,"grid%u_left_NS_center_b",grid->gn);
    patch->c[1] = GetParameterDoubleF_E(var);
    sprintf(var,"grid%u_left_NS_center_c",grid->gn);
    patch->c[2] = GetParameterDoubleF_E(var);
    
    /* filling min */
    patch->min[0] = -1;
    patch->min[1] = -1;
    patch->min[2] = 0;
    
    /* filling max */
    patch->max[0] = 1;
    patch->max[1] = 1;
    patch->max[2] = 1;
    
    /* filling flags */
    patch->coordsys = CubedSpherical;
    
    /* collocation */
    patch->collocation[0] = Chebyshev_Extrema;
    patch->collocation[1] = Chebyshev_Extrema;
    patch->collocation[2] = Chebyshev_Extrema;
    
    /* basis */
    patch->basis[0] = Chebyshev_Tn_BASIS;
    patch->basis[1] = Chebyshev_Tn_BASIS;
    patch->basis[2] = Chebyshev_Tn_BASIS;
  }
}

/* populating properties of patch for right NS */
static void populate_right_NS(Grid_T *const grid,const unsigned pn)
{
  unsigned p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R2 = add_field("surface_function",0,patch,NO);
    double *R2_array;
    Flag_T side = p-pn;
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    unsigned n,ijk;
    
    /* filling flags */
    patch->CoordSysInfo->CubedSphericalCoord->side = side;
    patch->CoordSysInfo->CubedSphericalCoord->type = NS_T_CS;
        
    /* filling grid */
    patch->grid = grid;
    
    /* filling patch number */
    patch->pn = p;
    
    /* filling inner boundary */
    patch->innerB = 0;
    
    /* filling n */
    patch->n[0] = (unsigned)GetParameterI("n_a");
    patch->n[1] = (unsigned)GetParameterI("n_b");
    patch->n[2] = (unsigned)GetParameterI("n_c");
    /* check for override */
    sprintf(var,"right_NS");
    make_keyword_parameter(&ret,var,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      abortEr("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      abortEr("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      abortEr("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);

    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,"grid%u_right_NS_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_up",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_right_centeral_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,"grid%u_right_NS_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_down",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_right_centeral_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,"grid%u_right_NS_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_left",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_right_centeral_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,"grid%u_right_NS_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_right",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_right_centeral_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,"grid%u_right_NS_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_back",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_right_centeral_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,"grid%u_right_NS_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_front",grid->gn);
        R2_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,"grid%u_right_centeral_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = GetParameterDoubleF_E(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      default:
        abortEr(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,"grid%u_right_NS_center_a",grid->gn);
    patch->c[0] = GetParameterDoubleF_E(var);
    sprintf(var,"grid%u_right_NS_center_b",grid->gn);
    patch->c[1] = GetParameterDoubleF_E(var);
    sprintf(var,"grid%u_right_NS_center_c",grid->gn);
    patch->c[2] = GetParameterDoubleF_E(var);
    
    /* filling min */
    patch->min[0] = -1;
    patch->min[1] = -1;
    patch->min[2] = 0;
    
    /* filling max */
    patch->max[0] = 1;
    patch->max[1] = 1;
    patch->max[2] = 1;
    
    /* filling flags */
    patch->coordsys = CubedSpherical;
    
    /* collocation */
    patch->collocation[0] = Chebyshev_Extrema;
    patch->collocation[1] = Chebyshev_Extrema;
    patch->collocation[2] = Chebyshev_Extrema;
    
    /* basis */
    patch->basis[0] = Chebyshev_Tn_BASIS;
    patch->basis[1] = Chebyshev_Tn_BASIS;
    patch->basis[2] = Chebyshev_Tn_BASIS;
  }
}

/* populating properties of patch for right NS surrounding */
static void populate_right_NS_surrounding(Grid_T *const grid,const unsigned pn)
{
  unsigned p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R1 = add_field("surface_function",0,patch,NO);
    double *R1_array;
    Flag_T side = p-pn;
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    unsigned n,ijk;
    
    /* filling flags */
    patch->CoordSysInfo->CubedSphericalCoord->side = side;
    patch->CoordSysInfo->CubedSphericalCoord->type = SR_T_CS;
    
    /* filling grid */
    patch->grid = grid;
    
    /* filling patch number */
    patch->pn = p;
    
    /* filling inner boundary */
    patch->innerB = 0;
    
    /* filling n */
    patch->n[0] = (unsigned)GetParameterI("n_a");
    patch->n[1] = (unsigned)GetParameterI("n_b");
    patch->n[2] = (unsigned)GetParameterI("n_c");
    /* check for override */
    sprintf(var,"right_NS");
    make_keyword_parameter(&ret,var,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      abortEr("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      abortEr("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      abortEr("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,"grid%u_right_NS_surrounding_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_up",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];

      break;
      case DOWN:
        /* filling name */
        sprintf(name,"grid%u_right_NS_surrounding_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_down",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,"grid%u_right_NS_surrounding_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_left",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,"grid%u_right_NS_surrounding_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_right",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,"grid%u_right_NS_surrounding_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_back",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,"grid%u_right_NS_surrounding_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_NS_surface_function_front",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      default:
        abortEr(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,"grid%u_right_NS_center_a",grid->gn);
    patch->c[0] = GetParameterDoubleF_E(var);
    sprintf(var,"grid%u_right_NS_center_b",grid->gn);
    patch->c[1] = GetParameterDoubleF_E(var);
    sprintf(var,"grid%u_right_NS_center_c",grid->gn);
    patch->c[2] = GetParameterDoubleF_E(var);
    
    /* filling min */
    patch->min[0] = -1;
    patch->min[1] = -1;
    patch->min[2] = 0;
    
    /* filling max */
    patch->max[0] = 1;
    patch->max[1] = 1;
    patch->max[2] = 1;
    
    /* filling flags */
    patch->coordsys = CubedSpherical;
    
    /* collocation */
    patch->collocation[0] = Chebyshev_Extrema;
    patch->collocation[1] = Chebyshev_Extrema;
    patch->collocation[2] = Chebyshev_Extrema;
    
    /* basis */
    patch->basis[0] = Chebyshev_Tn_BASIS;
    patch->basis[1] = Chebyshev_Tn_BASIS;
    patch->basis[2] = Chebyshev_Tn_BASIS;
  }
}

/* populating properties of patch for right BH surrounding */
static void populate_right_BH_surrounding(Grid_T *const grid,const unsigned pn)
{
  unsigned p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R1 = add_field("surface_function",0,patch,NO);
    double *R1_array;
    Flag_T side = p-pn;
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    unsigned n,ijk;
    
    /* filling flags */
    patch->CoordSysInfo->CubedSphericalCoord->side = side;
    patch->CoordSysInfo->CubedSphericalCoord->type = SR_T_CS;
    
    /* filling grid */
    patch->grid = grid;
    
    /* filling patch number */
    patch->pn = p;
    
    /* filling inner boundary */
    patch->innerB = 1;
    
    /* filling n */
    patch->n[0] = (unsigned)GetParameterI("n_a");
    patch->n[1] = (unsigned)GetParameterI("n_b");
    patch->n[2] = (unsigned)GetParameterI("n_c");
    /* check for override */
    sprintf(var,"right_BH");
    make_keyword_parameter(&ret,var,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      abortEr("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      abortEr("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      abortEr("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,"grid%u_right_BH_surrounding_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_BH_surface_function_up",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,"grid%u_right_BH_surrounding_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_BH_surface_function_down",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,"grid%u_right_BH_surrounding_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_BH_surface_function_left",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,"grid%u_right_BH_surrounding_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_BH_surface_function_right",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,"grid%u_right_BH_surrounding_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_BH_surface_function_back",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,"grid%u_right_BH_surrounding_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_right_BH_surface_function_front",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      default:
        abortEr(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,"grid%u_right_BH_center_a",grid->gn);
    patch->c[0] = GetParameterDoubleF_E(var);
    sprintf(var,"grid%u_right_BH_center_b",grid->gn);
    patch->c[1] = GetParameterDoubleF_E(var);
    sprintf(var,"grid%u_right_BH_center_c",grid->gn);
    patch->c[2] = GetParameterDoubleF_E(var);
    
    /* filling min */
    patch->min[0] = -1;
    patch->min[1] = -1;
    patch->min[2] = 0;
    
    /* filling max */
    patch->max[0] = 1;
    patch->max[1] = 1;
    patch->max[2] = 1;
    
    /* filling flags */
    patch->coordsys = CubedSpherical;
    
    /* collocation */
    patch->collocation[0] = Chebyshev_Extrema;
    patch->collocation[1] = Chebyshev_Extrema;
    patch->collocation[2] = Chebyshev_Extrema;
    
    /* basis */
    patch->basis[0] = Chebyshev_Tn_BASIS;
    patch->basis[1] = Chebyshev_Tn_BASIS;
    patch->basis[2] = Chebyshev_Tn_BASIS;
  }
}

/* populating properties of patch for left NS surrounding */
static void populate_left_NS_surrounding(Grid_T *const grid,const unsigned pn)
{
  unsigned p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R1 = add_field("surface_function",0,patch,NO);
    double *R1_array;
    Flag_T side = p-pn;
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    unsigned n,ijk;
    
    /* filling flags */
    patch->CoordSysInfo->CubedSphericalCoord->side = side;
    patch->CoordSysInfo->CubedSphericalCoord->type = SR_T_CS;
    
    /* filling grid */
    patch->grid = grid;
    
    /* filling patch number */
    patch->pn = p;
    
    /* filling inner boundary */
    patch->innerB = 0;
    
    /* filling n */
    patch->n[0] = (unsigned)GetParameterI("n_a");
    patch->n[1] = (unsigned)GetParameterI("n_b");
    patch->n[2] = (unsigned)GetParameterI("n_c");
    /* check for override */
    sprintf(var,"left_NS");
    make_keyword_parameter(&ret,var,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      abortEr("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      abortEr("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      abortEr("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,"grid%u_left_NS_surrounding_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_up",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,"grid%u_left_NS_surrounding_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_down",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,"grid%u_left_NS_surrounding_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_left",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,"grid%u_left_NS_surrounding_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_right",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,"grid%u_left_NS_surrounding_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_back",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,"grid%u_left_NS_surrounding_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,"grid%u_left_NS_surface_function_front",grid->gn);
        R1_array = GetParameterArrayF_E(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,"grid%u_surrounding_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = GetParameterDoubleF_E(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      default:
        abortEr(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,"grid%u_left_NS_center_a",grid->gn);
    patch->c[0] = GetParameterDoubleF_E(var);
    sprintf(var,"grid%u_left_NS_center_b",grid->gn);
    patch->c[1] = GetParameterDoubleF_E(var);
    sprintf(var,"grid%u_left_NS_center_c",grid->gn);
    patch->c[2] = GetParameterDoubleF_E(var);
    
    /* filling min */
    patch->min[0] = -1;
    patch->min[1] = -1;
    patch->min[2] = 0;
    
    /* filling max */
    patch->max[0] = 1;
    patch->max[1] = 1;
    patch->max[2] = 1;
    
    /* filling flags */
    patch->coordsys = CubedSpherical;
    
    /* collocation */
    patch->collocation[0] = Chebyshev_Extrema;
    patch->collocation[1] = Chebyshev_Extrema;
    patch->collocation[2] = Chebyshev_Extrema;
    
    /* basis */
    patch->basis[0] = Chebyshev_Tn_BASIS;
    patch->basis[1] = Chebyshev_Tn_BASIS;
    patch->basis[2] = Chebyshev_Tn_BASIS;
  }
}

/* populating properties of patch for right NS */
static void populate_outermost(Grid_T *const grid,const unsigned pn,const unsigned o)
{
  unsigned p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Flag_T side = p-pn;
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    unsigned n;
    
    /* filling flags */
    patch->CoordSysInfo->CubedSphericalCoord->side = side;
    
    /* filling grid */
    patch->grid = grid;
    
    /* filling patch number */
    patch->pn = p;
    
    /* filling inner boundary */
    patch->innerB = 0;
    
    if (o == 0)
    {
      patch->CoordSysInfo->CubedSphericalCoord->type = OT_T1_CS;
      
      switch(side)
      {
        case UP:
          sprintf(var,"grid%u_surrounding_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = GetParameterDoubleF_E(var); 
        break;
        case DOWN:
          sprintf(var,"grid%u_surrounding_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = -GetParameterDoubleF_E(var); 
        break;
        case LEFT:
          sprintf(var,"grid%u_surrounding_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = -GetParameterDoubleF_E(var);         
        break;
        case RIGHT:
          sprintf(var,"grid%u_surrounding_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = GetParameterDoubleF_E(var);         
        break;
        case BACK:
          sprintf(var,"grid%u_surrounding_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = -GetParameterDoubleF_E(var);           
        break;
        case FRONT:
          sprintf(var,"grid%u_surrounding_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = GetParameterDoubleF_E(var);         
        break;
        default:
          abortEr(NO_OPTION);
      }

      /* filling Rs */
      sprintf(var,"grid%u_outermost%u_R2",grid->gn,o);
      patch->CoordSysInfo->CubedSphericalCoord->R2 = GetParameterDoubleF_E(var);
    
    }
    else
    {
      patch->CoordSysInfo->CubedSphericalCoord->type = OT_T2_CS;
      /* filling Rs */
      sprintf(var,"grid%u_outermost%u_R1",grid->gn,o);
      patch->CoordSysInfo->CubedSphericalCoord->R1 = GetParameterDoubleF_E(var);
      sprintf(var,"grid%u_outermost%u_R2",grid->gn,o);
      patch->CoordSysInfo->CubedSphericalCoord->R2 = GetParameterDoubleF_E(var);
    
    }
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,"grid%u_outermost%u_up",grid->gn,o);
        patch->name = dup_s(name);
        
      break;
      case DOWN:
        /* filling name */
        sprintf(name,"grid%u_outermost%u_down",grid->gn,o);
        patch->name = dup_s(name);
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,"grid%u_outermost%u_left",grid->gn,o);
        patch->name = dup_s(name);
        
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,"grid%u_outermost%u_right",grid->gn,o);
        patch->name = dup_s(name);
        
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,"grid%u_outermost%u_back",grid->gn,o);
        patch->name = dup_s(name);
        
      break;
      case FRONT:
        /* filling name */
        sprintf(name,"grid%u_outermost%u_front",grid->gn,o);
        patch->name = dup_s(name);
      
      break;
      default:
        abortEr(NO_OPTION);
    }

    /* filling n */
    patch->n[0] = (unsigned)GetParameterI("n_a");
    patch->n[1] = (unsigned)GetParameterI("n_b");
    patch->n[2] = (unsigned)GetParameterI("n_c");
    /* check for override */
    sprintf(var,"Outermost%u",o);
    make_keyword_parameter(&ret,var,"n");
    n = (unsigned)GetParameterI(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (unsigned)GetParameterI(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (unsigned)GetParameterI(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      abortEr("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      abortEr("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      abortEr("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    /* filling center */
    patch->c[0] = 0;
    patch->c[1] = 0;
    patch->c[2] = 0;
    
    /* filling min */
    patch->min[0] = -1;
    patch->min[1] = -1;
    patch->min[2] = 0;
    
    /* filling max */
    patch->max[0] = 1;
    patch->max[1] = 1;
    patch->max[2] = 1;
    
    /* filling flags */
    patch->coordsys = CubedSpherical;
    
    /* collocation */
    patch->collocation[0] = Chebyshev_Extrema;
    patch->collocation[1] = Chebyshev_Extrema;
    patch->collocation[2] = Chebyshev_Extrema;
    
    /* basis */
    patch->basis[0] = Chebyshev_Tn_BASIS;
    patch->basis[1] = Chebyshev_Tn_BASIS;
    patch->basis[2] = Chebyshev_Tn_BASIS;
  }
}

/* memory alloc patches for BNS_CubedSpherical type */
void alloc_patches_BNS_CubedSpherical_grid(Grid_T *const grid)
{
  unsigned Np = 30;/* number of patches without outermost's 
                   4 sets of cubed sphere = 4*6
                   4 filling box
                   2 central box */
  unsigned outermost;
  unsigned i;
  
  outermost = (unsigned) GetParameterI("Number_of_Outermost_Split");
  if (outermost != (unsigned)INT_MAX)
    Np += 6*outermost;
  
  grid->patch = calloc((Np+1),sizeof(*grid->patch));
  pointerEr(grid->patch);
  
  for (i = 0; i < Np; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    pointerEr(grid->patch[i]);
  }
  
}

/* memory alloc patches for BBN_CubedSpherical type */
void alloc_patches_BBN_CubedSpherical_grid(Grid_T *const grid)
{
  unsigned Np = 23;/* number of patches without outermost's 
                   3 sets of cubed sphere = 3*6
                   4 filling boxex
                   1 central box */
  unsigned outermost;
  unsigned i;
  
  outermost = (unsigned) GetParameterI("Number_of_Outermost_Split");
  if (outermost != (unsigned)INT_MAX)
    Np += 6*outermost;
  
  grid->patch = calloc((Np+1),sizeof(*grid->patch));
  pointerEr(grid->patch);
  
  for (i = 0; i < Np; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    pointerEr(grid->patch[i]);
  }
  
}

/* memory alloc patches for single neutron star using cubed spherical + box grid */
void alloc_patches_SNS_CubedSpherical_Box_grid(Grid_T *const grid)
{
  unsigned Np = 18;/* number of patches without outermost's 
                   2 sets of cubed sphere = 2*6
                   4 filling boxes
                   2 central boxes */
  unsigned outermost;
  unsigned i;
  
  outermost = (unsigned) GetParameterI("Number_of_Outermost_Split");
  if (outermost != (unsigned)INT_MAX)
    Np += 6*outermost;
  
  grid->patch = calloc((Np+1),sizeof(*grid->patch));
  pointerEr(grid->patch);
  
  for (i = 0; i < Np; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    pointerEr(grid->patch[i]);
  }
  
}

/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_up
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 0,j = 1,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2,dR2_dx,dR2_dy,dR2_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dx->v[p];
      dxc2_dx  = dR2_dx/d1-R2*(X[0]*JT_NS_T_CS_up(patch,_a_,_x_,p)+X[1]*JT_NS_T_CS_up(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dy->v[p];
      dxc2_dy  = dR2_dy/d1-R2*(X[0]*JT_NS_T_CS_up(patch,_a_,_y_,p)+X[1]*JT_NS_T_CS_up(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dz->v[p];
      dxc2_dz  = dR2_dz/d1-R2*(X[0]*JT_NS_T_CS_up(patch,_a_,_z_,p)+X[1]*JT_NS_T_CS_up(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_down
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 1,j = 0,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2,dR2_dx,dR2_dy,dR2_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dx->v[p];
      dxc2_dx  = dR2_dx/d1-R2*(X[0]*JT_NS_T_CS_down(patch,_a_,_x_,p)+X[1]*JT_NS_T_CS_down(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dy->v[p];
      dxc2_dy  = dR2_dy/d1-R2*(X[0]*JT_NS_T_CS_down(patch,_a_,_y_,p)+X[1]*JT_NS_T_CS_down(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dz->v[p];
      dxc2_dz  = dR2_dz/d1-R2*(X[0]*JT_NS_T_CS_down(patch,_a_,_z_,p)+X[1]*JT_NS_T_CS_down(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_left
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 0,j = 2,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2,dR2_dx,dR2_dy,dR2_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dx->v[p];
      dxc2_dx  = dR2_dx/d1-R2*(X[0]*JT_NS_T_CS_left(patch,_a_,_x_,p)+X[1]*JT_NS_T_CS_left(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dy->v[p];
      dxc2_dy  = dR2_dy/d1-R2*(X[0]*JT_NS_T_CS_left(patch,_a_,_y_,p)+X[1]*JT_NS_T_CS_left(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dz->v[p];
      dxc2_dz  = dR2_dz/d1-R2*(X[0]*JT_NS_T_CS_left(patch,_a_,_z_,p)+X[1]*JT_NS_T_CS_left(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_right
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 2,j = 0,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2,dR2_dx,dR2_dy,dR2_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dx->v[p];
      dxc2_dx  = dR2_dx/d1-R2*(X[0]*JT_NS_T_CS_right(patch,_a_,_x_,p)+X[1]*JT_NS_T_CS_right(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dy->v[p];
      dxc2_dy  = dR2_dy/d1-R2*(X[0]*JT_NS_T_CS_right(patch,_a_,_y_,p)+X[1]*JT_NS_T_CS_right(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dz->v[p];
      dxc2_dz  = dR2_dz/d1-R2*(X[0]*JT_NS_T_CS_right(patch,_a_,_z_,p)+X[1]*JT_NS_T_CS_right(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_back
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 2,j = 1,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2,dR2_dx,dR2_dy,dR2_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dx->v[p];
      dxc2_dx  = dR2_dx/d1-R2*(X[0]*JT_NS_T_CS_back(patch,_a_,_x_,p)+X[1]*JT_NS_T_CS_back(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dy->v[p];
      dxc2_dy  = dR2_dy/d1-R2*(X[0]*JT_NS_T_CS_back(patch,_a_,_y_,p)+X[1]*JT_NS_T_CS_back(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dz->v[p];
      dxc2_dz  = dR2_dz/d1-R2*(X[0]*JT_NS_T_CS_back(patch,_a_,_z_,p)+X[1]*JT_NS_T_CS_back(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_front
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 1,j = 2,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2,dR2_dx,dR2_dy,dR2_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dx->v[p];
      dxc2_dx  = dR2_dx/d1-R2*(X[0]*JT_NS_T_CS_front(patch,_a_,_x_,p)+X[1]*JT_NS_T_CS_front(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dy->v[p];
      dxc2_dy  = dR2_dy/d1-R2*(X[0]*JT_NS_T_CS_front(patch,_a_,_y_,p)+X[1]*JT_NS_T_CS_front(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dR2_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR2_dz->v[p];
      dxc2_dz  = dR2_dz/d1-R2*(X[0]*JT_NS_T_CS_front(patch,_a_,_z_,p)+X[1]*JT_NS_T_CS_front(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_up
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 0,j = 1,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc2,
         xc1, dxc1_dx,dxc1_dy,dxc1_dz,
         R1,dR1_dx,dR1_dy,dR1_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dx->v[p];
      dxc1_dx  = dR1_dx/d1-R1*(X[0]*JT_SR_T_CS_up(patch,_a_,_x_,p)+X[1]*JT_SR_T_CS_up(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      J = (K[k==l]-dxc1_dx+(x[k]-xc1)*dxc1_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dy->v[p];
      dxc1_dy  = dR1_dy/d1-R1*(X[0]*JT_SR_T_CS_up(patch,_a_,_y_,p)+X[1]*JT_SR_T_CS_up(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      J = (K[k==l]-dxc1_dy+(x[k]-xc1)*dxc1_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dz->v[p];
      dxc1_dz  = dR1_dz/d1-R1*(X[0]*JT_SR_T_CS_up(patch,_a_,_z_,p)+X[1]*JT_SR_T_CS_up(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      J = (K[k==l]-dxc1_dz+(x[k]-xc1)*dxc1_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_down
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 1,j = 0,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc2,
         xc1, dxc1_dx,dxc1_dy,dxc1_dz,
         R1,dR1_dx,dR1_dy,dR1_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dx->v[p];
      dxc1_dx  = dR1_dx/d1-R1*(X[0]*JT_SR_T_CS_down(patch,_a_,_x_,p)+X[1]*JT_SR_T_CS_down(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      J = (K[k==l]-dxc1_dx+(x[k]-xc1)*dxc1_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dy->v[p];
      dxc1_dy  = dR1_dy/d1-R1*(X[0]*JT_SR_T_CS_down(patch,_a_,_y_,p)+X[1]*JT_SR_T_CS_down(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      J = (K[k==l]-dxc1_dy+(x[k]-xc1)*dxc1_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dz->v[p];
      dxc1_dz  = dR1_dz/d1-R1*(X[0]*JT_SR_T_CS_down(patch,_a_,_z_,p)+X[1]*JT_SR_T_CS_down(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      J = (K[k==l]-dxc1_dz+(x[k]-xc1)*dxc1_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_left
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 0,j = 2,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc2,
         xc1, dxc1_dx,dxc1_dy,dxc1_dz,
         R1,dR1_dx,dR1_dy,dR1_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dx->v[p];
      dxc1_dx  = dR1_dx/d1-R1*(X[0]*JT_SR_T_CS_left(patch,_a_,_x_,p)+X[1]*JT_SR_T_CS_left(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      J = (K[k==l]-dxc1_dx+(x[k]-xc1)*dxc1_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dy->v[p];
      dxc1_dy  = dR1_dy/d1-R1*(X[0]*JT_SR_T_CS_left(patch,_a_,_y_,p)+X[1]*JT_SR_T_CS_left(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      J = (K[k==l]-dxc1_dy+(x[k]-xc1)*dxc1_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dz->v[p];
      dxc1_dz  = dR1_dz/d1-R1*(X[0]*JT_SR_T_CS_left(patch,_a_,_z_,p)+X[1]*JT_SR_T_CS_left(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      J = (K[k==l]-dxc1_dz+(x[k]-xc1)*dxc1_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_right
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 2,j = 0,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc2,
         xc1, dxc1_dx,dxc1_dy,dxc1_dz,
         R1,dR1_dx,dR1_dy,dR1_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dx->v[p];
      dxc1_dx  = dR1_dx/d1-R1*(X[0]*JT_SR_T_CS_right(patch,_a_,_x_,p)+X[1]*JT_SR_T_CS_right(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      J = (K[k==l]-dxc1_dx+(x[k]-xc1)*dxc1_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dy->v[p];
      dxc1_dy  = dR1_dy/d1-R1*(X[0]*JT_SR_T_CS_right(patch,_a_,_y_,p)+X[1]*JT_SR_T_CS_right(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      J = (K[k==l]-dxc1_dy+(x[k]-xc1)*dxc1_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dz->v[p];
      dxc1_dz  = dR1_dz/d1-R1*(X[0]*JT_SR_T_CS_right(patch,_a_,_z_,p)+X[1]*JT_SR_T_CS_right(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      J = (K[k==l]-dxc1_dz+(x[k]-xc1)*dxc1_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_back
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 2,j = 1,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc2,
         xc1, dxc1_dx,dxc1_dy,dxc1_dz,
         R1,dR1_dx,dR1_dy,dR1_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dx->v[p];
      dxc1_dx  = dR1_dx/d1-R1*(X[0]*JT_SR_T_CS_back(patch,_a_,_x_,p)+X[1]*JT_SR_T_CS_back(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      J = (K[k==l]-dxc1_dx+(x[k]-xc1)*dxc1_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dy->v[p];
      dxc1_dy  = dR1_dy/d1-R1*(X[0]*JT_SR_T_CS_back(patch,_a_,_y_,p)+X[1]*JT_SR_T_CS_back(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      J = (K[k==l]-dxc1_dy+(x[k]-xc1)*dxc1_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dz->v[p];
      dxc1_dz  = dR1_dz/d1-R1*(X[0]*JT_SR_T_CS_back(patch,_a_,_z_,p)+X[1]*JT_SR_T_CS_back(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      J = (K[k==l]-dxc1_dz+(x[k]-xc1)*dxc1_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_front
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 1,j = 2,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc2,
         xc1, dxc1_dx,dxc1_dy,dxc1_dz,
         R1,dR1_dx,dR1_dy,dR1_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dx   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dx->v[p];
      dxc1_dx  = dR1_dx/d1-R1*(X[0]*JT_SR_T_CS_front(patch,_a_,_x_,p)+X[1]*JT_SR_T_CS_front(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      J = (K[k==l]-dxc1_dx+(x[k]-xc1)*dxc1_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dy   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dy->v[p];
      dxc1_dy  = dR1_dy/d1-R1*(X[0]*JT_SR_T_CS_front(patch,_a_,_y_,p)+X[1]*JT_SR_T_CS_front(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      J = (K[k==l]-dxc1_dy+(x[k]-xc1)*dxc1_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2;
      xc1 = S*R1/d1;
      dR1_dz   = patch->CoordSysInfo->CubedSphericalCoord->dR1_dz->v[p];
      dxc1_dz  = dR1_dz/d1-R1*(X[0]*JT_SR_T_CS_front(patch,_a_,_z_,p)+X[1]*JT_SR_T_CS_front(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      J = (K[k==l]-dxc1_dz+(x[k]-xc1)*dxc1_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_up
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 0,j = 1,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dx  = -R2*(X[0]*JT_OT_T1_CS_up(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_up(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dy  = -R2*(X[0]*JT_OT_T1_CS_up(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_up(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dz  = -R2*(X[0]*JT_OT_T1_CS_up(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_up(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_down
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 1,j = 0,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dx  = -R2*(X[0]*JT_OT_T1_CS_down(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_down(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dy  = -R2*(X[0]*JT_OT_T1_CS_down(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_down(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dz  = -R2*(X[0]*JT_OT_T1_CS_down(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_down(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_left
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 0,j = 2,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dx  = -R2*(X[0]*JT_OT_T1_CS_left(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_left(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dy  = -R2*(X[0]*JT_OT_T1_CS_left(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_left(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dz  = -R2*(X[0]*JT_OT_T1_CS_left(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_left(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_right
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 2,j = 0,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dx  = -R2*(X[0]*JT_OT_T1_CS_right(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_right(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dy  = -R2*(X[0]*JT_OT_T1_CS_right(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_right(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dz  = -R2*(X[0]*JT_OT_T1_CS_right(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_right(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_back
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 2,j = 1,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dx  = -R2*(X[0]*JT_OT_T1_CS_back(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_back(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dy  = -R2*(X[0]*JT_OT_T1_CS_back(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_back(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dz  = -R2*(X[0]*JT_OT_T1_CS_back(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_back(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_front
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 1,j = 2,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,
         xc2, dxc2_dx,dxc2_dy,dxc2_dz,
         R2;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dx  = -R2*(X[0]*JT_OT_T1_CS_front(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_front(patch,_b_,_x_,p))/d3;
      dxc2_dx *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dx/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dy  = -R2*(X[0]*JT_OT_T1_CS_front(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_front(patch,_b_,_y_,p))/d3;
      dxc2_dy *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dy/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      xc2 = S*R2/d1;
      dxc2_dz  = -R2*(X[0]*JT_OT_T1_CS_front(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_front(patch,_b_,_z_,p))/d3;
      dxc2_dz *= S;
      J = (K[k==l]-(x[k]-xc1)*dxc2_dz/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_up
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 0,j = 1,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,dxc1_dx,dxc1_dy,dxc1_dz,
         xc2,dxc2_dx,dxc2_dy,dxc2_dz,
         R1,R2;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dx  = -R1*(X[0]*JT_OT_T2_CS_up(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_up(patch,_b_,_x_,p))/d3;
      dxc2_dx  = -R2*(X[0]*JT_OT_T2_CS_up(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_up(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      dxc2_dx *= S;
      J = (K[k==l]-dxc1_dx -(x[k]-xc1)*(dxc2_dx-dxc1_dx)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dy  = -R1*(X[0]*JT_OT_T2_CS_up(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_up(patch,_b_,_y_,p))/d3;
      dxc2_dy  = -R2*(X[0]*JT_OT_T2_CS_up(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_up(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      dxc2_dy *= S;
      J = (K[k==l]-dxc1_dy -(x[k]-xc1)*(dxc2_dy-dxc1_dy)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dz  = -R1*(X[0]*JT_OT_T2_CS_up(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_up(patch,_b_,_z_,p))/d3;
      dxc2_dz  = -R2*(X[0]*JT_OT_T2_CS_up(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_up(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      dxc2_dz *= S;
      J = (K[k==l]-dxc1_dz -(x[k]-xc1)*(dxc2_dz-dxc1_dz)/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_down
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 1,j = 0,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,dxc1_dx,dxc1_dy,dxc1_dz,
         xc2,dxc2_dx,dxc2_dy,dxc2_dz,
         R1,R2;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dx  = -R1*(X[0]*JT_OT_T2_CS_down(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_down(patch,_b_,_x_,p))/d3;
      dxc2_dx  = -R2*(X[0]*JT_OT_T2_CS_down(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_down(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      dxc2_dx *= S;
      J = (K[k==l]-dxc1_dx -(x[k]-xc1)*(dxc2_dx-dxc1_dx)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dy  = -R1*(X[0]*JT_OT_T2_CS_down(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_down(patch,_b_,_y_,p))/d3;
      dxc2_dy  = -R2*(X[0]*JT_OT_T2_CS_down(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_down(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      dxc2_dy *= S;
      J = (K[k==l]-dxc1_dy -(x[k]-xc1)*(dxc2_dy-dxc1_dy)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dz  = -R1*(X[0]*JT_OT_T2_CS_down(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_down(patch,_b_,_z_,p))/d3;
      dxc2_dz  = -R2*(X[0]*JT_OT_T2_CS_down(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_down(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      dxc2_dz *= S;
      J = (K[k==l]-dxc1_dz -(x[k]-xc1)*(dxc2_dz-dxc1_dz)/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_left
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 0,j = 2,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,dxc1_dx,dxc1_dy,dxc1_dz,
         xc2,dxc2_dx,dxc2_dy,dxc2_dz,
         R1,R2;  
         
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dx  = -R1*(X[0]*JT_OT_T2_CS_left(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_left(patch,_b_,_x_,p))/d3;
      dxc2_dx  = -R2*(X[0]*JT_OT_T2_CS_left(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_left(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      dxc2_dx *= S;
      J = (K[k==l]-dxc1_dx -(x[k]-xc1)*(dxc2_dx-dxc1_dx)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dy  = -R1*(X[0]*JT_OT_T2_CS_left(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_left(patch,_b_,_y_,p))/d3;
      dxc2_dy  = -R2*(X[0]*JT_OT_T2_CS_left(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_left(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      dxc2_dy *= S;
      J = (K[k==l]-dxc1_dy -(x[k]-xc1)*(dxc2_dy-dxc1_dy)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dz  = -R1*(X[0]*JT_OT_T2_CS_left(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_left(patch,_b_,_z_,p))/d3;
      dxc2_dz  = -R2*(X[0]*JT_OT_T2_CS_left(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_left(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      dxc2_dz *= S;
      J = (K[k==l]-dxc1_dz -(x[k]-xc1)*(dxc2_dz-dxc1_dz)/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_right
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 2,j = 0,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,dxc1_dx,dxc1_dy,dxc1_dz,
         xc2,dxc2_dx,dxc2_dy,dxc2_dz,
         R1,R2;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dx  = -R1*(X[0]*JT_OT_T2_CS_right(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_right(patch,_b_,_x_,p))/d3;
      dxc2_dx  = -R2*(X[0]*JT_OT_T2_CS_right(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_right(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      dxc2_dx *= S;
      J = (K[k==l]-dxc1_dx -(x[k]-xc1)*(dxc2_dx-dxc1_dx)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dy  = -R1*(X[0]*JT_OT_T2_CS_right(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_right(patch,_b_,_y_,p))/d3;
      dxc2_dy  = -R2*(X[0]*JT_OT_T2_CS_right(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_right(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      dxc2_dy *= S;
      J = (K[k==l]-dxc1_dy -(x[k]-xc1)*(dxc2_dy-dxc1_dy)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dz  = -R1*(X[0]*JT_OT_T2_CS_right(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_right(patch,_b_,_z_,p))/d3;
      dxc2_dz  = -R2*(X[0]*JT_OT_T2_CS_right(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_right(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      dxc2_dz *= S;
      J = (K[k==l]-dxc1_dz -(x[k]-xc1)*(dxc2_dz-dxc1_dz)/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_back
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const unsigned i = 2,j = 1,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,dxc1_dx,dxc1_dy,dxc1_dz,
         xc2,dxc2_dx,dxc2_dy,dxc2_dz,
         R1,R2;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dx  = -R1*(X[0]*JT_OT_T2_CS_back(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_back(patch,_b_,_x_,p))/d3;
      dxc2_dx  = -R2*(X[0]*JT_OT_T2_CS_back(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_back(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      dxc2_dx *= S;
      J = (K[k==l]-dxc1_dx -(x[k]-xc1)*(dxc2_dx-dxc1_dx)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dy  = -R1*(X[0]*JT_OT_T2_CS_back(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_back(patch,_b_,_y_,p))/d3;
      dxc2_dy  = -R2*(X[0]*JT_OT_T2_CS_back(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_back(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      dxc2_dy *= S;
      J = (K[k==l]-dxc1_dy -(x[k]-xc1)*(dxc2_dy-dxc1_dy)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dz  = -R1*(X[0]*JT_OT_T2_CS_back(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_back(patch,_b_,_z_,p))/d3;
      dxc2_dz  = -R2*(X[0]*JT_OT_T2_CS_back(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_back(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      dxc2_dz *= S;
      J = (K[k==l]-dxc1_dz -(x[k]-xc1)*(dxc2_dz-dxc1_dz)/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_front
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const unsigned i = 1,j = 2,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  unsigned l;
  double xc1,dxc1_dx,dxc1_dy,dxc1_dz,
         xc2,dxc2_dx,dxc2_dy,dxc2_dz,
         R1,R2;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(K[i==l]-x[i]*K[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(K[j==l]-x[j]*K[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dx  = -R1*(X[0]*JT_OT_T2_CS_front(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_front(patch,_b_,_x_,p))/d3;
      dxc2_dx  = -R2*(X[0]*JT_OT_T2_CS_front(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_front(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      dxc2_dx *= S;
      J = (K[k==l]-dxc1_dx -(x[k]-xc1)*(dxc2_dx-dxc1_dx)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dy  = -R1*(X[0]*JT_OT_T2_CS_front(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_front(patch,_b_,_y_,p))/d3;
      dxc2_dy  = -R2*(X[0]*JT_OT_T2_CS_front(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_front(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      dxc2_dy *= S;
      J = (K[k==l]-dxc1_dy -(x[k]-xc1)*(dxc2_dy-dxc1_dy)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+SQR(X[0])+SQR(X[1]));
      d3 = Power3(d1);
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1 = S*R1/d1;
      xc2 = S*R2/d1;
      dxc1_dz  = -R1*(X[0]*JT_OT_T2_CS_front(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_front(patch,_b_,_z_,p))/d3;
      dxc2_dz  = -R2*(X[0]*JT_OT_T2_CS_front(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_front(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      dxc2_dz *= S;
      J = (K[k==l]-dxc1_dz -(x[k]-xc1)*(dxc2_dz-dxc1_dz)/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      abortEr("No such an enum!\n");
  }
  
  return J;
}

/* setting up sign and indices for cubde sphere based on side*/  
void SignAndIndex_permutation_CubedSphere(const Flag_T side,unsigned *const a,unsigned *const b,unsigned *const c,double *const s)
{
  switch (side)
  {
    case UP:
      *s = 1;
      *a = 0;
      *b = 1;
      *c = 2;
    break;
    case DOWN:
      *s = -1;
      *a = 1;
      *b = 0;
      *c = 2;
    break;
    case LEFT:
      *s = -1;
      *a = 0;
      *b = 2;
      *c = 1;
    break;
    case RIGHT:
      *s = 1;
      *a = 2;
      *b = 0;
      *c = 1;
    break;
    case BACK:
      *s = -1;
      *a = 2;
      *b = 1;
      *c = 0;
    break;
    case FRONT:
      *s = 1;
      *a = 1;
      *b = 2;
      *c = 0;
    break;
    default:
      abortEr(NO_OPTION);
  }
}

/* performing a test for CS coordinates:
// test includes:
// o. 2-D interpolations of radius.
// o. x_of_X function
// o. X_of_x function. */
void test_CubedSpherical_Coordinates(Grid_T *const grid)
{
  unsigned p;
  Flag_T flg = NONE;
  
  printf("Testing Cubed Spherical Coordinates:\n");
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch 	= grid->patch[p];
    const Flag_T type   = patch->CoordSysInfo->CubedSphericalCoord->type;
    Field_T *const R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f,
            *const R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
    double Xp[3],xp[3],R;
    unsigned n;
    
    for(n = 0; n < patch->nn; ++n)
    {
      const double *const X = patch->node[n]->X;
      const double *const x = patch->node[n]->x;
      
      /* test x_of_X */
      x_of_X(xp,X,patch);
      if (!EQL(L2_norm(3,xp,x),0))
      {
        printf("x_of_X failed.\n");
        flg = FOUND;
        break;
      }
      
      /* test X_of_x */
      X_of_x(Xp,x,patch);
      if (!EQL(L2_norm(3,Xp,X),0))
      {
        printf("X_of_x failed.\n");
        flg = FOUND;
        break;
      }
      
      /* test Radius related */
      switch (type)
      {
        case NS_T_CS:
          R = R_interpolation_CS(R2_f,X);
          if (!EQL(L2_norm(1,&R,&R2_f->v[n]),0))
          {
            printf("R interpolation failed.\n");
            flg = FOUND;
            break;
          }
          
        break;
        case SR_T_CS:
          R = R_interpolation_CS(R1_f,X);
          if (!EQL(L2_norm(1,&R,&R1_f->v[n]),0))
          {
            printf("R interpolation failed.\n");
            flg = FOUND;
            break;
          }
        break;
        default:
        break;
      }
  
    }/* end of for(n = 0; n < patch->nn; ++n) */
    if (flg == FOUND)
      break;
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
  if (flg != FOUND)
    printf("Testing Cubed Spherical Coordinates: [+].\n");
  else
    printf("Testing Cubed Spherical Coordinates: [-].\n");
}
