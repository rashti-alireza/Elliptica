/*
// Alireza Rashti
// April 2019
*/

#include "cubed_spherical_coordinate.h"

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

/* making Jacobian transformation for ProjectiveHemisphere coord.
void make_JacobianT_ProjectiveHemisphere_coord(Patch_T *const patch)
{
  patch->JacobianT->j      = JT_ProjectiveHemisphere;
  patch->JacobianT->dN0_dx = dN0_dx_ProjectiveHemisphere_patch;
  patch->JacobianT->dN0_dy = dN0_dy_ProjectiveHemisphere_patch;
  patch->JacobianT->dN0_dz = dN0_dz_ProjectiveHemisphere_patch;
  patch->JacobianT->dN1_dx = dN1_dx_ProjectiveHemisphere_patch;
  patch->JacobianT->dN1_dy = dN1_dy_ProjectiveHemisphere_patch;
  patch->JacobianT->dN1_dz = dN1_dz_ProjectiveHemisphere_patch;
  patch->JacobianT->dN2_dx = dN2_dx_ProjectiveHemisphere_patch;
  patch->JacobianT->dN2_dy = dN2_dy_ProjectiveHemisphere_patch;
  patch->JacobianT->dN2_dz = dN2_dz_ProjectiveHemisphere_patch;
  R1_R2_derivative(patch);
}
*/
/* preparing R1 and R2 derivatives of projective hemisphere coords.
// NOTE: one must remove dR_? after each updating of R.
static void R1_R2_derivative(Patch_T *const patch)
{
  Field_T *dR1_dX = add_field("dR1_dX",0,patch,NO),
          *dR1_dY = add_field("dR1_dY",0,patch,NO),
          *dR2_dX = add_field("dR2_dX",0,patch,NO),
          *dR2_dY = add_field("dR2_dY",0,patch,NO),
          *dR1_dx = add_field("dR1_dx",0,patch,YES),
          *dR1_dy = add_field("dR1_dy",0,patch,YES),
          *dR1_dz = add_field("dR1_dz",0,patch,YES),
          *dR2_dx = add_field("dR2_dx",0,patch,YES),
          *dR2_dy = add_field("dR2_dy",0,patch,YES),
          *dR2_dz = add_field("dR2_dz",0,patch,YES);
  Field_T *const R1 = patch->pool[Ind("R1_ProjectiveHemisphere")];
  Field_T *const R2 = patch->pool[Ind("R2_ProjectiveHemisphere")];
  const unsigned *const n = patch->n;
  unsigned i,j;
          
  dR1_dX->v = Partial_Derivative(R1,"a");
  dR1_dY->v = Partial_Derivative(R1,"b");
  dR2_dX->v = Partial_Derivative(R2,"a");
  dR2_dY->v = Partial_Derivative(R2,"b");
    
  for (i = 0; i < n[0]; ++i)
    for (j = 0; j < n[1]; ++j)
    {
      unsigned p = L(n,i,j,0);
      dR1_dx->v[p] = dR1_dX->v[p]*dq2_dq1(patch,_a_,_x_,p)+
                     dR1_dY->v[p]*dq2_dq1(patch,_b_,_x_,p);
      dR1_dy->v[p] = dR1_dX->v[p]*dq2_dq1(patch,_a_,_y_,p)+
                     dR1_dY->v[p]*dq2_dq1(patch,_b_,_y_,p);
      dR1_dz->v[p] = dR1_dX->v[p]*dq2_dq1(patch,_a_,_z_,p)+
                     dR1_dY->v[p]*dq2_dq1(patch,_b_,_z_,p);
      dR2_dx->v[p] = dR2_dX->v[p]*dq2_dq1(patch,_a_,_x_,p)+
                     dR2_dY->v[p]*dq2_dq1(patch,_b_,_x_,p);
      dR2_dy->v[p] = dR2_dX->v[p]*dq2_dq1(patch,_a_,_y_,p)+
                     dR2_dY->v[p]*dq2_dq1(patch,_b_,_y_,p);
      dR2_dz->v[p] = dR2_dX->v[p]*dq2_dq1(patch,_a_,_z_,p)+
                     dR2_dY->v[p]*dq2_dq1(patch,_b_,_z_,p);
      
    }
                      
  remove_field(dR1_dX);
  remove_field(dR1_dY);
  remove_field(dR2_dX);
  remove_field(dR2_dY);
  
  patch->CoordSysInfo->ProjectiveCoord->dR1_dx = dR1_dx;
  patch->CoordSysInfo->ProjectiveCoord->dR1_dy = dR1_dy;
  patch->CoordSysInfo->ProjectiveCoord->dR1_dz = dR1_dz;
  patch->CoordSysInfo->ProjectiveCoord->dR2_dx = dR2_dx;
  patch->CoordSysInfo->ProjectiveCoord->dR2_dy = dR2_dy;
  patch->CoordSysInfo->ProjectiveCoord->dR2_dz = dR2_dz; 
}*/

/* Jacobian transformation for dN/dX?.
// ->return value: dN/dX?
static double dN_dX(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e)
{
  double jN_X = 0;

  if (patch->collocation[q1_e%3] == Chebyshev_Extrema)
  {
    if (q2_e%3 == q1_e%3)
      jN_X = 2.0/(-patch->max[q1_e%3]+patch->min[q1_e%3]); 
    else
      jN_X = 0;
  }
  else if (patch->collocation[q1_e%3] == Chebyshev_Nodes)
  {
    if (q2_e%3 == q1_e%3)
      jN_X = 2.0/(-patch->max[q1_e%3]+patch->min[q1_e%3]); 
    else
      jN_X = 0;
  }
  else
  {
    abortEr(INCOMPLETE_FUNC);
  }

  return jN_X;
}
*/

/* interpolation of 2d field for Projective Hemisphere coordinate.
// Note in Projective Hemisphere 2d fields are f(X,Y,0).
// it mostly is used for radius or derivative of radius.
// convection:
// R = interesting field
// patch = patch that has R
// X = the curvilinear coord in which we want R(X).
// ->return value: R(X)
double interpolation_2d_PH(Field_T *const R, const Patch_T *const patch,const double *const X)
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
  
  UNUSED(patch);
  
  return interp;
}
*/

/* populating properties of patch for left NS */
static void populate_left_NS(Grid_T *const grid,const unsigned pn)
{
  unsigned p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R2 = add_field("NS_surface",0,patch,NO);
    double *R2_array;
    Flag_T side = p-pn;
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    unsigned n,i,j;
    
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
    Field_T *R2 = add_field("NS_surface",0,patch,NO);
    double *R2_array;
    Flag_T side = p-pn;
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    unsigned n,i,j;
    
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R2->v[ij0] = R2_array[ij0];
          }
      
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
    Field_T *R1 = add_field("NS_surface",0,patch,NO);
    double *R1_array;
    Flag_T side = p-pn;
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    unsigned n,i,j;
    
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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

/* populating properties of patch for left NS surrounding */
static void populate_left_NS_surrounding(Grid_T *const grid,const unsigned pn)
{
  unsigned p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R1 = add_field("NS_surface",0,patch,NO);
    double *R1_array;
    Flag_T side = p-pn;
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    unsigned n,i,j;
    
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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
        
        for (i = 0; i < patch->n[0]; ++i)
          for (j = 0; j < patch->n[1]; ++j)
          {
            unsigned ij0 = L(patch->n,i,j,0);
            R1->v[ij0] = R1_array[ij0];
          }
      
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

/* memory alloc patches for BNS_Projective type */
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

/* Jacobian transformation for ProjectiveHemisphereUp patch.
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1
*/
//double JT_ProjectiveHemisphere(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p)
//{
  /* ds/ds = 1 */
  //if (q2_e == q1_e)
    //return 1;
  /* dx/d? not defined */
 // else if (q2_e == _x_ || q2_e == _y_|| q2_e == _z_)
   // abortEr(INCOMPLETE_FUNC);
  /* dx/d? not defined */
 /* else if (q1_e == _a_ || q1_e == _b_|| q1_e == _c_)
    abortEr(INCOMPLETE_FUNC);
  
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const c = patch->c;
  const double 
         x = patch->node[p]->x[0]-c[0],
         y = patch->node[p]->x[1]-c[1],
         z = patch->node[p]->x[2]-c[2],
         r2 = SQR(x)+SQR(y)+SQR(z),
         r1 = sqrt(r2);
  const unsigned *const n = patch->n;
  double S,dS_dx,dS_dy,dS_dz,
         u,v,r3,
         du_dx,du_dy,du_dz,
         dv_dx,dv_dy,dv_dz,
         R1,R2,dR2_dx,dR2_dy,dR2_dz,
         dR1_dx,dR1_dy,dR1_dz;
  */
  /* if dZ_d? we need dR_d? thus: */
  /*
  if (q2_e == _c_)
  {
    unsigned i,j,k;
    IJK(p,n,&i,&j,&k);
    
    R1 = patch->CoordSysInfo->ProjectiveCoord->R1_f->v[L(n,i,j,0)];
    R2 = patch->CoordSysInfo->ProjectiveCoord->R2_f->v[L(n,i,j,0)];
    dR2_dx = patch->CoordSysInfo->ProjectiveCoord->dR2_dx->v[L(n,i,j,0)];
    dR2_dy = patch->CoordSysInfo->ProjectiveCoord->dR2_dy->v[L(n,i,j,0)];
    dR2_dz = patch->CoordSysInfo->ProjectiveCoord->dR2_dz->v[L(n,i,j,0)];
    dR1_dx = patch->CoordSysInfo->ProjectiveCoord->dR1_dx->v[L(n,i,j,0)];
    dR1_dy = patch->CoordSysInfo->ProjectiveCoord->dR1_dy->v[L(n,i,j,0)];
    dR1_dz = patch->CoordSysInfo->ProjectiveCoord->dR1_dz->v[L(n,i,j,0)];
  }
  else
  {
    r3 = Power3(r1);
    S = 1./r1;
    dS_dx = -x/r3;
    dS_dy = -y/r3;
    dS_dz = -z/r3;
    u = x/r1;
    v = y/r1;
  }
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      du_dx = x*dS_dx+S;
      dv_dx = y*dS_dx;
      J = dX_du(u,v)*du_dx+dX_dv(u,v)*dv_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dy:
      du_dy = x*dS_dy;
      dv_dy = y*dS_dy+S;
      J = dX_du(u,v)*du_dy+dX_dv(u,v)*dv_dy;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case da_dz:
      du_dz = x*dS_dz;
      dv_dz = y*dS_dz;
      J = dX_du(u,v)*du_dz+dX_dv(u,v)*dv_dz;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dx:
      du_dx = x*dS_dx+S;
      dv_dx = y*dS_dx;
      J = dY_du(u,v)*du_dx+dY_dv(u,v)*dv_dx;
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dy:
      du_dy = x*dS_dy;
      dv_dy = y*dS_dy+S;
      J = dY_du(u,v)*du_dy+dY_dv(u,v)*dv_dy;
            
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case db_dz:
      du_dz = x*dS_dz;
      dv_dz = y*dS_dz;
      J = dY_du(u,v)*du_dz+dY_dv(u,v)*dv_dz;
            
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dx:
      J = 2*(x/r1/(R2-R1) - (r1*(dR2_dx-dR1_dx)-R2*dR1_dx+R1*dR2_dx)/SQR(R2-R1));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
    break;
    case dc_dy:
      J = 2*(y/r1/(R2-R1) - (r1*(dR2_dy-dR1_dy)-R2*dR1_dy+R1*dR2_dy)/SQR(R2-R1));
          
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up."); 
    break;
    case dc_dz:
      J = 2*(z/r1/(R2-R1) - (r1*(dR2_dz-dR1_dz)-R2*dR1_dz+R1*dR2_dz)/SQR(R2-R1));
      
      if (fpclassify(J) == FP_NAN || fpclassify(J) == FP_INFINITE)
        abortEr("Jacobian is messed up.");
      
    break;
    default:
      abortEr("No such a enum!\n");
  }
  
  return J;
}*/

