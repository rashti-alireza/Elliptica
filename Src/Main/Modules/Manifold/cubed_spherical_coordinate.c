/*
// Alireza Rashti
// April 2019
*/

#include "cubed_spherical_coordinate.h"

/* filling with cubed spherical coordinate patches for SBH grid */
void fill_patches_SBH_CubedSpherical_grid(Grid_T *const grid)
{
  const Uint N_outermost_split = (Uint) Pgeti("Number_of_Outermost_Split");
  Uint i,pn;
  
  pn = 0;
  populate_central_BH_around(grid,pn);
  pn += 6; /* +6 cubed sphere */
  
  for (i = 0; i < N_outermost_split; i++)
  {
    populate_outermost(grid,pn,i);
    pn += 6; /* +6 cubed sphere */
  }

}

/* filling with cubed spherical coordinate patches for SNS grid */
void fill_patches_SNS_CubedSpherical_grid(Grid_T *const grid)
{
  const Uint N_outermost_split = (Uint) Pgeti("Number_of_Outermost_Split");
  Uint i,pn;
  
  pn = 0; /* patch number */
  populate_central_NS_central_box(grid,pn++);/* +1 */
  populate_central_NS(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_central_NS_around(grid,pn);
  pn += 6; /* +6 cubed sphere */
  
  for (i = 0; i < N_outermost_split; i++)
  {
    populate_outermost(grid,pn,i);
    pn += 6; /* +6 cubed sphere */
  }

}

/* filling with cubed spherical + box coordinate patches for SNS grid */
void fill_patches_SNS_CubedSpherical_Box_grid(Grid_T *const grid)
{
  const Uint N_outermost_split = (Uint) Pgeti("Number_of_Outermost_Split");
  Uint i,pn;
  
  pn = 0; /* patch number */
  populate_left_NS_central_box(grid,pn++);/* +1 */
  populate_left_NS(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_left_NS_around(grid,pn);
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
  const Uint N_outermost_split = (Uint) Pgeti("Number_of_Outermost_Split");
  Uint i,pn;
  
  pn = 0; /* patch number */
  populate_left_NS_central_box(grid,pn++);/* +1 */
  populate_left_NS(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_left_NS_around(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_right_NS_central_box(grid,pn++);/*+1 */
  populate_right_NS(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_right_NS_around(grid,pn);
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
  const Uint N_outermost_split = (Uint) Pgeti("Number_of_Outermost_Split");
  Uint i,pn;
  
  pn = 0; /* patch number */
  populate_left_NS_central_box(grid,pn++);/* +1 */
  populate_left_NS(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_left_NS_around(grid,pn);
  pn += 6; /* +6 cubed sphere */
  populate_right_BH_around(grid,pn);
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

/* allocating and filling grid->patch for split cubed spherical.
// when one patch allocated, grid->np also added by 1. */
void fill_patches_Split_CubedSpherical_grid(Grid_T *const grid)
{
  const double r_outermost = Pgetd("grid_outermost_radius");
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS)
  {
    Flag_T ns_side = NONE, bh_side = NONE;
    Flag_T bh_filled = NONE;
    
    if (Pcmps("grid_set_NS","left"))
    {
      ns_side = LEFT;
      bh_side = RIGHT;
    }
    else
    {
      bh_side = LEFT;
      ns_side = RIGHT;
    }
    if (strstr_i(Pgets("grid_set_BH"),"excised"))
    {
      bh_filled = NO;
    }
    else if (strstr_i(Pgets("grid_set_BH"),"filled"))
    {
      bh_filled = YES;
    }
    else
    {
      Error0(NO_OPTION);
    }
    
    /* boxes */ 
    populate_box_patch_SplitCS(grid,"central_box",ns_side,"NS");
    /* cubed sphericals */
    populate_CS_patch_SplitCS(grid,"NS",ns_side);
    populate_CS_patch_SplitCS(grid,"NS_around",ns_side);
    populate_CS_patch_SplitCS(grid,"BH_around",bh_side);
    
    #if USE_SCS_FILLING_BOX == 1
    populate_box_patch_SplitCS(grid,"filling_box",UP,"filling_box");
    populate_box_patch_SplitCS(grid,"filling_box",DOWN,"filling_box");
    populate_box_patch_SplitCS(grid,"filling_box",BACK,"filling_box");
    populate_box_patch_SplitCS(grid,"filling_box",FRONT,"filling_box");
    if (!EQL(r_outermost,0))
      populate_CS_patch_SplitCS(grid,"outermost",NONE);
    #endif
    
    if (!EQL(r_outermost,0))
    {
      populate_CS_patch_SplitCS(grid,"outermost",LEFT);
      populate_CS_patch_SplitCS(grid,"outermost",RIGHT);
    }
    
    /* NOTE: order matters */
    if (bh_filled == YES)
    {
      populate_CS_patch_SplitCS(grid,"BH",bh_side);
      populate_box_patch_SplitCS(grid,"central_box",bh_side ,"BH");
    }
    else
    {
      /* set innerB for BH_around */
      Uint nbh = 0,p;
      Patch_T **patches = 
        collect_patches(grid,"BH_around_IB",&nbh);
      
      for (p = 0; p < nbh; ++p)
      {
        Patch_T *patch = patches[p];
        patch->innerB = 1;
      }
      
      Free(patches);
    }
    
  }
  else if (grid->kind == Grid_SplitCubedSpherical_NSNS)
  {
    Flag_T ns_side1 = NONE, ns_side2 = NONE;
    
    if (Pcmps("grid_set_NS1","left"))
    {
      ns_side1 = LEFT;
      ns_side2 = RIGHT;
    }
    else
    {
      ns_side2 = LEFT;
      ns_side1 = RIGHT;
    }
    
    /* boxes */ 
    populate_box_patch_SplitCS(grid,"central_box",ns_side1 ,"NS1");
    populate_box_patch_SplitCS(grid,"central_box",ns_side2 ,"NS2");
   
    /* cubed sphericals */
    populate_CS_patch_SplitCS(grid,"NS1",ns_side1);
    populate_CS_patch_SplitCS(grid,"NS1_around",ns_side1);
    populate_CS_patch_SplitCS(grid,"NS2",ns_side2);
    populate_CS_patch_SplitCS(grid,"NS2_around",ns_side2);
    
    #if USE_SCS_FILLING_BOX == 1
    populate_box_patch_SplitCS(grid,"filling_box",UP   ,"filling_box");
    populate_box_patch_SplitCS(grid,"filling_box",DOWN ,"filling_box");
    populate_box_patch_SplitCS(grid,"filling_box",BACK ,"filling_box");
    populate_box_patch_SplitCS(grid,"filling_box",FRONT,"filling_box");
    if (!EQL(r_outermost,0))
      populate_CS_patch_SplitCS(grid,"outermost",NONE);
    #endif

    if (!EQL(r_outermost,0))
    {
      populate_CS_patch_SplitCS(grid,"outermost",LEFT);
      populate_CS_patch_SplitCS(grid,"outermost",RIGHT);
    }
 
  }
  else if (grid->kind == Grid_SplitCubedSpherical_BHBH)
  {
    Flag_T bh_side1 = NONE, bh_side2 = NONE;
    Flag_T bh_filled = NONE;
    
    if (strstr_i(Pgets("grid_set_BH1"),"left"))
    {
      bh_side1 = LEFT;
      bh_side2 = RIGHT;
    }
    else
    {
      bh_side2 = LEFT;
      bh_side1 = RIGHT;
    }
    if (strstr_i(Pgets("grid_set_BH1"),"excised") &&
        strstr_i(Pgets("grid_set_BH2"),"excised"))
    {
      bh_filled = NO;
    }
    else if (strstr_i(Pgets("grid_set_BH1"),"filled") &&
             strstr_i(Pgets("grid_set_BH2"),"filled"))
    {
      bh_filled = YES;
    }
    else
    {
      Error0(NO_OPTION);
    }
    
    /* cubed sphericals */
    populate_CS_patch_SplitCS(grid,"BH1_around",bh_side1);
    populate_CS_patch_SplitCS(grid,"BH2_around",bh_side2);
    
    #if USE_SCS_FILLING_BOX == 1
    /* boxes */ 
    populate_box_patch_SplitCS(grid,"filling_box",UP   ,"filling_box");
    populate_box_patch_SplitCS(grid,"filling_box",DOWN ,"filling_box");
    populate_box_patch_SplitCS(grid,"filling_box",BACK ,"filling_box");
    populate_box_patch_SplitCS(grid,"filling_box",FRONT,"filling_box");
    if (!EQL(r_outermost,0))
      populate_CS_patch_SplitCS(grid,"outermost",NONE);
    #endif
    
    if (!EQL(r_outermost,0))
    {
      populate_CS_patch_SplitCS(grid,"outermost",LEFT);
      populate_CS_patch_SplitCS(grid,"outermost",RIGHT);
    }
 
    /* NOTE: order matters */
    if (bh_filled == YES)
    {
      populate_CS_patch_SplitCS(grid,"BH1",bh_side1);
      populate_box_patch_SplitCS(grid,"central_box",bh_side1 ,"BH1");
      populate_CS_patch_SplitCS(grid,"BH2",bh_side2);
      populate_box_patch_SplitCS(grid,"central_box",bh_side2 ,"BH2");
    }
    else
    {
      /* set innerB for BH_around */
      Uint nbh = 0,p;
      Patch_T **patches = 0;
     
      /* BH1 */
      patches = collect_patches(grid,"BH1_around_IB",&nbh);
      for (p = 0; p < nbh; ++p)
      {
        Patch_T *patch = patches[p];
        patch->innerB = 1;
      }
      Free(patches);
      
      /* BH 2 */
      patches = collect_patches(grid,"BH2_around_IB",&nbh);
      for (p = 0; p < nbh; ++p)
      {
        Patch_T *patch = patches[p];
        patch->innerB = 1;
      }
      Free(patches);
    }
  }
  else if (grid->kind == Grid_SplitCubedSpherical_SBH)
  {
    Flag_T bh_side    = NONE;
    Flag_T bh_filled = NONE;
    
    if (strstr_i(Pgets("grid_set_BH"),"center"))
    {
      bh_side = CENTER;
    }
    else
    {
      Error0(NO_OPTION);
    }
    if (strstr_i(Pgets("grid_set_BH"),"excised"))
    {
      bh_filled = NO;
    }
    else if (strstr_i(Pgets("grid_set_BH"),"filled"))
    {
      bh_filled = YES;
    }
    else
    {
      Error0(NO_OPTION);
    }
    
    /* cubed sphericals */
    populate_CS_patch_SplitCS(grid,"BH_around",bh_side);
    
    if (!EQL(r_outermost,0))
      populate_CS_patch_SplitCS(grid,"outermost",NONE);
    
    /* NOTE: order matters */
    if (bh_filled == YES)
    {
      populate_CS_patch_SplitCS(grid,"BH",bh_side);
      populate_box_patch_SplitCS(grid,"central_box",bh_side ,"BH");
    }
    else
    {
      /* set innerB for BH_around */
      Uint nbh = 0,p;
      Patch_T **patches = 0;
     
      /* BH */
      patches = collect_patches(grid,"BH_around_IB",&nbh);
      for (p = 0; p < nbh; ++p)
      {
        Patch_T *patch = patches[p];
        patch->innerB = 1;
      }
      Free(patches);
    }
    
  }
  else if (grid->kind == Grid_SplitCubedSpherical_SNS)
  {
    Flag_T ns_side = NONE;
    
    if (Pcmps("grid_set_NS","center"))
    {
      ns_side = CENTER;
    }
    else
    {
      Error0(NO_OPTION);
    }
    
    /* cubed sphericals */
    populate_CS_patch_SplitCS(grid,"NS_around",ns_side);
    populate_CS_patch_SplitCS(grid,"NS",ns_side);
    /* box */
    populate_box_patch_SplitCS(grid,"central_box",ns_side ,"NS");
    
    if (!EQL(r_outermost,0))
      populate_CS_patch_SplitCS(grid,"outermost",NONE);
  }
  else
  {
    Error0(NO_OPTION);
  }
}

/* populating properties of a patch for a split cubed spherical object,
// like: BH, NS, arounds, outermost. */
void 
populate_CS_patch_SplitCS
  (
  Grid_T *const grid,
  const char *const obj0,/* NS, BH or etc. */
  const Flag_T dir0/* LEFT or RIGHT or CENTER or NONE */
  )
{
  const Uint NUMBER_OF_SIDES = 6;
  const Uint Nsd[3] = {(Uint)Pgeti(P_"SplitCS_Nsplit_a"),
                       (Uint)Pgeti(P_"SplitCS_Nsplit_b"),
                       (Uint)Pgeti(P_"SplitCS_Nsplit_c")};
  char parU[STR_SIZE3] = {'\0'};
  char parD[STR_SIZE3] = {'\0'};
  char par[STR_SIZE3]  = {'\0'};
  char obj[STR_SIZE1]  = {'\0'};
  char name[STR_SIZE3] = {'\0'};
  const char *dir      = 0;
  Flag_T type = UNDEFINED;
  Uint p,d0,d1,d2,ijk;
  
  /* object name */
  set_object_name_split_CS(obj,obj0);
  
  if (!strcmp(obj,"outermost"))
  {
    if (dir0 == NONE)/* for single objects */
      dir = "NA";
    else
      dir = StrSide[dir0];
    type = OT_T_SCS;
  }
  else
  {
    dir = StrSide[dir0];
    type = OB_T_SCS;
  }
  
  assert(dir0 == LEFT || dir0 == RIGHT || dir0 == CENTER || dir0 == NONE);
  
  for (d0 = 0; d0 < Nsd[0]; d0++)
  {
    for (d1 = 0; d1 <  Nsd[1]; d1++)
    {
      for (d2 = 0; d2 <  Nsd[2]; d2++)
      {
        for (p = 0; p < NUMBER_OF_SIDES; ++p)
        {
          Flag_T side = (Flag_T)(p);
          
          /* left outermost hemisphere, doesn't have right */
          if (dir0 == LEFT && type == OT_T_SCS && side == RIGHT)
            continue;
          
          /* right outermost hemisphere, doesn't have left */
          if (dir0 == RIGHT && type == OT_T_SCS && side == LEFT)
            continue;
          
          Patch_T *const patch = calloc(1,sizeof(*patch));
          IsNull(patch);
          grid->patch = 
            realloc(grid->patch,(grid->np+2)*sizeof(*grid->patch));
          IsNull(grid->patch);
          grid->patch[grid->np]   = patch;
          grid->patch[grid->np+1] = 0;
          
          /* filling grid */
          patch->grid = grid;
        
          /* filling patch number */
          patch->pn = grid->np;
        
          grid->np += 1;
          
          Field_T *R1 = add_field(SigmaD,0,patch,NO);
          Field_T *R2 = add_field(SigmaU,0,patch,NO);
          double *rU = 0, *rD = 0;
          
          assert(StrSide[side]);
          
          /* covering region */
          if (strcmp_i(obj,"NS")         || 
              strcmp_i(obj,"NS1")        ||
              strcmp_i(obj,"NS2")        || 
              strcmp_i(obj,"BH")         ||
              strcmp_i(obj,"BH1")        || 
              strcmp_i(obj,"BH2")        ||
              strcmp_i(obj,"NS_around")  ||
              strcmp_i(obj,"NS1_around") ||
              strcmp_i(obj,"NS2_around") ||
              strcmp_i(obj,"BH_around")  ||
              strcmp_i(obj,"BH1_around") ||
              strcmp_i(obj,"BH2_around") ||
              strcmp_i(obj,"outermost")
             )
          {
            
            /* if only one split => on both Outer and Inner Boundary */
            if (Nsd[2] == 1)
              sprintf(patch->CoordSysInfo->region,
                "(%s)(%s_OB)(%s_IB)",obj,obj,obj);
                
            else if (d2 == Nsd[2]-1)/* if on Outer Boundary */
              sprintf(patch->CoordSysInfo->region,
                "(%s)(%s_OB)",obj,obj);
            
            else if (d2 == 0)/* if on Inner Boundary */
              sprintf(patch->CoordSysInfo->region,
               "(%s)(%s_IB)",obj,obj);
            
            else
              sprintf(patch->CoordSysInfo->region,
                "(%s)",obj);
          }
          else
          {
            Error0(NO_OPTION);
          }
          
          // printf("patch region = %s\n",patch->CoordSysInfo->region);
          
          /* filling flags */
          patch->CoordSysInfo->CubedSphericalCoord->side = side;
          patch->CoordSysInfo->CubedSphericalCoord->type = type;
          
          /* filling n */
          patch->n[0] = (Uint)Pgeti(P_"SplitCS_n_a");
          patch->n[1] = (Uint)Pgeti(P_"SplitCS_n_b");
          patch->n[2] = (Uint)Pgeti(P_"SplitCS_n_c");
          
          /* filling nn */
          patch->nn = total_nodes_patch(patch);
          
          /* filling number of split */
          patch->nsplit[0] = Nsd[0];
          patch->nsplit[1] = Nsd[1];
          patch->nsplit[2] = Nsd[2];

          /* filling name */
          SCS_par_name(name);
          patch->name = dup_s(name);
          
          /* filling Rs */
          SCS_par_sigma(parU,SigmaU);
          SCS_par_sigma(parD,SigmaD);
          rU = Pgetdd(parU);
          rD = Pgetdd(parD);
          patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
          patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
          R1->v = alloc_double(patch->nn);
          R2->v = alloc_double(patch->nn);
          for (ijk = 0; ijk < patch->nn; ++ijk)
          {
            R1->v[ijk] = rD[ijk];
            R2->v[ijk] = rU[ijk];
          }
          
          /* set xc's */
          SCS_par_xc_length(par,"xc1");
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(par);
          
          SCS_par_xc_length(par,"xc2");
          patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(par);
          
          /* filling center */
          SCS_par_CS_center(par,"a");
          patch->c[0] = Pgetd(par);
          SCS_par_CS_center(par,"b");
          patch->c[1] = Pgetd(par);
          SCS_par_CS_center(par,"c");
          patch->c[2] = Pgetd(par);
          
           /* filling min */
          SCS_par_min(par,0);
          patch->min[0] = Pgetd(par);
          SCS_par_min(par,1);
          patch->min[1] = Pgetd(par);
          SCS_par_min(par,2);
          patch->min[2] = Pgetd(par);
          
          /* filling max */
          SCS_par_max(par,0);
          patch->max[0] = Pgetd(par);
          SCS_par_max(par,1);
          patch->max[1] = Pgetd(par);
          SCS_par_max(par,2);
          patch->max[2] = Pgetd(par);

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
        }/* for (p = 0; p < NUMBER_OF_SIDES; ++p) */
      }
    }
  }
}

/* making value of coords. it is a general function for cubed spherical type */
void make_nodes_CubedSpherical_coord(Patch_T *const patch)
{
  struct Collocation_s coll_s[3] = {0};
  const Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
  const Flag_T type = patch->CoordSysInfo->CubedSphericalCoord->type;
  double S = 0; /* sign */
  Uint a = 0, b = 0, c = 0;/* permuted indices */
  const Uint nn = patch->nn;
  const Uint *const n = patch->n;
  const Field_T *const R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f,
                *const R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
  const double xc1 = patch->CoordSysInfo->CubedSphericalCoord->xc1,
               xc2 = patch->CoordSysInfo->CubedSphericalCoord->xc2,
                R1 = patch->CoordSysInfo->CubedSphericalCoord->R1,
                R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
               
  const double *const C = patch->c;/* center of origine translated */
  Uint i,j,k,l;
  
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
        
        ijk_to_i_j_k(l,n,&i,&j,&k);
        X[0] = point_value(i,&coll_s[0]);
        X[1] = point_value(j,&coll_s[1]);
        X[2] = point_value(k,&coll_s[2]);
        d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
        patch->node[l]->X = X;
        x1 = xc1;
        x2 = S*R2_f->v[i_j_k_to_ijk(n,i,j,0)]/d;
        
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
        
        ijk_to_i_j_k(l,n,&i,&j,&k);
        X[0] = point_value(i,&coll_s[0]);
        X[1] = point_value(j,&coll_s[1]);
        X[2] = point_value(k,&coll_s[2]);
        d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
        patch->node[l]->X = X;
        x2 = xc2;
        x1 = S*R1_f->v[i_j_k_to_ijk(n,i,j,0)]/d;
        
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
        double x1,ratio,d;
        
        ijk_to_i_j_k(l,n,&i,&j,&k);
        X[0] = point_value(i,&coll_s[0]);
        X[1] = point_value(j,&coll_s[1]);
        X[2] = point_value(k,&coll_s[2]);
        d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
        patch->node[l]->X = X;
        x1 = xc1;
        ratio  = 1.-S*d*xc1/R2;
        
        x[c] = x1/(1.-ratio*X[2]);
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
        double x1,ratio,d;
        
        ijk_to_i_j_k(l,n,&i,&j,&k);
        X[0] = point_value(i,&coll_s[0]);
        X[1] = point_value(j,&coll_s[1]);
        X[2] = point_value(k,&coll_s[2]);
        d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
        patch->node[l]->X = X;
        
        x1 = S*R1/d;
        ratio  = 1.-R1/R2;
        
        x[c] = x1/(1.-ratio*X[2]);
        x[a] = S*x[c]*X[0];
        x[b] = S*x[c]*X[1];
        
        x[a]+= C[a];
        x[b]+= C[b];
        x[c]+= C[c];
      }
    break;
    case OB_T_SCS:
      for (l = 0; l < nn; l++)
      {
        double *X = alloc_double(3);
        double *x = patch->node[l]->x;
        double x1,x2,d;
        
        ijk_to_i_j_k(l,n,&i,&j,&k);
        X[0] = point_value(i,&coll_s[0]);
        X[1] = point_value(j,&coll_s[1]);
        X[2] = point_value(k,&coll_s[2]);
        d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
        patch->node[l]->X = X;
        
        x1 = S*(xc1 == DBL_MAX ? R1_f->v[i_j_k_to_ijk(n,i,j,0)]/d : xc1);
        x2 = S*(xc2 == DBL_MAX ? R2_f->v[i_j_k_to_ijk(n,i,j,0)]/d : xc2);
        
        x[c] = x1+(x2-x1)*X[2];
        x[a] = S*x[c]*X[0];
        x[b] = S*x[c]*X[1];
        
        x[a]+= C[a];
        x[b]+= C[b];
        x[c]+= C[c];
      }
      break;
    case OT_T_SCS:
      for (l = 0; l < nn; l++)
      {
        double *X = alloc_double(3);
        double *x = patch->node[l]->x;
        double x1,ratio,d;
        
        ijk_to_i_j_k(l,n,&i,&j,&k);
        X[0] = point_value(i,&coll_s[0]);
        X[1] = point_value(j,&coll_s[1]);
        X[2] = point_value(k,&coll_s[2]);
        d = sqrt(1+Pow2(X[0])+Pow2(X[1]));
        patch->node[l]->X = X;

        x1 = S*(xc1 == DBL_MAX ? R1_f->v[i_j_k_to_ijk(n,i,j,0)]/d : xc1);
        ratio  = 1.-R1_f->v[i_j_k_to_ijk(n,i,j,0)]/(R2_f->v[i_j_k_to_ijk(n,i,j,0)]);
        
        x[c] = x1/(1.-ratio*X[2]);
        x[a] = S*x[c]*X[0];
        x[b] = S*x[c]*X[1];
        
        x[a]+= C[a];
        x[b]+= C[b];
        x[c]+= C[c];
      }
      break;
    default:
      Error0(NO_OPTION);
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
        Error0(NO_JOB);
    }/* end of switch */
    R2_derivative(patch);/* surface function derivative */
    make_coeffs_2d(patch->CoordSysInfo->CubedSphericalCoord->R2_f,0,1);
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
        Error0(NO_JOB);
    }/* end of switch */
    R1_derivative(patch);/* surface function derivative */
    make_coeffs_2d(patch->CoordSysInfo->CubedSphericalCoord->R1_f,0,1);
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
        Error0(NO_JOB);
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
        Error0(NO_JOB);
    }/* end of switch */
  }/* end of else if (type == OT_T2_CS) */
  else if (type == OB_T_SCS)
  {
    /* transformation function */
    patch->JacobianT->j  = JT_OJ_T_SCS;
    
    /* set sign and permutations */
    double sign = 0;
    Uint iper = 0,jper = 0,kper = 0;
    SignAndIndex_permutation_CubedSphere(side,&iper,&jper,&kper,&sign);
    patch->JacobianT->SCS->sign = sign;
    patch->JacobianT->SCS->iper = iper;
    patch->JacobianT->SCS->jper = jper;
    patch->JacobianT->SCS->kper = kper;
    
    /* surface function derivative */
    R12_derivatives_SCS(patch);
    make_coeffs_2d(patch->CoordSysInfo->CubedSphericalCoord->R1_f,0,1);
    make_coeffs_2d(patch->CoordSysInfo->CubedSphericalCoord->R2_f,0,1);
    
  }/* end of if (type == OB_T_SCS) */
  else if (type == OT_T_SCS)
  {
    /* transformation function */
    patch->JacobianT->j  = JT_OT_T_SCS;
    
    /* set sign and permutations */
    double sign = 0;
    Uint iper = 0,jper = 0,kper = 0;
    SignAndIndex_permutation_CubedSphere(side,&iper,&jper,&kper,&sign);
    patch->JacobianT->SCS->sign = sign;
    patch->JacobianT->SCS->iper = iper;
    patch->JacobianT->SCS->jper = jper;
    patch->JacobianT->SCS->kper = kper;
    
    /* surface function derivative */
    R12_derivatives_SCS(patch);
    make_coeffs_2d(patch->CoordSysInfo->CubedSphericalCoord->R1_f,0,1);
    make_coeffs_2d(patch->CoordSysInfo->CubedSphericalCoord->R2_f,0,1);
  }/* end of if (type == OT_T_SCS) */
  else
    Error0(NO_OPTION);
  
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
  Field_T *const R1 = patch->fields[Ind("surface_function")];
  const Uint nn = patch->nn;
  Uint ijk;
          
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
  Field_T *const R2 = patch->fields[Ind("surface_function")];
  const Uint nn = patch->nn;
  Uint ijk;
          
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

/* preparing R1 and R2 derivatives (surface functions ) 
// of Split Cubed Spherical coords.
// NOTE: one must remove dR_? after each updating of R. */
static void R12_derivatives_SCS(Patch_T *const patch)
{
  Field_T *dR2_dX = add_field("dR2_dX",0,patch,NO),
          *dR2_dY = add_field("dR2_dY",0,patch,NO),
          *dR2_dx = add_field("dR2_dx",0,patch,YES),
          *dR2_dy = add_field("dR2_dy",0,patch,YES),
          *dR2_dz = add_field("dR2_dz",0,patch,YES);
  Field_T *const R2 = patch->fields[Ind(SigmaU)];
  
  Field_T *dR1_dX = add_field("dR1_dX",0,patch,NO),
          *dR1_dY = add_field("dR1_dY",0,patch,NO),
          *dR1_dx = add_field("dR1_dx",0,patch,YES),
          *dR1_dy = add_field("dR1_dy",0,patch,YES),
          *dR1_dz = add_field("dR1_dz",0,patch,YES);
  Field_T *const R1 = patch->fields[Ind(SigmaD)];
  
  const Uint nn = patch->nn;
  Uint ijk;
          
  dR2_dX->v = Partial_Derivative(R2,"a");
  dR2_dY->v = Partial_Derivative(R2,"b");
  
  dR1_dX->v = Partial_Derivative(R1,"a");
  dR1_dY->v = Partial_Derivative(R1,"b");
    
  for (ijk = 0; ijk < nn; ++ijk)
  {
    dR2_dx->v[ijk] = dR2_dX->v[ijk]*dq2_dq1(patch,_a_,_x_,ijk)+
                     dR2_dY->v[ijk]*dq2_dq1(patch,_b_,_x_,ijk);
    dR2_dy->v[ijk] = dR2_dX->v[ijk]*dq2_dq1(patch,_a_,_y_,ijk)+
                     dR2_dY->v[ijk]*dq2_dq1(patch,_b_,_y_,ijk);
    dR2_dz->v[ijk] = dR2_dX->v[ijk]*dq2_dq1(patch,_a_,_z_,ijk)+
                     dR2_dY->v[ijk]*dq2_dq1(patch,_b_,_z_,ijk);
                     
    dR1_dx->v[ijk] = dR1_dX->v[ijk]*dq2_dq1(patch,_a_,_x_,ijk)+
                     dR1_dY->v[ijk]*dq2_dq1(patch,_b_,_x_,ijk);
    dR1_dy->v[ijk] = dR1_dX->v[ijk]*dq2_dq1(patch,_a_,_y_,ijk)+
                     dR1_dY->v[ijk]*dq2_dq1(patch,_b_,_y_,ijk);
    dR1_dz->v[ijk] = dR1_dX->v[ijk]*dq2_dq1(patch,_a_,_z_,ijk)+
                     dR1_dY->v[ijk]*dq2_dq1(patch,_b_,_z_,ijk);
  }
                      
  remove_field(dR2_dX);
  remove_field(dR2_dY);
  
  remove_field(dR1_dX);
  remove_field(dR1_dY);
  
  patch->CoordSysInfo->CubedSphericalCoord->dR2_dx = dR2_dx;
  patch->CoordSysInfo->CubedSphericalCoord->dR2_dy = dR2_dy;
  patch->CoordSysInfo->CubedSphericalCoord->dR2_dz = dR2_dz;
  
  patch->CoordSysInfo->CubedSphericalCoord->dR1_dx = dR1_dx;
  patch->CoordSysInfo->CubedSphericalCoord->dR1_dy = dR1_dy;
  patch->CoordSysInfo->CubedSphericalCoord->dR1_dz = dR1_dz;
  
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
static void populate_left_NS(Grid_T *const grid,const Uint pn)
{
  Uint p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R2 = add_field("surface_function",0,patch,NO);
    double *R2_array;
    Flag_T side = (Flag_T)(p-pn);
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    Uint n,ijk;
    
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
    patch->n[0] = (Uint)PgetiEZ("n_a");
    patch->n[1] = (Uint)PgetiEZ("n_b");
    patch->n[2] = (Uint)PgetiEZ("n_c");
    /* check for override */
    sprintf(var,"left_NS");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      Error0("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      Error0("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      Error0("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_up",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"left_central_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_down",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"left_central_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_left",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"left_central_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];

      break;
      case RIGHT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_right",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"left_central_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_back",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"left_central_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_front",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"left_central_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      default:
        Error0(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_a",grid->gn);
    patch->c[0] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_b",grid->gn);
    patch->c[1] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_c",grid->gn);
    patch->c[2] = Pgetd(var);
    
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

/* populating properties of patch for right BH */
void populate_right_BH(Grid_T *const grid,const Uint pn)
{
  Uint p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R2 = add_field("surface_function",0,patch,NO);
    double *R2_array;
    Flag_T side = (Flag_T)(p-pn);
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    Uint n,ijk;
    
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
    patch->n[0] = (Uint)PgetiEZ("n_a");
    patch->n[1] = (Uint)PgetiEZ("n_b");
    patch->n[2] = (Uint)PgetiEZ("n_c");
    /* check for override */
    sprintf(var,"right_BH");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      Error0("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      Error0("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      Error0("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_up",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_down",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_left",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];

      break;
      case RIGHT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_right",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_back",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_front",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      default:
        Error0(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,PATCH_NAME_PRT_P_"right_BH_center_a",grid->gn);
    patch->c[0] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"right_BH_center_b",grid->gn);
    patch->c[1] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"right_BH_center_c",grid->gn);
    patch->c[2] = Pgetd(var);
    
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

/* populating properties of patch for central NS */
static void populate_central_NS(Grid_T *const grid,const Uint pn)
{
  Uint p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R2 = add_field("surface_function",0,patch,NO);
    double *R2_array;
    Flag_T side = (Flag_T)(p-pn);
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    Uint n,ijk;
    
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
    patch->n[0] = (Uint)PgetiEZ("n_a");
    patch->n[1] = (Uint)PgetiEZ("n_b");
    patch->n[2] = (Uint)PgetiEZ("n_c");
    /* check for override */
    sprintf(var,"NS");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      Error0("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      Error0("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      Error0("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_up",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"central_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_down",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"central_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_left",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"central_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];

      break;
      case RIGHT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_right",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"central_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_back",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"central_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_front",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"central_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      default:
        Error0(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,PATCH_NAME_PRT_P_"NS_center_a",grid->gn);
    patch->c[0] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"NS_center_b",grid->gn);
    patch->c[1] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"NS_center_c",grid->gn);
    patch->c[2] = Pgetd(var);
    
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
static void populate_right_NS(Grid_T *const grid,const Uint pn)
{
  Uint p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R2 = add_field("surface_function",0,patch,NO);
    double *R2_array;
    Flag_T side = (Flag_T)(p-pn);
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    Uint n,ijk;
    
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
    patch->n[0] = (Uint)PgetiEZ("n_a");
    patch->n[1] = (Uint)PgetiEZ("n_b");
    patch->n[2] = (Uint)PgetiEZ("n_c");
    /* check for override */
    sprintf(var,"right_NS");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      Error0("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      Error0("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      Error0("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);

    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_up",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_down",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_c",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_left",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_right",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_b",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_back",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_front",grid->gn);
        R2_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R2_f = R2;
        sprintf(var,PATCH_NAME_PRT_P_"right_central_box_size_a",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var)/2.;
        R2->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R2->v[ijk] = R2_array[ijk];
      
      break;
      default:
        Error0(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_a",grid->gn);
    patch->c[0] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_b",grid->gn);
    patch->c[1] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_c",grid->gn);
    patch->c[2] = Pgetd(var);
    
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

/* populating properties of patch for right NS around */
static void populate_right_NS_around(Grid_T *const grid,const Uint pn)
{
  Uint p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R1 = add_field("surface_function",0,patch,NO);
    double *R1_array;
    Flag_T side = (Flag_T)(p-pn);
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    Uint n,ijk;
    
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
    patch->n[0] = (Uint)PgetiEZ("n_a");
    patch->n[1] = (Uint)PgetiEZ("n_b");
    patch->n[2] = (Uint)PgetiEZ("n_c");
    /* check for override */
    sprintf(var,"right_NS");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      Error0("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      Error0("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      Error0("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_around_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_up",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];

      break;
      case DOWN:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_around_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_down",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_around_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_left",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_around_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_right",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_around_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_back",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_NS_around_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_NS_surface_function_front",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      default:
        Error0(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_a",grid->gn);
    patch->c[0] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_b",grid->gn);
    patch->c[1] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"right_NS_center_c",grid->gn);
    patch->c[2] = Pgetd(var);
    
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

/* populating properties of patch for right BH around */
static void populate_right_BH_around(Grid_T *const grid,const Uint pn)
{
  Uint p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R1 = add_field("surface_function",0,patch,NO);
    double *R1_array;
    Flag_T side = (Flag_T)(p-pn);
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    Uint n,ijk;
    
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
    patch->n[0] = (Uint)PgetiEZ("n_a");
    patch->n[1] = (Uint)PgetiEZ("n_b");
    patch->n[2] = (Uint)PgetiEZ("n_c");
    /* check for override */
    sprintf(var,"right_BH");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      Error0("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      Error0("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      Error0("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_around_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_up",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_around_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_down",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_around_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_left",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_around_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_right",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_around_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_back",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"right_BH_around_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"right_BH_surface_function_front",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      default:
        Error0(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,PATCH_NAME_PRT_P_"right_BH_center_a",grid->gn);
    patch->c[0] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"right_BH_center_b",grid->gn);
    patch->c[1] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"right_BH_center_c",grid->gn);
    patch->c[2] = Pgetd(var);
    
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

/* populating properties of patch for left NS around */
static void populate_left_NS_around(Grid_T *const grid,const Uint pn)
{
  Uint p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R1 = add_field("surface_function",0,patch,NO);
    double *R1_array;
    Flag_T side = (Flag_T)(p-pn);
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    Uint n,ijk;
    
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
    patch->n[0] = (Uint)PgetiEZ("n_a");
    patch->n[1] = (Uint)PgetiEZ("n_b");
    patch->n[2] = (Uint)PgetiEZ("n_c");
    /* check for override */
    sprintf(var,"left_NS");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      Error0("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      Error0("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      Error0("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_around_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_up",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_around_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_down",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_around_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_left",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_around_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_right",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_around_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_back",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"left_NS_around_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"left_NS_surface_function_front",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var)/2.;
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      default:
        Error0(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_a",grid->gn);
    patch->c[0] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_b",grid->gn);
    patch->c[1] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"left_NS_center_c",grid->gn);
    patch->c[2] = Pgetd(var);
    
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

/* populating properties of patch for central BH around */
static void populate_central_BH_around(Grid_T *const grid,const Uint pn)
{
  Uint p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R1 = add_field("surface_function",0,patch,NO);
    double *R1_array;
    Flag_T side = (Flag_T)(p-pn);
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    Uint n,ijk;
    
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
    patch->n[0] = (Uint)PgetiEZ("n_a");
    patch->n[1] = (Uint)PgetiEZ("n_b");
    patch->n[2] = (Uint)PgetiEZ("n_c");
    /* check for override */
    sprintf(var,"BH");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      Error0("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      Error0("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      Error0("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"BH_around_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"BH_surface_function_up",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"BH_around_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"BH_surface_function_down",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"BH_around_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"BH_surface_function_left",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"BH_around_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"BH_surface_function_right",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"BH_around_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"BH_surface_function_back",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"BH_around_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"BH_surface_function_front",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      default:
        Error0(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,PATCH_NAME_PRT_P_"BH_center_a",grid->gn);
    patch->c[0] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"BH_center_b",grid->gn);
    patch->c[1] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"BH_center_c",grid->gn);
    patch->c[2] = Pgetd(var);
    
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

/* populating properties of patch for central NS around */
static void populate_central_NS_around(Grid_T *const grid,const Uint pn)
{
  Uint p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Field_T *R1 = add_field("surface_function",0,patch,NO);
    double *R1_array;
    Flag_T side = (Flag_T)(p-pn);
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    Uint n,ijk;
    
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
    patch->n[0] = (Uint)PgetiEZ("n_a");
    patch->n[1] = (Uint)PgetiEZ("n_b");
    patch->n[2] = (Uint)PgetiEZ("n_c");
    /* check for override */
    sprintf(var,"NS");
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      Error0("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      Error0("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      Error0("n_c could not be set.\n");
    
    /* filling nn */
    patch->nn = total_nodes_patch(patch);
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_around_up",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_up",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case DOWN:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_around_down",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_down",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_around_left",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_left",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_around_right",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_right",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_around_back",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_back",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = -Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      case FRONT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"NS_around_front",grid->gn);
        patch->name = dup_s(name);
        
        /* filling Rs */
        sprintf(var,PATCH_NAME_PRT_P_"NS_surface_function_front",grid->gn);
        R1_array = Pgetdd(var);
        patch->CoordSysInfo->CubedSphericalCoord->R1_f = R1;
        sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
        patch->CoordSysInfo->CubedSphericalCoord->xc2 = Pgetd(var);
        R1->v = alloc_double(patch->nn);
        
        for (ijk = 0; ijk < patch->nn; ++ijk)
          R1->v[ijk] = R1_array[ijk];
      
      break;
      default:
        Error0(NO_OPTION);
    }
    
    /* filling center */
    sprintf(var,PATCH_NAME_PRT_P_"NS_center_a",grid->gn);
    patch->c[0] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"NS_center_b",grid->gn);
    patch->c[1] = Pgetd(var);
    sprintf(var,PATCH_NAME_PRT_P_"NS_center_c",grid->gn);
    patch->c[2] = Pgetd(var);
    
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
static void populate_outermost(Grid_T *const grid,const Uint pn,const Uint o)
{
  Uint p;/* patch */
  
  for (p = pn; p < pn+6; p++)
  {
    Patch_T *const patch = grid->patch[p];
    Flag_T side = (Flag_T)(p-pn);
    char name[100] = {'\0'};
    char var[100] = {'\0'};
    struct Ret_S ret;
    Uint n;
    
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
          sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var); 
        break;
        case DOWN:
          sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var); 
        break;
        case LEFT:
          sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var);         
        break;
        case RIGHT:
          sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var);         
        break;
        case BACK:
          sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = -Pgetd(var);           
        break;
        case FRONT:
          sprintf(var,PATCH_NAME_PRT_P_"around_box_length",grid->gn);
          patch->CoordSysInfo->CubedSphericalCoord->xc1 = Pgetd(var);         
        break;
        default:
          Error0(NO_OPTION);
      }

      /* filling Rs */
      sprintf(var,PATCH_NAME_PRT_P_"outermost%u_R2",grid->gn,o);
      patch->CoordSysInfo->CubedSphericalCoord->R2 = Pgetd(var);
    
    }
    else
    {
      patch->CoordSysInfo->CubedSphericalCoord->type = OT_T2_CS;
      /* filling Rs */
      sprintf(var,PATCH_NAME_PRT_P_"outermost%u_R1",grid->gn,o);
      patch->CoordSysInfo->CubedSphericalCoord->R1 = Pgetd(var);
      sprintf(var,PATCH_NAME_PRT_P_"outermost%u_R2",grid->gn,o);
      patch->CoordSysInfo->CubedSphericalCoord->R2 = Pgetd(var);
    
    }
    
    switch(side)
    {
      case UP:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"outermost%u_up",grid->gn,o);
        patch->name = dup_s(name);
        
      break;
      case DOWN:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"outermost%u_down",grid->gn,o);
        patch->name = dup_s(name);
      
      break;
      case LEFT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"outermost%u_left",grid->gn,o);
        patch->name = dup_s(name);
        
      
      break;
      case RIGHT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"outermost%u_right",grid->gn,o);
        patch->name = dup_s(name);
        
      
      break;
      case BACK:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"outermost%u_back",grid->gn,o);
        patch->name = dup_s(name);
        
      break;
      case FRONT:
        /* filling name */
        sprintf(name,PATCH_NAME_PRT_P_"outermost%u_front",grid->gn,o);
        patch->name = dup_s(name);
      
      break;
      default:
        Error0(NO_OPTION);
    }

    /* filling n */
    patch->n[0] = (Uint)PgetiEZ("n_a");
    patch->n[1] = (Uint)PgetiEZ("n_b");
    patch->n[2] = (Uint)PgetiEZ("n_c");
    /* check for override */
    sprintf(var,"Outermost%u",o);
    make_keyword_parameter(&ret,var,"n");
    n = (Uint)PgetiEZ(ret.s0);
    if (n != INT_MAX)	patch->n[0] = n;
    n = (Uint)PgetiEZ(ret.s1);
    if (n != INT_MAX)	patch->n[1] = n;
    n = (Uint)PgetiEZ(ret.s2);
    if (n != INT_MAX)	patch->n[2] = n;
      
    if(patch->n[0] == INT_MAX)
      Error0("n_a could not be set.\n");
    if(patch->n[1] == INT_MAX)
      Error0("n_b could not be set.\n");
    if(patch->n[2] == INT_MAX)
      Error0("n_c could not be set.\n");
    
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
  Uint Np = 30;/* number of patches without outermost's 
                   4 sets of cubed sphere = 4*6
                   4 filling box
                   2 central box */
  Uint outermost;
  Uint i;
  
  outermost = (Uint) PgetiEZ("Number_of_Outermost_Split");
  if (outermost != (Uint)INT_MAX)
    Np += 6*outermost;
  
  grid->patch = calloc((Np+1),sizeof(*grid->patch));
  IsNull(grid->patch);
  
  for (i = 0; i < Np; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    IsNull(grid->patch[i]);
  }
  
}

/* memory alloc patches for BBN_CubedSpherical type */
void alloc_patches_BBN_CubedSpherical_grid(Grid_T *const grid)
{
  Uint Np = 23;/* number of patches without outermost's 
                   3 sets of cubed sphere = 3*6
                   4 filling boxex
                   1 central box */
  Uint outermost;
  Uint i;
  
  outermost = (Uint) PgetiEZ("Number_of_Outermost_Split");
  if (outermost != (Uint)INT_MAX)
    Np += 6*outermost;
  
  grid->patch = calloc((Np+1),sizeof(*grid->patch));
  IsNull(grid->patch);
  
  for (i = 0; i < Np; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    IsNull(grid->patch[i]);
  }
  
}

/* given type of the object (lower case or capital) if exists
// it writes the object name correctly, otherwise gives error. */
void set_object_name_split_CS(char *const obj,const char *const type)
{
  assert(obj);
  
  Uint i = 0;
  while (SCS_ObjType[i])
  {
    if (strcmp_i(SCS_ObjType[i],type))
    {
      sprintf(obj,"%s",SCS_ObjType[i]);
      return;
    }
    i++;
  }
  
  Error0("No such object name!\n");
}

/* set parameters of split cubed spherical, number of splits,
// surface functions etc.
// this must be called in characteristic function.
// NOTE: this also sets parameters of those regions
// which are not used in grid, (for example, BH patches which scooped out
// from grid); thus, one can later use these parameters if wants to
// cover those regions for instance in BH filler. 
// NOTE: separate parts of this function have written in braces 
// thus in future one can use this components for similar purposes. */
void set_params_of_split_cubed_spherical_grid(Grid_Char_T *const grid_char)
{
  FUNC_TIC
  
  Grid_T *const grid = grid_char->grid;/* this is the new grid */
  const int Verbose  = Pcmps("grid_verbose","yes");
  Uint n[3] = {0};
  Uint i,j,k,d0,d1,d2;
  
  assert(grid);
  
  /* first the kind of grid */
  grid->kind = set_grid_kind(Pgets("grid_kind"));
  
  /* find out the splits and resolutions */
  /* { */
  /* resolution in each dir. */
  const Uint givenN[3] = {(Uint)Pgeti("n_a"),
                          (Uint)Pgeti("n_b"),
                          (Uint)Pgeti("n_c")};
  const Uint maxN[3] = {(Uint)Pgeti("grid_SplitCS_max_n_a"),
                        (Uint)Pgeti("grid_SplitCS_max_n_b"),
                        (Uint)Pgeti("grid_SplitCS_max_n_c")};
  Uint newN[3] = {0};/* new number of nodes in each direction */
  Uint Nsd[3]  = {0};/* number of splits in each dir. */
  Uint Nns[3]  = {0};/* number of nodes in each split patch */
  const double r_outermost = Pgetd("grid_outermost_radius");
  
  assert(maxN[0]>2 && maxN[1]>2 && maxN[2]>2);
  
  /* for each direction divide until the max resolution */
  for (i = 0; i < 3; ++i)
  {
    Uint d = 1;/* divide */
    n[i] = givenN[i]/d;
    while (n[i] > maxN[i])
    {
      d++;
      n[i] = givenN[i]/d;
      assert(d < givenN[i]);
    }
    Nsd[i] = d;
    /* adjust givenN */
    newN[i] = d*n[i];
    Nns[i]  = n[i];
  }
  /* adjust parameters */
  Pseti(P_"SplitCS_n_a",(int)Nns[0]);
  Pseti(P_"SplitCS_n_b",(int)Nns[1]);
  Pseti(P_"SplitCS_n_c",(int)Nns[2]);
  /* number of splits in each direction */
  Pseti(P_"SplitCS_Nsplit_a",(int)Nsd[0]);
  Pseti(P_"SplitCS_Nsplit_b",(int)Nsd[1]);
  Pseti(P_"SplitCS_Nsplit_c",(int)Nsd[2]);
  /* current number of splits in each direction (global parameter) */
  Pseti("grid_SplitCS_Nsplit_a",(int)Nsd[0]);
  Pseti("grid_SplitCS_Nsplit_b",(int)Nsd[1]);
  Pseti("grid_SplitCS_Nsplit_c",(int)Nsd[2]);
  
  /* test */
  if(Verbose)
  printf(Pretty0"given params (n_a,n_b,n_c) = (%u,%u,%u)\n"
         Pretty0"modified     (n_a,n_b,n_c) = (%u,%u,%u)\n"
         Pretty0"max allowed  (n_a,n_b,n_c) = (%u,%u,%u)\n"
         Pretty0"each patch   (n_a,n_b,n_c) = (%u,%u,%u)\n"
         Pretty0"split        (  a,  b,  c) = (%u,%u,%u)\n",
         givenN[0],givenN[1],givenN[2],
         newN[0],newN[1],newN[2],
         maxN[0],maxN[1],maxN[2],
         Nns[0],Nns[1],Nns[2],
         Nsd[0],Nsd[1],Nsd[2]);
  /* } */
  
  /* set parameters for all patches for Binary system */
  if (Pcmps("grid_kind","SplitCubedSpherical(BH+NS)") ||
      Pcmps("grid_kind","SplitCubedSpherical(NS+NS)") ||
      Pcmps("grid_kind","SplitCubedSpherical(BH+BH)"))
  
  { 
  const Uint Num_Obj = 2;
  const double S  = grid_char->S;/* separation */
  Uint obj_n;/* BH or NS */
  
  /* two different directions and lowercase */
  if (strstr_i(grid_char->params[0]->dir,"left"))
    grid_char->params[0]->dir = "left";
  else if (strstr_i(grid_char->params[0]->dir,"right"))
    grid_char->params[0]->dir = "right";
  else
    Error0("Bad argument.");
    
  if (strstr_i(grid_char->params[1]->dir,"left"))
    grid_char->params[1]->dir = "left";
  else if (strstr_i(grid_char->params[1]->dir,"right"))
    grid_char->params[1]->dir = "right";
  else
    Error0("Bad argument.");
    
  assert(!strstr(grid_char->params[0]->dir,grid_char->params[1]->dir));
  
  if(S < 0)
    Error0("The distance between the two compact objects "
           "must be positive.\n");
  /* first populate parameters only for objects NS/BH */
  for (obj_n = 0; obj_n < Num_Obj; ++obj_n)
  {
    /* note (X,Y,Z) in [-1,1]x[-1,1]x[0,1]. */
    const double Xm = -1,XM = 1;
    const double Ym = -1,YM = 1;
    const double Zm = 0 ,ZM = 1;
    /* step in each direction, note X in [-1,1]x[-1,1]x[0,1]. 
    // also NOTE that Z min and max are 0 and 1 in cubed spherical. */
    double step[3] = {(XM-Xm)/Nsd[0],(YM-Ym)/Nsd[1],DBL_MAX};
    double min[3] = {0},max[3] = {0};
    double rup,rdown;
    double th = 0,ph = 0,X[3] = {0};
    double *rU = 0, *rD = 0;
    Patch_T patch[1] = {0};
    struct Collocation_s coll_s[2] = {0};
    char parU[STR_SIZE3] = {'\0'};
    char parD[STR_SIZE3] = {'\0'};
    char par[STR_SIZE3]  = {'\0'};
    char obj[STR_SIZE1]  = {'\0'};
    Uint N_total,p;
    
    const char *obj0 = grid_char->params[obj_n]->obj;
    const char *dir  = grid_char->params[obj_n]->dir;
    double l = grid_char->params[obj_n]->l;
    double w = grid_char->params[obj_n]->w;
    double h = grid_char->params[obj_n]->h;
    const double *reClm = grid_char->params[obj_n]->relClm;
    const double *imClm = grid_char->params[obj_n]->imgClm;
    Uint lmax = grid_char->params[obj_n]->lmax;
    /* find r step */
    double diag = sqrt(Pow2(l)+Pow2(w)+Pow2(h))/2.;
    double rmin = diag;
    double rmax = grid_char->params[obj_n]->r_min;
    double rstep = (rmax-rmin)/Nsd[2];
    
    set_object_name_split_CS(obj,obj0);
    
    /* some checks */
    assert(rmax > 0);
    assert(rstep > 0);
    assert(l > 0 && w > 0 && h > 0);
    
    /* must be lower case letters */
    if(strcmp(dir,"left") && strcmp(dir,"right"))
      Error0("Direction must be lower case!");
    
    if(rmax < 0)
      Errors("%s object must have positive radius.\n",dir);
    if(2*rmax > S)
      Errors("%s object radius is too big.\n",dir);

    for (d0 = 0; d0 < Nsd[0]; d0++)
    {
      min[0] = Xm + step[0]*d0;
      max[0] = Xm + step[0]*(d0+1);

      for (d1 = 0; d1 <  Nsd[1]; d1++)
      {
        min[1] = Ym + step[1]*d1;
        max[1] = Ym + step[1]*(d1+1);
        
        for (d2 = 0; d2 <  Nsd[2]; d2++)
        {
          min[2] = Zm;
          max[2] = ZM;
          
          rdown  = rmin + rstep*d2;
          rup    = rmin + rstep*(d2+1);
          
          assert(rup > diag);
          
          /* set xc to default which is DBL_MAX.
          // NOTE: it's important to set it as a double.
          // since i don't know what happens if you set it as 
          // a string par and convert it into double! */
          for (p = 0; p < 6; ++p)
          {
            Flag_T side = (Flag_T)(p);
            SCS_par_xc_length(par,"xc2");
            cs_Psetd(par,DBL_MAX);
            SCS_par_xc_length(par,"xc1");
            cs_Psetd(par,DBL_MAX);
          }
          
          /* set min and max parameters */
          for (p = 0; p < 6; ++p)
          {
            Flag_T side = (Flag_T)(p);
            SCS_par_min(par,0);
            cs_Psetd(par,min[0]);
            SCS_par_min(par,1);
            cs_Psetd(par,min[1]);
            SCS_par_min(par,2);
            cs_Psetd(par,min[2]);
            SCS_par_max(par,0);
            cs_Psetd(par,max[0]);
            SCS_par_max(par,1);
            cs_Psetd(par,max[1]);
            SCS_par_max(par,2);
            cs_Psetd(par,max[2]);
          }
          
          /* filling min */
          patch->min[0] = min[0];
          patch->min[1] = min[1];
          patch->min[2] = min[2];

          /* filling max */
          patch->max[0] = max[0];
          patch->max[1] = max[1];
          patch->max[2] = max[2];

          /* collocation */
          patch->collocation[0] = Chebyshev_Extrema;
          patch->collocation[1] = Chebyshev_Extrema;
          patch->collocation[2] = Chebyshev_Extrema;

          /* basis */
          patch->basis[0] = Chebyshev_Tn_BASIS;
          patch->basis[1] = Chebyshev_Tn_BASIS;
          patch->basis[2] = Chebyshev_Tn_BASIS;
          
          /* n */
          patch->n[0] = Nns[0];
          patch->n[1] = Nns[1];
          patch->n[2] = Nns[2];
          
          initialize_collocation_struct(patch,&coll_s[0],0);
          initialize_collocation_struct(patch,&coll_s[1],1);
          
          N_total = Nns[0]*Nns[1]*Nns[2];
          rU = alloc_double(N_total);
          rD = alloc_double(N_total);
          
          /* note: order matters */
          /* if top level and d2 > 0 */
          if ( d2 != 0 && d2 == Nsd[2]-1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);              
              
              X[2] = 1;
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  theta_phi_of_XY_CS(&th,&ph,X,side);
                  double r = interpolation_Ylm(reClm,imClm,lmax,th,ph);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = r;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                    assert(r>rdown);
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if bottom level and more than 1 split, one has flat surf. */
          else if (d2 == 0 && Nsd[2] > 1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);
              
              X[2] = DBL_MAX;/* catch error */
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set small xc (used in object) as opposed to up
              // (used in around) to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc1");
              cs_Psetd(par,xc);
              
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if only one level => one of them has flat surface */
          else if (d2 == 0 && Nsd[2] == 1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);

              X[2] = 1;
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  theta_phi_of_XY_CS(&th,&ph,X,side);
                  double r = interpolation_Ylm(reClm,imClm,lmax,th,ph);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = r;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set small xc (used in object) as opposed to up 
              // (used in around) to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc1");
              cs_Psetd(par,xc);
              
            }/* for (p = 0; p < 6; ++p) */
          }
          /* between the above cases they have perfect S2 surface */
          else
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              for (i = 0; i < Nns[0]; ++i)
              {
                for (j = 0; j < Nns[1]; ++j)
                {
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          Free(rU);
          Free(rD);
        }/* for (d2 = 0; d2 <  Nsd[2]; d2++) */
      }/* for (d1 = 0; d1 <  Nsd[1]; d1++) */
    }/* for (d0 = 0; d0 < Nsd[0]; d0++) */
    
    /* set center of patch */
    /* assuming objects are on y-axis */
    for (p = 0; p < 6; ++p)
    {
      Flag_T side = (Flag_T)(p);
      SCS_par_CS_center(par,"a");
      cs_Psetd(par,0.0);
    
      SCS_par_CS_center(par,"b");
      if (!strcmp(dir,"left"))
        cs_Psetd(par,-S/2.);
      else if (!strcmp(dir,"right"))
        cs_Psetd(par,S/2.);
      else
        Error0(NO_OPTION);
      
      SCS_par_CS_center(par,"c");
      cs_Psetd(par,0.0);
    }
    
    /* set parameter for centeral box */
    obj0 = "central_box";
    step[0] = l/Nsd[0];
    step[1] = w/Nsd[1];
    step[2] = h/Nsd[2];
    set_object_name_split_CS(obj,obj0);
    
    if (!strcmp(dir,"left"))
    {
      Flag_T side = LEFT;
      double cen[3] = {0};
      /* orign center (0,-S/2.,0) */
      cen[0] = -l/2.+step[0]/2.;
      for (d0 = 0; d0 < Nsd[0]; d0++)
      {
        cen[1] = -S/2.-w/2.+step[1]/2.;
        for (d1 = 0; d1 <  Nsd[1]; d1++)
        {
          cen[2] = -h/2.+step[2]/2.;
          for (d2 = 0; d2 <  Nsd[2]; d2++)
          {
            /* assuming objects are on y-axis */
            SCS_par_box_center(par,"a");
            cs_Psetd(par,cen[0]);
            SCS_par_box_center(par,"b");
            cs_Psetd(par,cen[1]);
            SCS_par_box_center(par,"c");
            cs_Psetd(par,cen[2]);
            
            cen[2] += step[2];
          }
          cen[1] += step[1];
        }
        cen[0] += step[0];
      }
    }
    else if (!strcmp(dir,"right"))
    {
      Flag_T side = RIGHT;
      double cen[3] = {0};
      /* orign center (0,S/2.,0) */
      cen[0] = -l/2.+step[0]/2.;
      for (d0 = 0; d0 < Nsd[0]; d0++)
      {
        cen[1] = S/2.-w/2.+step[1]/2.;
        for (d1 = 0; d1 <  Nsd[1]; d1++)
        {
          cen[2] = -h/2.+step[2]/2.;
          for (d2 = 0; d2 <  Nsd[2]; d2++)
          {
            /* assuming objects are on y-axis */
            SCS_par_box_center(par,"a");
            cs_Psetd(par,cen[0]);
            SCS_par_box_center(par,"b");
            cs_Psetd(par,cen[1]);
            SCS_par_box_center(par,"c");
            cs_Psetd(par,cen[2]);
            
            cen[2] += step[2];
          }
          cen[1] += step[1];
        }
        cen[0] += step[0];
      }
    }
    else
      Error0(NO_OPTION);
    
    /* set lengths */
    for (d0 = 0; d0 < Nsd[0]; d0++)
    {
      for (d1 = 0; d1 <  Nsd[1]; d1++)
      {
        for (d2 = 0; d2 <  Nsd[2]; d2++)
        {
          SCS_par_box_length(par,"l");
          cs_Psetd(par,step[0]);
         
          SCS_par_box_length(par,"w");
          cs_Psetd(par,step[1]);
         
          SCS_par_box_length(par,"h");
          cs_Psetd(par,step[2]);
        }/* for (d2 = 0; d2 <  Nsd[2]; d2++) */
      }/* for (d1 = 0; d1 <  Nsd[1]; d1++) */
    }/* for (d0 = 0; d0 < Nsd[0]; d0++) */
  }
  
  /* populate parameters arounds of objects. */
  for (obj_n = 0; obj_n < Num_Obj; ++obj_n)
  {
    /* note (X,Y,Z) in [-1,1]x[-1,1]x[0,1]. */
    const double Xm = -1,XM = 1;
    const double Ym = -1,YM = 1;
    const double Zm = 0 ,ZM = 1;
    /* step in each direction, note X in [-1,1]x[-1,1]x[0,1]. 
    // also NOTE that Z min and max are 0 and 1 in cubed spherical. */
    const double step[3] = {(XM-Xm)/Nsd[0],(YM-Ym)/Nsd[1],DBL_MAX};
    double min[3] = {0},max[3] = {0};
    double rup,rdown;
    double th = 0,ph = 0,X[3] = {0};
    double *rU = 0, *rD = 0;
    Patch_T patch[1] = {0};
    struct Collocation_s coll_s[2] = {0};
    char parU[STR_SIZE3] = {'\0'};
    char parD[STR_SIZE3] = {'\0'};
    char par[STR_SIZE3]  = {'\0'};
    char obj[STR_SIZE1]  = {'\0'};
    Uint N_total,p;
    
    /* find r step */
    const char *objstem = grid_char->params[obj_n]->obj;
    const char *dir = grid_char->params[obj_n]->dir;
    double l = S;
    double w = S;
    double h = S;
    const double *reClm = grid_char->params[obj_n]->relClm;
    const double *imClm = grid_char->params[obj_n]->imgClm;
    Uint lmax = grid_char->params[obj_n]->lmax;
    double diag = sqrt(Pow2(l)+Pow2(w)+Pow2(h))/2.;
    double rmax = l/2 > MaxMag_d(w/2,h/2) ? l/2 : MaxMag_d(w/2,h/2);
    double rmin = grid_char->params[obj_n]->r_min;
    double rstep = (rmax-rmin)/Nsd[2];
    
    sprintf(obj,"%s_around",objstem);
    set_object_name_split_CS(obj,obj);
    
    /* some checks */
    assert(l > 0 && w > 0 && h > 0);
    /* must be lower case letters */
    if(strcmp(dir,"left") && strcmp(dir,"right"))
      Error0("Direction must be lower case!");
    assert(rstep > 0);
    
    for (d0 = 0; d0 < Nsd[0]; d0++)
    {
      min[0] = Xm + step[0]*d0;
      max[0] = Xm + step[0]*(d0+1);

      for (d1 = 0; d1 <  Nsd[1]; d1++)
      {
        min[1] = Ym + step[1]*d1;
        max[1] = Ym + step[1]*(d1+1);
        
        for (d2 = 0; d2 <  Nsd[2]; d2++)
        {
          min[2] = Zm;
          max[2] = ZM;
          
          rdown  = rmin + rstep*d2;
          rup    = rmin + rstep*(d2+1);
          
          assert(rup < diag);
          
          /* set xc to default which is DBL_MAX.
          // NOTE: it's important to set it as a double.
          // since i don't know what happens if you set it as 
          // a string par and convert it into double! */
          for (p = 0; p < 6; ++p)
          {
            Flag_T side = (Flag_T)(p);
            SCS_par_xc_length(par,"xc2");
            cs_Psetd(par,DBL_MAX);
            SCS_par_xc_length(par,"xc1");
            cs_Psetd(par,DBL_MAX);
          }
          
          
          /* set min and max parameters */
          for (p = 0; p < 6; ++p)
          {
            Flag_T side = (Flag_T)(p);
            SCS_par_min(par,0);
            cs_Psetd(par,min[0]);
            SCS_par_min(par,1);
            cs_Psetd(par,min[1]);
            SCS_par_min(par,2);
            cs_Psetd(par,min[2]);
            SCS_par_max(par,0);
            cs_Psetd(par,max[0]);
            SCS_par_max(par,1);
            cs_Psetd(par,max[1]);
            SCS_par_max(par,2);
            cs_Psetd(par,max[2]);
          }
          
          /* filling min */
          patch->min[0] = min[0];
          patch->min[1] = min[1];
          patch->min[2] = min[2];

          /* filling max */
          patch->max[0] = max[0];
          patch->max[1] = max[1];
          patch->max[2] = max[2];

          /* collocation */
          patch->collocation[0] = Chebyshev_Extrema;
          patch->collocation[1] = Chebyshev_Extrema;
          patch->collocation[2] = Chebyshev_Extrema;

          /* basis */
          patch->basis[0] = Chebyshev_Tn_BASIS;
          patch->basis[1] = Chebyshev_Tn_BASIS;
          patch->basis[2] = Chebyshev_Tn_BASIS;
          
          /* n */
          patch->n[0] = Nns[0];
          patch->n[1] = Nns[1];
          patch->n[2] = Nns[2];
          
          initialize_collocation_struct(patch,&coll_s[0],0);
          initialize_collocation_struct(patch,&coll_s[1],1);
          
          N_total = Nns[0]*Nns[1]*Nns[2];
          rU = alloc_double(N_total);
          rD = alloc_double(N_total);
          
          /* note: order matters */
          /* if top level and d2 > 0, r_up has flat surface  */
          if ( d2 != 0 && d2 == Nsd[2]-1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);              
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);
              
              X[2] = 1;
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set large xc (used in around) as opposed to down
              // (used in object) to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc2");
              cs_Psetd(par,xc);
              
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if bottom level and more than 1 split. */
          else if (d2 == 0 && Nsd[2] > 1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              X[2] = 0;
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  theta_phi_of_XY_CS(&th,&ph,X,side);
                  double r = interpolation_Ylm(reClm,imClm,lmax,th,ph);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = r;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if only one level one of them has flat surface */
          else if (d2 == 0 && Nsd[2] == 1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);

              X[2] = 0;
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  theta_phi_of_XY_CS(&th,&ph,X,side);
                  double r = interpolation_Ylm(reClm,imClm,lmax,th,ph);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = r;
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set large xc (used in around) as opposed to down
              // (used in object) to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc2");
              cs_Psetd(par,xc);
              
            }/* for (p = 0; p < 6; ++p) */
          }
          /* between the above cases they have perfect S2 surface */
          else
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              for (i = 0; i < Nns[0]; ++i)
              {
                for (j = 0; j < Nns[1]; ++j)
                {
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          Free(rU);
          Free(rD);
        }/* for (d2 = 0; d2 <  Nsd[2]; d2++) */
      }/* for (d1 = 0; d1 <  Nsd[1]; d1++) */
    }/* for (d0 = 0; d0 < Nsd[0]; d0++) */
    
    for (p = 0; p < 6; ++p)
    {
      Flag_T side = (Flag_T)(p);
            
      /* set center of patch */
      /* assuming objects are on y-axis */
      SCS_par_CS_center(par,"a");
      cs_Psetd(par,0.0);
    
      SCS_par_CS_center(par,"b");
      if (!strcmp(dir,"left"))
        cs_Psetd(par,-S/2.);
      else if (!strcmp(dir,"right"))
        cs_Psetd(par,S/2.);
      else
        Error0(NO_OPTION);
      
      SCS_par_CS_center(par,"c");
      cs_Psetd(par,0.0);
    }
  }
  
  /* populate parameters for hemisphere outermost patches
  // notation: obj_n == 0 => "left" side and obj == 1 => "right" side. */
  if (!EQL(r_outermost,0))/* if there is any outermost patch */
  for (obj_n = 0; obj_n < 2; ++obj_n)
  {
    const char *dir = (obj_n == 0 ? "left" : "right");
    char par[STR_SIZE3]  = {'\0'};
    char obj[STR_SIZE1]  = {'\0'};
    
    for (Uint p = 0; p < 6; ++p)
    {
      Flag_T side = (Flag_T)(p);
      
      /* left outermost hemisphere, doesn't have right */
      if (!strcmp(dir,"left") && side == RIGHT)
        continue;
      
      /* right outermost hemisphere, doesn't have left */
      if (!strcmp(dir,"right") && side == LEFT)
        continue;
      
      /* setting intervals */  
      /* note (X,Y,Z) in [?,?]x[?,?]x[0,1]. */
      const double Zm = 0. ,ZM = 1.;
      double Xm,XM;
      double Ym,YM;
      
      switch(side)
      {
        case UP:
          Xm = -1.;
          XM = 1.;
          Ym = so_ud_scale*(obj_n == 0 ? -1. : 0.);
          YM = so_ud_scale*(obj_n == 0 ?  0. : 1.);
        break;
        case DOWN:
          Ym = -1.;
          YM = 1.;
          Xm = so_ud_scale*(obj_n == 0 ? -1. : 0.);
          XM = so_ud_scale*(obj_n == 0 ?  0. : 1.);
        break;
        
        case LEFT:
        case RIGHT:
          Xm = -1.*so_lr_scale;
          XM = 1.*so_lr_scale;
          Ym = -1.*so_lr_scale;
          YM = 1.*so_lr_scale;
        break;
        
        case BACK:
          Xm = -1.;
          XM = 1.;
          Ym = so_bf_scale*(obj_n == 0 ? -1. : 0.);
          YM = so_bf_scale*(obj_n == 0 ?  0. : 1.);
        break;
        case FRONT:
          Xm = so_bf_scale*(obj_n == 0 ? -1. : 0.);
          XM = so_bf_scale*(obj_n == 0 ?  0. : 1.);
          Ym = -1.;
          YM = 1.;
        break;
        
        default:
        Error0(NO_OPTION);
      }
      
      /* step in each direction.
      // NOTE that Z min and max are 0 and 1 in cubed spherical. */
      double step[3] = {(XM-Xm)/Nsd[0],(YM-Ym)/Nsd[1],DBL_MAX};
      double min[3] = {0},max[3] = {0};
      double rup,rdown;
      double X[3] = {0};
      double *rU = 0, *rD = 0;
      Patch_T patch[1] = {0};
      struct Collocation_s coll_s[2] = {0};
      char parU[STR_SIZE3] = {'\0'};
      char parD[STR_SIZE3] = {'\0'};
      Uint N_total;
      
      /* find r step */
      const char *obj0 = "outermost";
      double l = S;
      double w = S;
      double h = S;
      double rmin = sqrt(Pow2(l)+Pow2(w)+Pow2(h))/2.;
      double rmax = Pgetd("grid_outermost_radius");
      double rstep = (rmax-rmin)/Nsd[2];
      
      set_object_name_split_CS(obj,obj0);
      /* some checks */
      assert(l > 0 && w > 0 && h > 0);
      if(S/2. > rmax)
        Error0("Outermost radius is too small\n");
      assert(rstep > 0);
      
      for (d0 = 0; d0 < Nsd[0]; d0++)
      {
        min[0] = Xm + step[0]*d0;
        max[0] = Xm + step[0]*(d0+1);

        for (d1 = 0; d1 <  Nsd[1]; d1++)
        {
          min[1] = Ym + step[1]*d1;
          max[1] = Ym + step[1]*(d1+1);
          
          for (d2 = 0; d2 <  Nsd[2]; d2++)
          {
            min[2] = Zm;
            max[2] = ZM;
            
            rdown  = rmin + rstep*d2;
            rup    = rmin + rstep*(d2+1);
            
            /* set xc to default which is DBL_MAX.
            // NOTE: it's important to set it as a double.
            // since i don't know what happens if you set it as 
            // a string par and convert it into double! */
            {
              SCS_par_xc_length(par,"xc2");
              cs_Psetd(par,DBL_MAX);
              SCS_par_xc_length(par,"xc1");
              cs_Psetd(par,DBL_MAX);
            }
            
            /* set min and max parameters */
            {
              SCS_par_min(par,0);
              cs_Psetd(par,min[0]);
              SCS_par_min(par,1);
              cs_Psetd(par,min[1]);
              SCS_par_min(par,2);
              cs_Psetd(par,min[2]);
              SCS_par_max(par,0);
              cs_Psetd(par,max[0]);
              SCS_par_max(par,1);
              cs_Psetd(par,max[1]);
              SCS_par_max(par,2);
              cs_Psetd(par,max[2]);
            }
              
            
            /* filling min */
            patch->min[0] = min[0];
            patch->min[1] = min[1];
            patch->min[2] = min[2];

            /* filling max */
            patch->max[0] = max[0];
            patch->max[1] = max[1];
            patch->max[2] = max[2];

            /* collocation */
            patch->collocation[0] = Chebyshev_Extrema;
            patch->collocation[1] = Chebyshev_Extrema;
            patch->collocation[2] = Chebyshev_Extrema;

            /* basis */
            patch->basis[0] = Chebyshev_Tn_BASIS;
            patch->basis[1] = Chebyshev_Tn_BASIS;
            patch->basis[2] = Chebyshev_Tn_BASIS;
            
            /* n */
            patch->n[0] = Nns[0];
            patch->n[1] = Nns[1];
            patch->n[2] = Nns[2];
            
            initialize_collocation_struct(patch,&coll_s[0],0);
            initialize_collocation_struct(patch,&coll_s[1],1);
            
            N_total = Nns[0]*Nns[1]*Nns[2];
            rU = alloc_double(N_total);
            rD = alloc_double(N_total);
            
            /* note: order matters */
            /* if top level and d2 > 0  */
            if ( d2 != 0 && d2 == Nsd[2]-1)
            {
              {
                SCS_par_sigma(parU,SigmaU)
                SCS_par_sigma(parD,SigmaD)              
                
                for (i = 0; i < Nns[0]; ++i)
                {
                  for (j = 0; j < Nns[1]; ++j)
                  {
                    for (k = 0; k < Nns[2]; ++k)
                    {
                      rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                      rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                    }
                  }
                }
                update_parameter_array(parU,rU,N_total);
                update_parameter_array(parD,rD,N_total);
              }/* for (p = 0; p < 6; ++p) */
            }
            /* if bottom level and more than 1 split, flat surface */
            else if (d2 == 0 && Nsd[2] > 1)
            {
              {
                SCS_par_sigma(parU,SigmaU)
                SCS_par_sigma(parD,SigmaD)
                
                double xc;
                if (side == UP || side == DOWN)
                  xc = h/2.;
                else if (side == LEFT || side == RIGHT)
                  xc = w;
                else if (side == BACK || side == FRONT)
                  xc = l/2.;
                else
                  Error0(NO_OPTION);
                
                X[2] = DBL_MAX;/* catch error */
                for (i = 0; i < Nns[0]; ++i)
                {
                  X[0] = point_value(i,&coll_s[0]);
                  for (j = 0; j < Nns[1]; ++j)
                  {
                    X[1] = point_value(j,&coll_s[1]);
                    for (k = 0; k < Nns[2]; ++k)
                    {
                      rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                      rD[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));;
                    }
                  }
                }
                update_parameter_array(parU,rU,N_total);
                update_parameter_array(parD,rD,N_total);
                
                /* set small xc to increase accuracy in interpolations 
                // and derivatives. */
                SCS_par_xc_length(par,"xc1");
                cs_Psetd(par,xc);
                
              }/* for (p = 0; p < 6; ++p) */
            }
            /* if only one level one of them has flat surface */
            else if (d2 == 0 && Nsd[2] == 1)
            {
              {
                SCS_par_sigma(parU,SigmaU)
                SCS_par_sigma(parD,SigmaD)
                
                double xc;
                if (side == UP || side == DOWN)
                  xc = h/2.;
                else if (side == LEFT || side == RIGHT)
                  xc = w;
                else if (side == BACK || side == FRONT)
                  xc = l/2.;
                else
                  Error0(NO_OPTION);

                X[2] = DBL_MAX;/* catch error */
                for (i = 0; i < Nns[0]; ++i)
                {
                  X[0] = point_value(i,&coll_s[0]);
                  for (j = 0; j < Nns[1]; ++j)
                  {
                    X[1] = point_value(j,&coll_s[1]);
                    for (k = 0; k < Nns[2]; ++k)
                    {
                      rU[i_j_k_to_ijk(Nns,i,j,k)] = rmax;
                      rD[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));
                    }
                  }
                }
                update_parameter_array(parU,rU,N_total);
                update_parameter_array(parD,rD,N_total);
                
                /* set small xc to increase accuracy in interpolations 
                // and derivatives. */
                SCS_par_xc_length(par,"xc1");
                cs_Psetd(par,xc);

              }
            }
            /* between the above cases they have perfect S2 surface */
            else
            {
              {
                SCS_par_sigma(parU,SigmaU)
                SCS_par_sigma(parD,SigmaD)
                
                for (i = 0; i < Nns[0]; ++i)
                {
                  for (j = 0; j < Nns[1]; ++j)
                  {
                    for (k = 0; k < Nns[2]; ++k)
                    {
                      rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                      rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                    }
                  }
                }
                update_parameter_array(parU,rU,N_total);
                update_parameter_array(parD,rD,N_total);
              }
            }
            Free(rU);
            Free(rD);
          }/* for (d2 = 0; d2 <  Nsd[2]; d2++) */
        }/* for (d1 = 0; d1 <  Nsd[1]; d1++) */
      }/* for (d0 = 0; d0 < Nsd[0]; d0++) */
    }/* for (p = 0; p < 6; ++p) */
    
    for (Uint p = 0; p < 6; ++p)
    {
      Flag_T side = (Flag_T)(p);
            
      /* set center of patch */
      /* assuming objects are on y-axis */
      SCS_par_CS_center(par,"a");
      cs_Psetd(par,0.0);
    
      SCS_par_CS_center(par,"b");
      cs_Psetd(par,0.0);
      
      SCS_par_CS_center(par,"c");
      cs_Psetd(par,0.0);
    }
  }
  
  
  #if USE_SCS_FILLING_BOX == 1
  /* populate parameters for outermost patches */
  if (!EQL(r_outermost,0))/* if there is any outermost patch */
  for (obj_n = 0; obj_n < 1; ++obj_n)
  {
    const char *const dir = "NA";
    char par[STR_SIZE3]  = {'\0'};
    char obj[STR_SIZE1]  = {'\0'};
    
    for (Uint p = 0; p < 6; ++p)
    {
    Flag_T side = (Flag_T)(p);
    
    /* note (X,Y,Z) in [-1,1]x[-1,1]x[0,1]. */
    double Xm = -1,XM = 1;
    double Ym = -1,YM = 1;
    double Zm = 0 ,ZM = 1;
    /* step in each direction, note X in [-1,1]x[-1,1]x[0,1]. 
    // also NOTE that Z min and max are 0 and 1 in cubed spherical. */
    double step[3] = {(XM-Xm)/Nsd[0],(YM-Ym)/Nsd[1],DBL_MAX};
    double min[3] = {0},max[3] = {0};
    double rup,rdown;
    double X[3] = {0};
    double *rU = 0, *rD = 0;
    Patch_T patch[1] = {0};
    struct Collocation_s coll_s[2] = {0};
    char parU[STR_SIZE3] = {'\0'};
    char parD[STR_SIZE3] = {'\0'};
    Uint N_total;
    //Uint p;
    
    /* find r step */
    const char *obj0 = "outermost";
    double l = 2*S;
    double w = 2*S;
    double h = 2*S;
    double rmin = sqrt(Pow2(l)+Pow2(w)+Pow2(h))/2.;
    double rmax = Pgetd("grid_outermost_radius");
    double rstep = (rmax-rmin)/Nsd[2];
    
    set_object_name_split_CS(obj,obj0);
    /* some checks */
    assert(l > 0 && w > 0 && h > 0);
    if(2*S > rmax)
      Error0("Outermost radius is too small\n");
    assert(rstep > 0);
    
    for (d0 = 0; d0 < Nsd[0]; d0++)
    {
      min[0] = Xm + step[0]*d0;
      max[0] = Xm + step[0]*(d0+1);

      for (d1 = 0; d1 <  Nsd[1]; d1++)
      {
        min[1] = Ym + step[1]*d1;
        max[1] = Ym + step[1]*(d1+1);
        
        for (d2 = 0; d2 <  Nsd[2]; d2++)
        {
          min[2] = Zm;
          max[2] = ZM;
          
          rdown  = rmin + rstep*d2;
          rup    = rmin + rstep*(d2+1);
          
          /* set xc to default which is DBL_MAX.
          // NOTE: it's important to set it as a double.
          // since i don't know what happens if you set it as 
          // a string par and convert it into double! */
         // for (p = 0; p < 6; ++p)
          {
            //Flag_T side = (Flag_T)(p);
            SCS_par_xc_length(par,"xc2");
            cs_Psetd(par,DBL_MAX);
            SCS_par_xc_length(par,"xc1");
            cs_Psetd(par,DBL_MAX);
          }
          
          /* set min and max parameters */
          //for (p = 0; p < 6; ++p)
          {
            //Flag_T side = (Flag_T)(p);
            SCS_par_min(par,0);
            cs_Psetd(par,min[0]);
            SCS_par_min(par,1);
            cs_Psetd(par,min[1]);
            SCS_par_min(par,2);
            cs_Psetd(par,min[2]);
            SCS_par_max(par,0);
            cs_Psetd(par,max[0]);
            SCS_par_max(par,1);
            cs_Psetd(par,max[1]);
            SCS_par_max(par,2);
            cs_Psetd(par,max[2]);
          }
            
          
          /* filling min */
          patch->min[0] = min[0];
          patch->min[1] = min[1];
          patch->min[2] = min[2];

          /* filling max */
          patch->max[0] = max[0];
          patch->max[1] = max[1];
          patch->max[2] = max[2];

          /* collocation */
          patch->collocation[0] = Chebyshev_Extrema;
          patch->collocation[1] = Chebyshev_Extrema;
          patch->collocation[2] = Chebyshev_Extrema;

          /* basis */
          patch->basis[0] = Chebyshev_Tn_BASIS;
          patch->basis[1] = Chebyshev_Tn_BASIS;
          patch->basis[2] = Chebyshev_Tn_BASIS;
          
          /* n */
          patch->n[0] = Nns[0];
          patch->n[1] = Nns[1];
          patch->n[2] = Nns[2];
          
          initialize_collocation_struct(patch,&coll_s[0],0);
          initialize_collocation_struct(patch,&coll_s[1],1);
          
          N_total = Nns[0]*Nns[1]*Nns[2];
          rU = alloc_double(N_total);
          rD = alloc_double(N_total);
          
          /* note: order matters */
          /* if top level and d2 > 0  */
          if ( d2 != 0 && d2 == Nsd[2]-1)
          {
            //for (p = 0; p < 6; ++p)
            {
              //Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU)
              SCS_par_sigma(parD,SigmaD)              
              
              for (i = 0; i < Nns[0]; ++i)
              {
                for (j = 0; j < Nns[1]; ++j)
                {
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if bottom level and more than 1 split, flat surface */
          else if (d2 == 0 && Nsd[2] > 1)
          {
            //for (p = 0; p < 6; ++p)
            {
              //Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU)
              SCS_par_sigma(parD,SigmaD)
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);
              
              X[2] = DBL_MAX;/* catch error */
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set small xc to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc1");
              cs_Psetd(par,xc);
              
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if only one level one of them has flat surface */
          else if (d2 == 0 && Nsd[2] == 1)
          {
            //for (p = 0; p < 6; ++p)
            {
              //Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU)
              SCS_par_sigma(parD,SigmaD)
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);

              X[2] = DBL_MAX;/* catch error */
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rmax;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set small xc to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc1");
              cs_Psetd(par,xc);

            }/* for (p = 0; p < 6; ++p) */
          }
          /* between the above cases they have perfect S2 surface */
          else
          {
            //for (p = 0; p < 6; ++p)
            {
              //Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU)
              SCS_par_sigma(parD,SigmaD)
              
              for (i = 0; i < Nns[0]; ++i)
              {
                for (j = 0; j < Nns[1]; ++j)
                {
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          Free(rU);
          Free(rD);
        }/* for (d2 = 0; d2 <  Nsd[2]; d2++) */
      }/* for (d1 = 0; d1 <  Nsd[1]; d1++) */
    }/* for (d0 = 0; d0 < Nsd[0]; d0++) */
    }/* for (p = 0; p < 6; ++p) */
    
    for (Uint p = 0; p < 6; ++p)
    {
      Flag_T side = (Flag_T)(p);
            
      /* set center of patch */
      /* assuming objects are on y-axis */
      SCS_par_CS_center(par,"a");
      cs_Psetd(par,0.0);
    
      SCS_par_CS_center(par,"b");
      cs_Psetd(par,0.0);
      
      SCS_par_CS_center(par,"c");
      cs_Psetd(par,0.0);
    }
  }
  
  /* populate parameters for filling boxes */
  const Flag_T fbox[] = {UP,DOWN,BACK,FRONT,UNDEFINED};
  obj_n = 0;
  while(fbox[obj_n] != UNDEFINED)
  {
    double step[3] = {0};
    char par[STR_SIZE3]  = {'\0'};
    char obj[STR_SIZE1]  = {'\0'};
    const char *obj0 = "filling_box";
    const char *dir;
    double l = 0, w = 0, h = 0;
    Flag_T side = fbox[obj_n];
    
    set_object_name_split_CS(obj,obj0);
    
    /* set length and centers */
    double cen[3] = {0};
    switch(fbox[obj_n])
    {
      case UP:
        dir = "up";
        l = S;
        w = 2*S;
        h = 0.5*S;
        step[0] = l/Nsd[0];
        step[1] = w/Nsd[1];
        step[2] = h/Nsd[2];
        /* orign center (0,0,3./4.*S) */
        cen[0] = -l/2.+step[0]/2.;
        for (d0 = 0; d0 < Nsd[0]; d0++)
        {
          cen[1]  = -w/2.+step[1]/2.;
          for (d1 = 0; d1 <  Nsd[1]; d1++)
          {
            cen[2] = (3./4.)*S -h/2.+step[2]/2.;
            for (d2 = 0; d2 <  Nsd[2]; d2++)
            {
              /* assuming objects are on y-axis */
              SCS_par_box_center(par,"a");
              cs_Psetd(par,cen[0]);
              SCS_par_box_center(par,"b");
              cs_Psetd(par,cen[1]);
              SCS_par_box_center(par,"c");
              cs_Psetd(par,cen[2]);
              
              cen[2] += step[2];
            }
            cen[1] += step[1];
          }
          cen[0] += step[0];
        }
      break;
      case DOWN:
        dir = "down";
        l = S;
        w = 2*S;
        h = 0.5*S;
        step[0] = l/Nsd[0];
        step[1] = w/Nsd[1];
        step[2] = h/Nsd[2];
        /* orign center (0,0,-3./4.*S) */
        cen[0] = -l/2.+step[0]/2.;
        for (d0 = 0; d0 < Nsd[0]; d0++)
        {
          cen[1]  = -w/2.+step[1]/2.;
          for (d1 = 0; d1 <  Nsd[1]; d1++)
          {
            cen[2] = (-3./4.)*S -h/2.+step[2]/2.;
            for (d2 = 0; d2 <  Nsd[2]; d2++)
            {
              /* assuming objects are on y-axis */
              SCS_par_box_center(par,"a");
              cs_Psetd(par,cen[0]);
              SCS_par_box_center(par,"b");
              cs_Psetd(par,cen[1]);
              SCS_par_box_center(par,"c");
              cs_Psetd(par,cen[2]);
              
              cen[2] += step[2];
            }
            cen[1] += step[1];
          }
          cen[0] += step[0];
        }
      break;
      case BACK:
        dir = "back";
        l = 0.5*S;
        w = 2*S;
        h = 2*S;
        step[0] = l/Nsd[0];
        step[1] = w/Nsd[1];
        step[2] = h/Nsd[2];
        /* orign center (-3./4.*S,0,0) */
        cen[0] = (-3./4.)*S-l/2.+step[0]/2.;
        for (d0 = 0; d0 < Nsd[0]; d0++)
        {
          cen[1]  = -w/2.+step[1]/2.;
          for (d1 = 0; d1 <  Nsd[1]; d1++)
          {
            cen[2] = -h/2.+step[2]/2.;
            for (d2 = 0; d2 <  Nsd[2]; d2++)
            {
              /* assuming objects are on y-axis */
              SCS_par_box_center(par,"a");
              cs_Psetd(par,cen[0]);
              SCS_par_box_center(par,"b");
              cs_Psetd(par,cen[1]);
              SCS_par_box_center(par,"c");
              cs_Psetd(par,cen[2]);
              
              cen[2] += step[2];
            }
            cen[1] += step[1];
          }
          cen[0] += step[0];
        }
      break;
      case FRONT:
        dir = "front";
        l = 0.5*S;
        w = 2*S;
        h = 2*S;
        step[0] = l/Nsd[0];
        step[1] = w/Nsd[1];
        step[2] = h/Nsd[2];
        /* orign center (3./4.*S,0,0) */
        cen[0] = (3./4.)*S-l/2.+step[0]/2.;
        for (d0 = 0; d0 < Nsd[0]; d0++)
        {
          cen[1]  = -w/2.+step[1]/2.;
          for (d1 = 0; d1 <  Nsd[1]; d1++)
          {
            cen[2] = -h/2.+step[2]/2.;
            for (d2 = 0; d2 <  Nsd[2]; d2++)
            {
              /* assuming objects are on y-axis */
              SCS_par_box_center(par,"a");
              cs_Psetd(par,cen[0]);
              SCS_par_box_center(par,"b");
              cs_Psetd(par,cen[1]);
              SCS_par_box_center(par,"c");
              cs_Psetd(par,cen[2]);
              
              cen[2] += step[2];
            }
            cen[1] += step[1];
          }
          cen[0] += step[0];
        }
      break;
      default:
        Error0(NO_OPTION);
    }
    step[0] = l/Nsd[0];
    step[1] = w/Nsd[1];
    step[2] = h/Nsd[2];
    
    /* some checks */
    assert(l > 0 && w > 0 && h > 0);
    
    for (d0 = 0; d0 < Nsd[0]; d0++)
    {
      for (d1 = 0; d1 <  Nsd[1]; d1++)
      {
        for (d2 = 0; d2 <  Nsd[2]; d2++)
        {
           SCS_par_box_length(par,"l");
           cs_Psetd(par,step[0]);
           
           SCS_par_box_length(par,"w");
           cs_Psetd(par,step[1]);
           
           SCS_par_box_length(par,"h");
           cs_Psetd(par,step[2]);
        }/* for (d2 = 0; d2 <  Nsd[2]; d2++) */
      }/* for (d1 = 0; d1 <  Nsd[1]; d1++) */
    }/* for (d0 = 0; d0 < Nsd[0]; d0++) */
    obj_n++;
  }
  #endif
  
  }/* binary */
  /* set parameters for all patches a single system */
  else if (Pcmps("grid_kind","SplitCubedSpherical(NS)") ||
           Pcmps("grid_kind","SplitCubedSpherical(BH)")   )
  {
  const Uint Num_Obj = 1;
  const double S  = grid_char->S;/* size of box around the single object */
  Uint obj_n;/* BH or NS */
  
  assert(strstr_i(grid_char->params[0]->dir,"center"));
  grid_char->params[0]->dir = "center";
  
  if(S < 0)
    Error0("Size must be positive.\n");
    
  /* first populate parameters only for objects NS/BH */
  for (obj_n = 0; obj_n < Num_Obj; ++obj_n)
  {
    /* note (X,Y,Z) in [-1,1]x[-1,1]x[0,1]. */
    const double Xm = -1,XM = 1;
    const double Ym = -1,YM = 1;
    const double Zm = 0 ,ZM = 1;
    /* step in each direction, note X in [-1,1]x[-1,1]x[0,1]. 
    // also NOTE that Z min and max are 0 and 1 in cubed spherical. */
    double step[3] = {(XM-Xm)/Nsd[0],(YM-Ym)/Nsd[1],DBL_MAX};
    double min[3] = {0},max[3] = {0};
    double rup,rdown;
    double th = 0,ph = 0,X[3] = {0};
    double *rU = 0, *rD = 0;
    Patch_T patch[1] = {0};
    struct Collocation_s coll_s[2] = {0};
    char parU[STR_SIZE3] = {'\0'};
    char parD[STR_SIZE3] = {'\0'};
    char par[STR_SIZE3]  = {'\0'};
    char obj[STR_SIZE1]  = {'\0'};
    Uint N_total,p;
    
    const char *obj0 = grid_char->params[obj_n]->obj;
    const char *dir  = grid_char->params[obj_n]->dir;
    double l = grid_char->params[obj_n]->l;
    double w = grid_char->params[obj_n]->w;
    double h = grid_char->params[obj_n]->h;
    const double *reClm = grid_char->params[obj_n]->relClm;
    const double *imClm = grid_char->params[obj_n]->imgClm;
    Uint lmax = grid_char->params[obj_n]->lmax;
    /* find r step */
    double diag = sqrt(Pow2(l)+Pow2(w)+Pow2(h))/2.;
    double rmin = diag;
    double rmax = grid_char->params[obj_n]->r_min;
    double rstep = (rmax-rmin)/Nsd[2];
    
    set_object_name_split_CS(obj,obj0);
    
    /* some checks */
    assert(rmax > 0);
    assert(rstep > 0);
    assert(l > 0 && w > 0 && h > 0);
    
    if(rmax < 0)
      Errors("%s object must have positive radius.\n",dir);
    if(2*rmax > S)
      Errors("%s object radius is too big.\n",dir);

    for (d0 = 0; d0 < Nsd[0]; d0++)
    {
      min[0] = Xm + step[0]*d0;
      max[0] = Xm + step[0]*(d0+1);

      for (d1 = 0; d1 <  Nsd[1]; d1++)
      {
        min[1] = Ym + step[1]*d1;
        max[1] = Ym + step[1]*(d1+1);
        
        for (d2 = 0; d2 <  Nsd[2]; d2++)
        {
          min[2] = Zm;
          max[2] = ZM;
          
          rdown  = rmin + rstep*d2;
          rup    = rmin + rstep*(d2+1);
          
          assert(rup > diag);
          
          /* set xc to default which is DBL_MAX.
          // NOTE: it's important to set it as a double.
          // since i don't know what happens if you set it as 
          // a string par and convert it into double! */
          for (p = 0; p < 6; ++p)
          {
            Flag_T side = (Flag_T)(p);
            SCS_par_xc_length(par,"xc2");
            cs_Psetd(par,DBL_MAX);
            SCS_par_xc_length(par,"xc1");
            cs_Psetd(par,DBL_MAX);
          }
          
          /* set min and max parameters */
          for (p = 0; p < 6; ++p)
          {
            Flag_T side = (Flag_T)(p);
            SCS_par_min(par,0);
            cs_Psetd(par,min[0]);
            SCS_par_min(par,1);
            cs_Psetd(par,min[1]);
            SCS_par_min(par,2);
            cs_Psetd(par,min[2]);
            SCS_par_max(par,0);
            cs_Psetd(par,max[0]);
            SCS_par_max(par,1);
            cs_Psetd(par,max[1]);
            SCS_par_max(par,2);
            cs_Psetd(par,max[2]);
          }
          
          /* filling min */
          patch->min[0] = min[0];
          patch->min[1] = min[1];
          patch->min[2] = min[2];

          /* filling max */
          patch->max[0] = max[0];
          patch->max[1] = max[1];
          patch->max[2] = max[2];

          /* collocation */
          patch->collocation[0] = Chebyshev_Extrema;
          patch->collocation[1] = Chebyshev_Extrema;
          patch->collocation[2] = Chebyshev_Extrema;

          /* basis */
          patch->basis[0] = Chebyshev_Tn_BASIS;
          patch->basis[1] = Chebyshev_Tn_BASIS;
          patch->basis[2] = Chebyshev_Tn_BASIS;
          
          /* n */
          patch->n[0] = Nns[0];
          patch->n[1] = Nns[1];
          patch->n[2] = Nns[2];
          
          initialize_collocation_struct(patch,&coll_s[0],0);
          initialize_collocation_struct(patch,&coll_s[1],1);
          
          N_total = Nns[0]*Nns[1]*Nns[2];
          rU = alloc_double(N_total);
          rD = alloc_double(N_total);
          
          /* note: order matters */
          /* if top level and d2 > 0 */
          if ( d2 != 0 && d2 == Nsd[2]-1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);              
              
              X[2] = 1;
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  theta_phi_of_XY_CS(&th,&ph,X,side);
                  double r = interpolation_Ylm(reClm,imClm,lmax,th,ph);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = r;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                    assert(r>rdown);
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if bottom level and more than 1 split, one has flat surf. */
          else if (d2 == 0 && Nsd[2] > 1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);
              
              X[2] = DBL_MAX;/* catch error */
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set small xc to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc1");
              cs_Psetd(par,xc);
              
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if only one level => one of them has flat surface */
          else if (d2 == 0 && Nsd[2] == 1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);

              X[2] = 1;
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  theta_phi_of_XY_CS(&th,&ph,X,side);
                  double r = interpolation_Ylm(reClm,imClm,lmax,th,ph);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = r;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set small xc to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc1");
              cs_Psetd(par,xc);
              
            }/* for (p = 0; p < 6; ++p) */
          }
          /* between the above cases they have perfect S2 surface */
          else
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              for (i = 0; i < Nns[0]; ++i)
              {
                for (j = 0; j < Nns[1]; ++j)
                {
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          Free(rU);
          Free(rD);
        }/* for (d2 = 0; d2 <  Nsd[2]; d2++) */
      }/* for (d1 = 0; d1 <  Nsd[1]; d1++) */
    }/* for (d0 = 0; d0 < Nsd[0]; d0++) */
    
    /* set center of patch */
    /* assuming objects are on y-axis */
    for (p = 0; p < 6; ++p)
    {
      Flag_T side = (Flag_T)(p);
      SCS_par_CS_center(par,"a");
      cs_Psetd(par,0.0);
    
      SCS_par_CS_center(par,"b");
      cs_Psetd(par,0.0);
      
      SCS_par_CS_center(par,"c");
      cs_Psetd(par,0.0);
    }
    
    /* set parameter for centeral box */
    obj0 = "central_box";
    step[0] = l/Nsd[0];
    step[1] = w/Nsd[1];
    step[2] = h/Nsd[2];
    set_object_name_split_CS(obj,obj0);
    
    if (!strcmp(dir,"center"))
    {
      Flag_T side = CENTER;
      double cen[3] = {0};
      /* orign center (0,0,0) */
      cen[0] = -l/2.+step[0]/2.;
      for (d0 = 0; d0 < Nsd[0]; d0++)
      {
        cen[1] = -w/2.+step[1]/2.;
        for (d1 = 0; d1 <  Nsd[1]; d1++)
        {
          cen[2] = -h/2.+step[2]/2.;
          for (d2 = 0; d2 <  Nsd[2]; d2++)
          {
            /* assuming objects are on y-axis */
            SCS_par_box_center(par,"a");
            cs_Psetd(par,cen[0]);
            SCS_par_box_center(par,"b");
            cs_Psetd(par,cen[1]);
            SCS_par_box_center(par,"c");
            cs_Psetd(par,cen[2]);
            
            cen[2] += step[2];
          }
          cen[1] += step[1];
        }
        cen[0] += step[0];
      }
    }
    else
      Error0(NO_OPTION);
    
    /* set lengths */
    for (d0 = 0; d0 < Nsd[0]; d0++)
    {
      for (d1 = 0; d1 <  Nsd[1]; d1++)
      {
        for (d2 = 0; d2 <  Nsd[2]; d2++)
        {
          SCS_par_box_length(par,"l");
          cs_Psetd(par,step[0]);
         
          SCS_par_box_length(par,"w");
          cs_Psetd(par,step[1]);
         
          SCS_par_box_length(par,"h");
          cs_Psetd(par,step[2]);
        }/* for (d2 = 0; d2 <  Nsd[2]; d2++) */
      }/* for (d1 = 0; d1 <  Nsd[1]; d1++) */
    }/* for (d0 = 0; d0 < Nsd[0]; d0++) */
  }
  
  /* populate parameters arounds of objects. */
  for (obj_n = 0; obj_n < Num_Obj; ++obj_n)
  {
    /* note (X,Y,Z) in [-1,1]x[-1,1]x[0,1]. */
    const double Xm = -1,XM = 1;
    const double Ym = -1,YM = 1;
    const double Zm = 0 ,ZM = 1;
    /* step in each direction, note X in [-1,1]x[-1,1]x[0,1]. 
    // also NOTE that Z min and max are 0 and 1 in cubed spherical. */
    const double step[3] = {(XM-Xm)/Nsd[0],(YM-Ym)/Nsd[1],DBL_MAX};
    double min[3] = {0},max[3] = {0};
    double rup,rdown;
    double th = 0,ph = 0,X[3] = {0};
    double *rU = 0, *rD = 0;
    Patch_T patch[1] = {0};
    struct Collocation_s coll_s[2] = {0};
    char parU[STR_SIZE3] = {'\0'};
    char parD[STR_SIZE3] = {'\0'};
    char par[STR_SIZE3]  = {'\0'};
    char obj[STR_SIZE1]  = {'\0'};
    Uint N_total,p;
    
    /* find r step */
    const char *objstem = grid_char->params[obj_n]->obj;
    const char *dir = grid_char->params[obj_n]->dir;
    double l = S;
    double w = S;
    double h = S;
    const double *reClm = grid_char->params[obj_n]->relClm;
    const double *imClm = grid_char->params[obj_n]->imgClm;
    Uint lmax = grid_char->params[obj_n]->lmax;
    double diag = sqrt(Pow2(l)+Pow2(w)+Pow2(h))/2.;
    double rmax = l/2 > MaxMag_d(w/2,h/2) ? l/2 : MaxMag_d(w/2,h/2);
    double rmin = grid_char->params[obj_n]->r_min;
    double rstep = (rmax-rmin)/Nsd[2];
    
    sprintf(obj,"%s_around",objstem);
    set_object_name_split_CS(obj,obj);
    
    /* some checks */
    assert(l > 0 && w > 0 && h > 0);
    assert(rstep > 0);
    
    for (d0 = 0; d0 < Nsd[0]; d0++)
    {
      min[0] = Xm + step[0]*d0;
      max[0] = Xm + step[0]*(d0+1);

      for (d1 = 0; d1 <  Nsd[1]; d1++)
      {
        min[1] = Ym + step[1]*d1;
        max[1] = Ym + step[1]*(d1+1);
        
        for (d2 = 0; d2 <  Nsd[2]; d2++)
        {
          min[2] = Zm;
          max[2] = ZM;
          
          rdown  = rmin + rstep*d2;
          rup    = rmin + rstep*(d2+1);
          
          assert(rup < diag);
          
          /* set xc to default which is DBL_MAX.
          // NOTE: it's important to set it as a double.
          // since i don't know what happens if you set it as 
          // a string par and convert it into double! */
          for (p = 0; p < 6; ++p)
          {
            Flag_T side = (Flag_T)(p);
            SCS_par_xc_length(par,"xc2");
            cs_Psetd(par,DBL_MAX);
            SCS_par_xc_length(par,"xc1");
            cs_Psetd(par,DBL_MAX);
          }
          
          
          /* set min and max parameters */
          for (p = 0; p < 6; ++p)
          {
            Flag_T side = (Flag_T)(p);
            SCS_par_min(par,0);
            cs_Psetd(par,min[0]);
            SCS_par_min(par,1);
            cs_Psetd(par,min[1]);
            SCS_par_min(par,2);
            cs_Psetd(par,min[2]);
            SCS_par_max(par,0);
            cs_Psetd(par,max[0]);
            SCS_par_max(par,1);
            cs_Psetd(par,max[1]);
            SCS_par_max(par,2);
            cs_Psetd(par,max[2]);
          }
          
          /* filling min */
          patch->min[0] = min[0];
          patch->min[1] = min[1];
          patch->min[2] = min[2];

          /* filling max */
          patch->max[0] = max[0];
          patch->max[1] = max[1];
          patch->max[2] = max[2];

          /* collocation */
          patch->collocation[0] = Chebyshev_Extrema;
          patch->collocation[1] = Chebyshev_Extrema;
          patch->collocation[2] = Chebyshev_Extrema;

          /* basis */
          patch->basis[0] = Chebyshev_Tn_BASIS;
          patch->basis[1] = Chebyshev_Tn_BASIS;
          patch->basis[2] = Chebyshev_Tn_BASIS;
          
          /* n */
          patch->n[0] = Nns[0];
          patch->n[1] = Nns[1];
          patch->n[2] = Nns[2];
          
          initialize_collocation_struct(patch,&coll_s[0],0);
          initialize_collocation_struct(patch,&coll_s[1],1);
          
          N_total = Nns[0]*Nns[1]*Nns[2];
          rU = alloc_double(N_total);
          rD = alloc_double(N_total);
          
          /* note: order matters */
          /* if top level and d2 > 0, r_up has flat surface  */
          if ( d2 != 0 && d2 == Nsd[2]-1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);              
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);
              
              X[2] = 1;
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set large xc to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc2");
              cs_Psetd(par,xc);
              
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if bottom level and more than 1 split. */
          else if (d2 == 0 && Nsd[2] > 1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              X[2] = 0;
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  theta_phi_of_XY_CS(&th,&ph,X,side);
                  double r = interpolation_Ylm(reClm,imClm,lmax,th,ph);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = r;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if only one level one of them has flat surface */
          else if (d2 == 0 && Nsd[2] == 1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);

              X[2] = 0;
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  theta_phi_of_XY_CS(&th,&ph,X,side);
                  double r = interpolation_Ylm(reClm,imClm,lmax,th,ph);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = r;
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set large xc to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc2");
              cs_Psetd(par,xc);
              
            }/* for (p = 0; p < 6; ++p) */
          }
          /* between the above cases they have perfect S2 surface */
          else
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU);
              SCS_par_sigma(parD,SigmaD);
              
              for (i = 0; i < Nns[0]; ++i)
              {
                for (j = 0; j < Nns[1]; ++j)
                {
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          Free(rU);
          Free(rD);
        }/* for (d2 = 0; d2 <  Nsd[2]; d2++) */
      }/* for (d1 = 0; d1 <  Nsd[1]; d1++) */
    }/* for (d0 = 0; d0 < Nsd[0]; d0++) */
    
    for (p = 0; p < 6; ++p)
    {
      Flag_T side = (Flag_T)(p);
            
      /* set center of patch */
      /* assuming objects are on y-axis */
      SCS_par_CS_center(par,"a");
      cs_Psetd(par,0.0);
    
      SCS_par_CS_center(par,"b");
      cs_Psetd(par,0.0);
      
      SCS_par_CS_center(par,"c");
      cs_Psetd(par,0.0);
    }
  }
  
  /* populate parameters for outermost patches */
  if (!EQL(r_outermost,0))/* if there is any outermost patch */
  for (obj_n = 0; obj_n < 1; ++obj_n)
  {
    /* note (X,Y,Z) in [-1,1]x[-1,1]x[0,1]. */
    const double Xm = -1,XM = 1;
    const double Ym = -1,YM = 1;
    const double Zm = 0 ,ZM = 1;
    /* step in each direction, note X in [-1,1]x[-1,1]x[0,1]. 
    // also NOTE that Z min and max are 0 and 1 in cubed spherical. */
    const double step[3] = {(XM-Xm)/Nsd[0],(YM-Ym)/Nsd[1],DBL_MAX};
    double min[3] = {0},max[3] = {0};
    double rup,rdown;
    double X[3] = {0};
    double *rU = 0, *rD = 0;
    const char *const dir = "NA";
    Patch_T patch[1] = {0};
    struct Collocation_s coll_s[2] = {0};
    char parU[STR_SIZE3] = {'\0'};
    char parD[STR_SIZE3] = {'\0'};
    char par[STR_SIZE3]  = {'\0'};
    char obj[STR_SIZE1]  = {'\0'};
    Uint N_total,p;
    
    /* find r step */
    const char *obj0 = "outermost";
    double l = S;
    double w = S;
    double h = S;
    double rmin = sqrt(Pow2(l)+Pow2(w)+Pow2(h))/2.;
    double rmax = Pgetd("grid_outermost_radius");
    double rstep = (rmax-rmin)/Nsd[2];
    
    set_object_name_split_CS(obj,obj0);
    /* some checks */
    assert(l > 0 && w > 0 && h > 0);
    if(2*S > rmax)
      Error0("Outermost radius is too small\n");
    assert(rstep > 0);
    
    for (d0 = 0; d0 < Nsd[0]; d0++)
    {
      min[0] = Xm + step[0]*d0;
      max[0] = Xm + step[0]*(d0+1);

      for (d1 = 0; d1 <  Nsd[1]; d1++)
      {
        min[1] = Ym + step[1]*d1;
        max[1] = Ym + step[1]*(d1+1);
        
        for (d2 = 0; d2 <  Nsd[2]; d2++)
        {
          min[2] = Zm;
          max[2] = ZM;
          
          rdown  = rmin + rstep*d2;
          rup    = rmin + rstep*(d2+1);
          
          /* set xc to default which is DBL_MAX.
          // NOTE: it's important to set it as a double.
          // since i don't know what happens if you set it as 
          // a string par and convert it into double! */
          for (p = 0; p < 6; ++p)
          {
            Flag_T side = (Flag_T)(p);
            SCS_par_xc_length(par,"xc2");
            cs_Psetd(par,DBL_MAX);
            SCS_par_xc_length(par,"xc1");
            cs_Psetd(par,DBL_MAX);
          }
                    
          /* set min and max parameters */
          for (p = 0; p < 6; ++p)
          {
            Flag_T side = (Flag_T)(p);
            SCS_par_min(par,0);
            cs_Psetd(par,min[0]);
            SCS_par_min(par,1);
            cs_Psetd(par,min[1]);
            SCS_par_min(par,2);
            cs_Psetd(par,min[2]);
            SCS_par_max(par,0);
            cs_Psetd(par,max[0]);
            SCS_par_max(par,1);
            cs_Psetd(par,max[1]);
            SCS_par_max(par,2);
            cs_Psetd(par,max[2]);
          }
            
          
          /* filling min */
          patch->min[0] = min[0];
          patch->min[1] = min[1];
          patch->min[2] = min[2];

          /* filling max */
          patch->max[0] = max[0];
          patch->max[1] = max[1];
          patch->max[2] = max[2];

          /* collocation */
          patch->collocation[0] = Chebyshev_Extrema;
          patch->collocation[1] = Chebyshev_Extrema;
          patch->collocation[2] = Chebyshev_Extrema;

          /* basis */
          patch->basis[0] = Chebyshev_Tn_BASIS;
          patch->basis[1] = Chebyshev_Tn_BASIS;
          patch->basis[2] = Chebyshev_Tn_BASIS;
          
          /* n */
          patch->n[0] = Nns[0];
          patch->n[1] = Nns[1];
          patch->n[2] = Nns[2];
          
          initialize_collocation_struct(patch,&coll_s[0],0);
          initialize_collocation_struct(patch,&coll_s[1],1);
          
          N_total = Nns[0]*Nns[1]*Nns[2];
          rU = alloc_double(N_total);
          rD = alloc_double(N_total);
          
          /* note: order matters */
          /* if top level and d2 > 0  */
          if ( d2 != 0 && d2 == Nsd[2]-1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU)
              SCS_par_sigma(parD,SigmaD)              
              
              for (i = 0; i < Nns[0]; ++i)
              {
                for (j = 0; j < Nns[1]; ++j)
                {
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if bottom level and more than 1 split, flat surface */
          else if (d2 == 0 && Nsd[2] > 1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU)
              SCS_par_sigma(parD,SigmaD)
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);
              
              X[2] = DBL_MAX;/* catch error */
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set small xc to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc1");
              cs_Psetd(par,xc);
              
            }/* for (p = 0; p < 6; ++p) */
          }
          /* if only one level one of them has flat surface */
          else if (d2 == 0 && Nsd[2] == 1)
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU)
              SCS_par_sigma(parD,SigmaD)
              
              double xc;
              if (side == UP || side == DOWN)
                xc = h/2.;
              else if (side == LEFT || side == RIGHT)
                xc = w/2.;
              else if (side == BACK || side == FRONT)
                xc = l/2.;
              else
                Error0(NO_OPTION);

              X[2] = DBL_MAX;/* catch error */
              for (i = 0; i < Nns[0]; ++i)
              {
                X[0] = point_value(i,&coll_s[0]);
                for (j = 0; j < Nns[1]; ++j)
                {
                  X[1] = point_value(j,&coll_s[1]);
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rmax;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = xc*sqrt(1+Pow2(X[0])+Pow2(X[1]));
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
              
              /* set small xc to increase accuracy in interpolations 
              // and derivatives. */
              SCS_par_xc_length(par,"xc1");
              cs_Psetd(par,xc);
              
            }/* for (p = 0; p < 6; ++p) */
          }
          /* between the above cases they have perfect S2 surface */
          else
          {
            for (p = 0; p < 6; ++p)
            {
              Flag_T side = (Flag_T)(p);
              SCS_par_sigma(parU,SigmaU)
              SCS_par_sigma(parD,SigmaD)
              
              for (i = 0; i < Nns[0]; ++i)
              {
                for (j = 0; j < Nns[1]; ++j)
                {
                  for (k = 0; k < Nns[2]; ++k)
                  {
                    rU[i_j_k_to_ijk(Nns,i,j,k)] = rup;
                    rD[i_j_k_to_ijk(Nns,i,j,k)] = rdown;
                  }
                }
              }
              update_parameter_array(parU,rU,N_total);
              update_parameter_array(parD,rD,N_total);
            }/* for (p = 0; p < 6; ++p) */
          }
          Free(rU);
          Free(rD);
        }/* for (d2 = 0; d2 <  Nsd[2]; d2++) */
      }/* for (d1 = 0; d1 <  Nsd[1]; d1++) */
    }/* for (d0 = 0; d0 < Nsd[0]; d0++) */
    
    for (p = 0; p < 6; ++p)
    {
      Flag_T side = (Flag_T)(p);
            
      /* set center of patch */
      /* assuming objects are on y-axis */
      SCS_par_CS_center(par,"a");
      cs_Psetd(par,0.0);
    
      SCS_par_CS_center(par,"b");
      cs_Psetd(par,0.0);
      
      SCS_par_CS_center(par,"c");
      cs_Psetd(par,0.0);
    }
  }
  }
  else
  {
    Error0(NO_OPTION);
  }
  
  FUNC_TOC  
} 

/* memory alloc patches for BBN_Split_CubedSpherical type 
// this function is deprecated. */
void alloc_patches_Split_CubedSpherical_grid(Grid_T *const grid)
{
  const Uint Np = (Uint)Pgeti(P_"SplitCS_Npatches");
  Uint i;
  
  grid->patch = calloc((Np+1),sizeof(*grid->patch));
  IsNull(grid->patch);
  
  for (i = 0; i < Np; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    IsNull(grid->patch[i]);
  }
  
}

/* memory alloc patches for single neutron star using cubed spherical + box grid */
void alloc_patches_SNS_CubedSpherical_Box_grid(Grid_T *const grid)
{
  Uint Np = 18;/* number of patches without outermost's 
                   2 sets of cubed sphere = 2*6
                   4 filling boxes
                   2 central boxes */
  Uint outermost;
  Uint i;
  
  outermost = (Uint) PgetiEZ("Number_of_Outermost_Split");
  if (outermost != (Uint)INT_MAX)
    Np += 6*outermost;
  
  grid->patch = calloc((Np+1),sizeof(*grid->patch));
  IsNull(grid->patch);
  
  for (i = 0; i < Np; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    IsNull(grid->patch[i]);
  }
  
}

/* memory alloc patches for single neutron star using cubed spherical grid */
void alloc_patches_SNS_CubedSpherical_grid(Grid_T *const grid)
{
  Uint Np = 13;/* number of patches without outermost's 
                   2 sets of cubed sphere = 2*6
                   1 central boxes */
  Uint outermost;
  Uint i;
  
  outermost = (Uint) PgetiEZ("Number_of_Outermost_Split");
  if (outermost != (Uint)INT_MAX)
    Np += 6*outermost;
  
  grid->patch = calloc((Np+1),sizeof(*grid->patch));
  IsNull(grid->patch);
  
  for (i = 0; i < Np; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    IsNull(grid->patch[i]);
  }
  
}

/* memory alloc patches for single black hole using cubed spherical grid */
void alloc_patches_SBH_CubedSpherical_grid(Grid_T *const grid)
{
  Uint Np = 6;/* number of patches without outermost's 
                   1 sets of cubed sphere = 1*6 */
  Uint outermost;
  Uint i;
  
  outermost = (Uint) PgetiEZ("Number_of_Outermost_Split");
  if (outermost != (Uint)INT_MAX)
    Np += 6*outermost;
  
  grid->patch = calloc((Np+1),sizeof(*grid->patch));
  IsNull(grid->patch);
  
  for (i = 0; i < Np; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    IsNull(grid->patch[i]);
  }
  
}

/* Jacobian transformation for split cubed spherical patch.type: OB_T_SCS
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OJ_T_SCS(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double kd[2] = {0.,1.};/* delta Kronecker */
  const double S   = patch->JacobianT->SCS->sign; /* sign */
  const Uint i = patch->JacobianT->SCS->iper;
  const Uint j = patch->JacobianT->SCS->jper;
  const Uint k = patch->JacobianT->SCS->kper;
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
  double xc1,xc2, R1,R2,
         dxc1_dx,dxc1_dy,dxc1_dz,
         dxc2_dx,dxc2_dy,dxc2_dz,
         dR1_dx,dR1_dy,dR1_dz,
         dR2_dx,dR2_dy,dR2_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(kd[i==l]-x[i]*kd[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(kd[i==l]-x[i]*kd[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(kd[i==l]-x[i]*kd[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(kd[j==l]-x[j]*kd[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(kd[j==l]-x[j]*kd[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(kd[j==l]-x[j]*kd[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      d3 = Power3(d1);
      
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      
      xc1 = S*(patch->CoordSysInfo->CubedSphericalCoord->xc1 == DBL_MAX ?
               R1/d1 :
               patch->CoordSysInfo->CubedSphericalCoord->xc1
               );
      xc2 = S*(patch->CoordSysInfo->CubedSphericalCoord->xc2 == DBL_MAX ?
               R2/d1 :
               patch->CoordSysInfo->CubedSphericalCoord->xc2
               );
      
      dR1_dx = patch->CoordSysInfo->CubedSphericalCoord->dR1_dx->v[p];
      dR2_dx = patch->CoordSysInfo->CubedSphericalCoord->dR2_dx->v[p];
      
      dxc1_dx  = (patch->CoordSysInfo->CubedSphericalCoord->xc1 == DBL_MAX ?
                  dR1_dx/d1-R1*(X[0]*JT_OJ_T_SCS(patch,_a_,_x_,p)+X[1]*JT_OJ_T_SCS(patch,_b_,_x_,p))/d3 :
                  0.);
      dxc1_dx *= S;
      
      dxc2_dx  = (patch->CoordSysInfo->CubedSphericalCoord->xc2 == DBL_MAX ?
                  dR2_dx/d1-R2*(X[0]*JT_OJ_T_SCS(patch,_a_,_x_,p)+X[1]*JT_OJ_T_SCS(patch,_b_,_x_,p))/d3 :
                  0.);
      dxc2_dx *= S;
      
      J = ((kd[k==l]-dxc1_dx)-(x[k]-xc1)*(dxc2_dx-dxc1_dx)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      d3 = Power3(d1);
      
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      
      xc1 = S*(patch->CoordSysInfo->CubedSphericalCoord->xc1 == DBL_MAX ?
               R1/d1 :
               patch->CoordSysInfo->CubedSphericalCoord->xc1
               );
      xc2 = S*(patch->CoordSysInfo->CubedSphericalCoord->xc2 == DBL_MAX ?
               R2/d1 :
               patch->CoordSysInfo->CubedSphericalCoord->xc2
               );
      
      dR1_dy = patch->CoordSysInfo->CubedSphericalCoord->dR1_dy->v[p];
      dR2_dy = patch->CoordSysInfo->CubedSphericalCoord->dR2_dy->v[p];
      
      dxc1_dy  = (patch->CoordSysInfo->CubedSphericalCoord->xc1 == DBL_MAX ?
                  dR1_dy/d1-R1*(X[0]*JT_OJ_T_SCS(patch,_a_,_y_,p)+X[1]*JT_OJ_T_SCS(patch,_b_,_y_,p))/d3 :
                  0.);
      dxc1_dy *= S;
      
      dxc2_dy  = (patch->CoordSysInfo->CubedSphericalCoord->xc2 == DBL_MAX ?
                  dR2_dy/d1-R2*(X[0]*JT_OJ_T_SCS(patch,_a_,_y_,p)+X[1]*JT_OJ_T_SCS(patch,_b_,_y_,p))/d3 :
                  0.);
      dxc2_dy *= S;
      
      J = ((kd[k==l]-dxc1_dy)-(x[k]-xc1)*(dxc2_dy-dxc1_dy)/(xc2-xc1))/(xc2-xc1);
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      d3 = Power3(d1);
      
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      
      xc1 = S*(patch->CoordSysInfo->CubedSphericalCoord->xc1 == DBL_MAX ?
               R1/d1 :
               patch->CoordSysInfo->CubedSphericalCoord->xc1);
      xc2 = S*(patch->CoordSysInfo->CubedSphericalCoord->xc2 == DBL_MAX ?
               R2/d1 :
               patch->CoordSysInfo->CubedSphericalCoord->xc2);
      
      dR1_dz = patch->CoordSysInfo->CubedSphericalCoord->dR1_dz->v[p];
      dR2_dz = patch->CoordSysInfo->CubedSphericalCoord->dR2_dz->v[p];
      
      dxc1_dz  = (patch->CoordSysInfo->CubedSphericalCoord->xc1 == DBL_MAX ?
                  dR1_dz/d1-R1*(X[0]*JT_OJ_T_SCS(patch,_a_,_z_,p)+X[1]*JT_OJ_T_SCS(patch,_b_,_z_,p))/d3 :
                  0.);
      dxc1_dz *= S;
      
      dxc2_dz  = (patch->CoordSysInfo->CubedSphericalCoord->xc2 == DBL_MAX ?
                  dR2_dz/d1-R2*(X[0]*JT_OJ_T_SCS(patch,_a_,_z_,p)+X[1]*JT_OJ_T_SCS(patch,_b_,_z_,p))/d3 :
                  0.);
      dxc2_dz *= S;
      
      J = ((kd[k==l]-dxc1_dz)-(x[k]-xc1)*(dxc2_dz-dxc1_dz)/(xc2-xc1))/(xc2-xc1);
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for split cubed spherical patch.type: OT_T_SCS
// NOTE: R2 is constant => dR2/d? = 0.
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T_SCS(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double kd[2] = {0.,1.};/* delta Kronecker */
  const double S   = patch->JacobianT->SCS->sign; /* sign */
  const Uint i = patch->JacobianT->SCS->iper;
  const Uint j = patch->JacobianT->SCS->jper;
  const Uint k = patch->JacobianT->SCS->kper;
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
  double xc1,ratio,R1,R2,
         dxc1_dx,dxc1_dy,dxc1_dz,
         dratio_dx,dratio_dy,dratio_dz,
         dR1_dx,dR1_dy,dR1_dz,
         dR2_dx,dR2_dy,dR2_dz;
  
  dA_da = get_dA_da(q2_e,q1_e);
  switch(dA_da)
  {
    case da_dx:
      l = 0;
      J = S*(kd[i==l]-x[i]*kd[k==l]/x[k])/x[k];
    break;
    case da_dy:
      l = 1;
      J = S*(kd[i==l]-x[i]*kd[k==l]/x[k])/x[k];
    break;
    case da_dz:
      l = 2;
      J = S*(kd[i==l]-x[i]*kd[k==l]/x[k])/x[k];
    break;
    case db_dx:
      l = 0;
      J = S*(kd[j==l]-x[j]*kd[k==l]/x[k])/x[k];
    break;
    case db_dy:
      l = 1;
      J = S*(kd[j==l]-x[j]*kd[k==l]/x[k])/x[k];
    break;
    case db_dz:
      l = 2;
      J = S*(kd[j==l]-x[j]*kd[k==l]/x[k])/x[k];
    break;
    case dc_dx:
      l = 0;
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      d3 = Power3(d1);
      
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      
      ratio = 1-R1/R2;
      
      //dR1_dx = patch->CoordSysInfo->CubedSphericalCoord->dR1_dx->v[p];
      //dR2_dx = patch->CoordSysInfo->CubedSphericalCoord->dR2_dx->v[p];
      
      dR1_dx = (patch->CoordSysInfo->CubedSphericalCoord->xc1 == DBL_MAX ?
                0. :
                patch->CoordSysInfo->CubedSphericalCoord->xc1 *
                (X[0]*JT_OT_T_SCS(patch,_a_,_x_,p)+X[1]*JT_OT_T_SCS(patch,_b_,_x_,p))/d1);
      dR2_dx = 0.;
      
      dratio_dx = (-dR1_dx+dR2_dx*R1/R2)/R2;
      
      xc1 = S*R1/d1;
      
      dxc1_dx  = dR1_dx/d1-R1*(X[0]*JT_OT_T_SCS(patch,_a_,_x_,p)+X[1]*JT_OT_T_SCS(patch,_b_,_x_,p))/d3;
      dxc1_dx *= S;
      
      J = ((-dxc1_dx+xc1*kd[k==l]/x[k])/x[k]-X[2]*dratio_dx)/ratio;
    break;
    case dc_dy:
      l  = 1;
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      d3 = Power3(d1);
      
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      
      ratio = 1-R1/R2;
      
      //dR1_dy = patch->CoordSysInfo->CubedSphericalCoord->dR1_dy->v[p];
      //dR2_dy = patch->CoordSysInfo->CubedSphericalCoord->dR2_dy->v[p];
      
      dR1_dy = (patch->CoordSysInfo->CubedSphericalCoord->xc1 == DBL_MAX ? 
                0. :
                patch->CoordSysInfo->CubedSphericalCoord->xc1 *
                (X[0]*JT_OT_T_SCS(patch,_a_,_y_,p)+X[1]*JT_OT_T_SCS(patch,_b_,_y_,p))/d1);
      dR2_dy = 0.;
      
      dratio_dy = (-dR1_dy+dR2_dy*R1/R2)/R2;
      
      xc1 = S*R1/d1;
      
      dxc1_dy  = dR1_dy/d1-R1*(X[0]*JT_OT_T_SCS(patch,_a_,_y_,p)+X[1]*JT_OT_T_SCS(patch,_b_,_y_,p))/d3;
      dxc1_dy *= S;
      
      J = ((-dxc1_dy+xc1*kd[k==l]/x[k])/x[k]-X[2]*dratio_dy)/ratio;
    break;
    case dc_dz:
      l  = 2;
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      d3 = Power3(d1);
      
      R1  = patch->CoordSysInfo->CubedSphericalCoord->R1_f->v[p];
      R2  = patch->CoordSysInfo->CubedSphericalCoord->R2_f->v[p];
      
      ratio = 1-R1/R2;
      
      //dR1_dz = patch->CoordSysInfo->CubedSphericalCoord->dR1_dz->v[p];
      //dR2_dz = patch->CoordSysInfo->CubedSphericalCoord->dR2_dz->v[p];
      dR1_dz = (patch->CoordSysInfo->CubedSphericalCoord->xc1 == DBL_MAX ? 
                0. :
                patch->CoordSysInfo->CubedSphericalCoord->xc1 *
                (X[0]*JT_OT_T_SCS(patch,_a_,_z_,p)+X[1]*JT_OT_T_SCS(patch,_b_,_z_,p))/d1);
      dR2_dz = 0.;
      
      dratio_dz = (-dR1_dz+dR2_dz*R1/R2)/R2;
      
      xc1 = S*R1/d1;
      
      dxc1_dz  = dR1_dz/d1-R1*(X[0]*JT_OT_T_SCS(patch,_a_,_z_,p)+X[1]*JT_OT_T_SCS(patch,_b_,_z_,p))/d3;
      dxc1_dz *= S;
      
      J = ((-dxc1_dz+xc1*kd[k==l]/x[k])/x[k]-X[2]*dratio_dz)/ratio;
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}


/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_up
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 0,j = 1,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_down
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 1,j = 0,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_left
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 0,j = 2,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_right
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 2,j = 0,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_back
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 2,j = 1,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : NS_T_CS_front
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_NS_T_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 1,j = 2,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_up
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 0,j = 1,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_down
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 1,j = 0,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_left
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 0,j = 2,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_right
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 2,j = 0,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_back
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 2,j = 1,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : SR_T_CS_front
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_SR_T_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 1,j = 2,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d3;
  Uint l;
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      d1 = sqrt(1+Pow2(X[0])+Pow2(X[1]));
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
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_up
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 0,j = 1,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  Uint l;
  double d1,L,xc1,dL_dx,dL_dy,dL_dz,R2;
  
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
      l      = 0;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dx  = (X[0]*JT_OT_T1_CS_up(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_up(patch,_b_,_x_,p))/d1;
      dL_dx *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dx)/L;
    break;
    case dc_dy:
      l      = 1;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dy  = (X[0]*JT_OT_T1_CS_up(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_up(patch,_b_,_y_,p))/d1;
      dL_dy *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dy)/L;
    break;
    case dc_dz:
      l      = 2;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dz  = (X[0]*JT_OT_T1_CS_up(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_up(patch,_b_,_z_,p))/d1;
      dL_dz *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dz)/L;
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_down
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 1,j = 0,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  Uint l;
  double d1,L,xc1,dL_dx,dL_dy,dL_dz,R2;
  
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
      l      = 0;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dx  = (X[0]*JT_OT_T1_CS_down(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_down(patch,_b_,_x_,p))/d1;
      dL_dx *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dx)/L;
    break;
    case dc_dy:
      l      = 1;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dy  = (X[0]*JT_OT_T1_CS_down(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_down(patch,_b_,_y_,p))/d1;
      dL_dy *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dy)/L;
    break;
    case dc_dz:
      l      = 2;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dz  = (X[0]*JT_OT_T1_CS_down(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_down(patch,_b_,_z_,p))/d1;
      dL_dz *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dz)/L;
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_left
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 0,j = 2,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  Uint l;
  double d1,L,xc1,dL_dx,dL_dy,dL_dz,R2;
  
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
      l      = 0;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dx  = (X[0]*JT_OT_T1_CS_left(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_left(patch,_b_,_x_,p))/d1;
      dL_dx *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dx)/L;
    break;
    case dc_dy:
      l      = 1;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dy  = (X[0]*JT_OT_T1_CS_left(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_left(patch,_b_,_y_,p))/d1;
      dL_dy *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dy)/L;
    break;
    case dc_dz:
      l      = 2;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dz  = (X[0]*JT_OT_T1_CS_left(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_left(patch,_b_,_z_,p))/d1;
      dL_dz *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dz)/L;
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_right
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 2,j = 0,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  Uint l;
  double d1,L,xc1,dL_dx,dL_dy,dL_dz,R2;
  
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
      l      = 0;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dx  = (X[0]*JT_OT_T1_CS_right(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_right(patch,_b_,_x_,p))/d1;
      dL_dx *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dx)/L;
    break;
    case dc_dy:
      l      = 1;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dy  = (X[0]*JT_OT_T1_CS_right(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_right(patch,_b_,_y_,p))/d1;
      dL_dy *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dy)/L;
    break;
    case dc_dz:
      l      = 2;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dz  = (X[0]*JT_OT_T1_CS_right(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_right(patch,_b_,_z_,p))/d1;
      dL_dz *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dz)/L;
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_back
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 2,j = 1,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  Uint l;
  double d1,L,xc1,dL_dx,dL_dy,dL_dz,R2;
  
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
      l      = 0;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dx  = (X[0]*JT_OT_T1_CS_back(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_back(patch,_b_,_x_,p))/d1;
      dL_dx *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dx)/L;
    break;
    case dc_dy:
      l      = 1;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dy  = (X[0]*JT_OT_T1_CS_back(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_back(patch,_b_,_y_,p))/d1;
      dL_dy *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dy)/L;
    break;
    case dc_dz:
      l      = 2;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dz  = (X[0]*JT_OT_T1_CS_back(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_back(patch,_b_,_z_,p))/d1;
      dL_dz *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dz)/L;
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T1_CS_front
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T1_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 1,j = 2,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  Uint l;
  double d1,L,xc1,dL_dx,dL_dy,dL_dz,R2;
  
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
      l      = 0;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dx  = (X[0]*JT_OT_T1_CS_front(patch,_a_,_x_,p)+X[1]*JT_OT_T1_CS_front(patch,_b_,_x_,p))/d1;
      dL_dx *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dx)/L;
    break;
    case dc_dy:
      l      = 1;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dy  = (X[0]*JT_OT_T1_CS_front(patch,_a_,_y_,p)+X[1]*JT_OT_T1_CS_front(patch,_b_,_y_,p))/d1;
      dL_dy *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dy)/L;
    break;
    case dc_dz:
      l      = 2;
      d1     = sqrt(1+Pow2(X[0])+Pow2(X[1]));
      R2     = patch->CoordSysInfo->CubedSphericalCoord->R2;
      xc1    = patch->CoordSysInfo->CubedSphericalCoord->xc1;
      L      = 1-S*d1*xc1/R2;
      dL_dz  = (X[0]*JT_OT_T1_CS_front(patch,_a_,_z_,p)+X[1]*JT_OT_T1_CS_front(patch,_b_,_z_,p))/d1;
      dL_dz *= -S*xc1/R2;
      J      = (xc1*K[k==l]/Pow2(x[k])-X[2]*dL_dz)/L;
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_up
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_up(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 0,j = 1,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d2;
  Uint l;
  double dd_dx,dd_dy,dd_dz,
         R1,R2,dR;
  
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
      l  = 0;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dx = X[0]*JT_OT_T2_CS_up(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_up(patch,_b_,_x_,p);
      J     = (K[k==l]/x[k]+dd_dx/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dy:
      l  = 1;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dy = X[0]*JT_OT_T2_CS_up(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_up(patch,_b_,_y_,p);
      J     = (K[k==l]/x[k]+dd_dy/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dz:
      l  = 2;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dz = X[0]*JT_OT_T2_CS_up(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_up(patch,_b_,_z_,p);
      J     = (K[k==l]/x[k]+dd_dz/d2)*S*R1/(dR*x[k]*d1);
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_down
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_down(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 1,j = 0,k = 2;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d2;
  Uint l;
  double dd_dx,dd_dy,dd_dz,
         R1,R2,dR;
    
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
      l  = 0;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dx = X[0]*JT_OT_T2_CS_down(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_down(patch,_b_,_x_,p);
      J     = (K[k==l]/x[k]+dd_dx/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dy:
      l  = 1;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dy = X[0]*JT_OT_T2_CS_down(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_down(patch,_b_,_y_,p);
      J     = (K[k==l]/x[k]+dd_dy/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dz:
      l  = 2;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dz = X[0]*JT_OT_T2_CS_down(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_down(patch,_b_,_z_,p);
      J     = (K[k==l]/x[k]+dd_dz/d2)*S*R1/(dR*x[k]*d1);
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_left
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_left(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 0,j = 2,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d2;
  Uint l;
  double dd_dx,dd_dy,dd_dz,
         R1,R2,dR;
           
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
      l  = 0;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dx = X[0]*JT_OT_T2_CS_left(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_left(patch,_b_,_x_,p);
      J     = (K[k==l]/x[k]+dd_dx/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dy:
      l  = 1;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dy = X[0]*JT_OT_T2_CS_left(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_left(patch,_b_,_y_,p);
      J     = (K[k==l]/x[k]+dd_dy/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dz:
      l  = 2;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dz = X[0]*JT_OT_T2_CS_left(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_left(patch,_b_,_z_,p);
      J     = (K[k==l]/x[k]+dd_dz/d2)*S*R1/(dR*x[k]*d1);
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_right
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_right(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 2,j = 0,k = 1;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d2;
  Uint l;
  double dd_dx,dd_dy,dd_dz,
         R1,R2,dR;
    
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
      l  = 0;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dx = X[0]*JT_OT_T2_CS_right(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_right(patch,_b_,_x_,p);
      J     = (K[k==l]/x[k]+dd_dx/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dy:
      l  = 1;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dy = X[0]*JT_OT_T2_CS_right(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_right(patch,_b_,_y_,p);
      J     = (K[k==l]/x[k]+dd_dy/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dz:
      l  = 2;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dz = X[0]*JT_OT_T2_CS_right(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_right(patch,_b_,_z_,p);
      J     = (K[k==l]/x[k]+dd_dz/d2)*S*R1/(dR*x[k]*d1);
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_back
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_back(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = -1; /* sign */
  const Uint i = 2,j = 1,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d2;
  Uint l;
  double dd_dx,dd_dy,dd_dz,
         R1,R2,dR;
    
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
      l  = 0;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dx = X[0]*JT_OT_T2_CS_back(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_back(patch,_b_,_x_,p);
      J     = (K[k==l]/x[k]+dd_dx/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dy:
      l  = 1;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dy = X[0]*JT_OT_T2_CS_back(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_back(patch,_b_,_y_,p);
      J     = (K[k==l]/x[k]+dd_dy/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dz:
      l  = 2;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dz = X[0]*JT_OT_T2_CS_back(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_back(patch,_b_,_z_,p);
      J     = (K[k==l]/x[k]+dd_dz/d2)*S*R1/(dR*x[k]*d1);
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* Jacobian transformation for cubed spherical patch.type : OT_T2_CS_front
// convention:
// _a_ = X, _b_ = Y, _c_ = Z
// ->return value: dq2/dq1 */
double JT_OT_T2_CS_front(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const Uint p)
{
  /* ds/ds = 1 */
  if (q2_e == q1_e)
    return 1;
    
  double J = 0;
  enum enum_dA_da dA_da = dA_da_UNDEFINED;
  const double *const C = patch->c;
  const double K[2] = {0.,1.};/* delta Kronecker */
  const double S = 1; /* sign */
  const Uint i = 1,j = 2,k = 0;/* permuted indices */
  const double x[3] = {patch->node[p]->x[0]-C[0],
                       patch->node[p]->x[1]-C[1],
                       patch->node[p]->x[2]-C[2]};
  const double *const X = patch->node[p]->X;
  double d1,d2;
  Uint l;
  double dd_dx,dd_dy,dd_dz,
         R1,R2,dR;
    
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
      l  = 0;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dx = X[0]*JT_OT_T2_CS_front(patch,_a_,_x_,p)+X[1]*JT_OT_T2_CS_front(patch,_b_,_x_,p);
      J     = (K[k==l]/x[k]+dd_dx/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dy:
      l  = 1;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dy = X[0]*JT_OT_T2_CS_front(patch,_a_,_y_,p)+X[1]*JT_OT_T2_CS_front(patch,_b_,_y_,p);
      J     = (K[k==l]/x[k]+dd_dy/d2)*S*R1/(dR*x[k]*d1);
    break;
    case dc_dz:
      l  = 2;
      R1 = patch->CoordSysInfo->CubedSphericalCoord->R1;
      R2 = patch->CoordSysInfo->CubedSphericalCoord->R2;
      dR = 1.-R1/R2;
      d2 = 1+Pow2(X[0])+Pow2(X[1]);
      d1 = sqrt(d2);
      dd_dz = X[0]*JT_OT_T2_CS_front(patch,_a_,_z_,p)+X[1]*JT_OT_T2_CS_front(patch,_b_,_z_,p);
      J     = (K[k==l]/x[k]+dd_dz/d2)*S*R1/(dR*x[k]*d1);
    break;
    default:
      Error0("No such an enum!\n");
  }
  
  return J;
}

/* setting up sign and indices for cubde sphere based on side*/  
void SignAndIndex_permutation_CubedSphere(const Flag_T side,Uint *const a,Uint *const b,Uint *const c,double *const s)
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
      Error0(NO_OPTION);
  }
}

/* performing a test for CS coordinates:
// test includes:
// o. 2-D interpolations of radius.
// o. x_of_X function
// o. X_of_x function. */
void test_CubedSpherical_Coordinates(Grid_T *const grid)
{
  Uint p;
  Flag_T flg = NONE;
  
  printf(Pretty0"Testing Cubed Spherical Coordinates:\n");
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch 	= grid->patch[p];
    const Flag_T type   = patch->CoordSysInfo->CubedSphericalCoord->type;
    Field_T *const R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f,
            *const R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
    double Xp[3],xp[3],R;
    Uint n;
    
    for(n = 0; n < patch->nn; ++n)
    {
      const double *const X = patch->node[n]->X;
      const double *const x = patch->node[n]->x;
      
      /* test x_of_X */
      x_of_X(xp,X,patch);
      if (!EQL(root_square(3,xp,x),0))
      {
        printf("x_of_X failed! difference = %e\n",root_square(3,xp,x));
        flg = FOUND;
      }
      
      /* test X_of_x */
      X_of_x(Xp,x,patch);
      if (!EQL(root_square(3,Xp,X),0))
      {
        printf("X_of_x failed! difference = %e\n",root_square(3,Xp,X));
        flg = FOUND;
      }
      
      /* test Radius related */
      switch (type)
      {
        case OB_T_SCS:
        case OT_T_SCS:
          R = R_interpolation_CS(R2_f,X);
          if (!EQL(root_square(1,&R,&R2_f->v[n]),0))
          {
            printf("R2 interpolation failed! difference = %e\n",root_square(1,&R,&R2_f->v[n]));
            flg = FOUND;
          }
          R = R_interpolation_CS(R1_f,X);
          if (!EQL(root_square(1,&R,&R1_f->v[n]),0))
          {
            printf("R1 interpolation failed! difference = %e\n",root_square(1,&R,&R1_f->v[n]));
            flg = FOUND;
          }
        break;
        
        case NS_T_CS:
          R = R_interpolation_CS(R2_f,X);
          if (!EQL(root_square(1,&R,&R2_f->v[n]),0))
          {
            printf("R2 interpolation failed! difference = %e\n",root_square(1,&R,&R2_f->v[n]));
            flg = FOUND;
          }
        break;
        
        case SR_T_CS:
          R = R_interpolation_CS(R1_f,X);
          if (!EQL(root_square(1,&R,&R1_f->v[n]),0))
          {
            printf("R1 interpolation failed! difference = %e\n",root_square(1,&R,&R1_f->v[n]));
            flg = FOUND;
          }
        break;
        
        default:
        break;
      }
  
    }/* end of for(n = 0; n < patch->nn; ++n) */
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
  if (flg != FOUND)
    printf("Testing Cubed Spherical Coordinates: [+].\n");
  else
    printf("Testing Cubed Spherical Coordinates: [-].\n");
}


