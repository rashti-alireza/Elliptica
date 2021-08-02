/*
// Alireza Rashti
// June 2018
*/

#include "interface.h"

static const double   Huge_Eps        = 1;
static const double   ScaleFactor     = 1E-5;/* scale factor */
static const double   PrecisionFactor = 1.0;/* precision factor */
static const Uint NFaces      = 6;/* total number of faces */
/* schematic names for surfaces correspond to FACE_T enum. */
static const char *FaceName[] = {"X == 0","X == 1","Y == 0","Y == 1",
                                 "Z == 0","Z == 1",0};
/* 
// realizing interfaces of each patch to set boundary condition
// at elliptic solver based on the type of this interface, such as
// Drichlet or Neumman.
// ->return value: EXIT_SUCCESS. */
int realize_interfaces(Grid_T *const grid)
{
  FUNC_TIC
  
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_BHBH ||
      grid->kind == Grid_SplitCubedSpherical_SNS  ||
      grid->kind == Grid_SplitCubedSpherical_SBH
     )
  {
    ri_split_cubed_spherical(grid);
  }
  else
  {
    ri_general_method(grid);
  }  
  
  FUNC_TOC
  
  return EXIT_SUCCESS;
} 

/* realizing interfaces for split cubed spherical grid kind.
// this function is required cos ri_general is very slow specially
// when the number of points and interface's are many, like split 
// cubed spherical (SCS), and also ri_general is failing for this 
// kind of grid. in this algo. we use the fact that in SCS 
// all interfaces either touch each other or reach innerB or outerB.
// furthermore, there is no overlapping patch. */
static void ri_split_cubed_spherical(Grid_T *const grid)
{
  FUNC_TIC
  
  /* keep track of counted points; 1 means counted, 0 means not. */
  Uint **point_flag = calloc(grid->np,sizeof(*point_flag));
  IsNull(point_flag);
  
  Uint p;

  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *const patch = grid->patch[p];
    
    point_flag[p] = calloc(patch->nn,sizeof(*point_flag[p]));
    IsNull(point_flag[p]);
    
    /* allocating interfaces */
    alloc_interface(patch);
    /* filling basic elements */
    fill_basics(patch);
    /* filling point[?]->N */
    fill_N(patch);
    /* set all of housK flag to zero */
    flush_houseK(patch);
  }
  
  /* find all subfaces but how they connected */  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    find_adjacent_scs(grid,grid->patch[p],point_flag[p]);
  }
  
  /* not thread safe */
  for (p = 0; p < grid->np; ++p)
  {
    set_subfaces_scs(grid,grid->patch[p]);
  }
  
  /* set df_dn flags and pair subfaces. */
  pair_subfaces_and_set_bc(grid);
  
  /* taking some precaution and adding more info at your whim */
  misc(grid);
  
  /* testing */
  test_subfaces(grid);
  
  /* printing boundary for test purposes */
  if(test_print(PRINT_INTERFACES))
    pr_interfaces(grid);
  
  /* freeing */
  free_points(grid);
  free_2d_mem(point_flag,grid->np);
  
  FUNC_TOC
}

/* realize interfaces a general method (works for many kind of grid) */
static void ri_general_method(Grid_T *const grid)
{
  FUNC_TIC
  /* keep track of counted points; 1 means counted, 0 means not. */
  Uint **point_flag = calloc(grid->np,sizeof(*point_flag));
  IsNull(point_flag);
  Uint p;

  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *const patch = grid->patch[p];
    
    point_flag[p] = calloc(patch->nn,sizeof(*point_flag[p]));
    IsNull(point_flag[p]);
    
    alloc_interface(patch);/* allocating interfaces */
    fill_basics(patch);/* filling basic elements */
    fill_N(patch);/* filling point[?]->N */
    flush_houseK(patch);/* set all of housK flag to zero */
  }  
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *const patch = grid->patch[p];
    Uint f,nf;
    
    /* populate PointSet_T  */
    nf = countf(patch->interface);
    for (f = 0; f < nf; ++f)
    {
      Interface_T *interface = patch->interface[f];
      PointSet_T **innerP, **edgeP;
      Uint i;
      
      init_Points(interface,&innerP,&edgeP);
      interface->innerP = innerP;
      interface->edgeP  = edgeP;
      
      FOR_ALL(i,innerP)
      {
        find_adjPnt(innerP[i]);
      }
      FOR_ALL(i,edgeP)
      {
        find_adjPnt(edgeP[i]);
      }
    }
  }
  
  /* find various geometry of each point */
  fill_geometry(grid,point_flag);
  
  /* freeing Point_T */
  free_points(grid);
  
  /* set df_dn flags */
  if(0) set_df_dn_and_pair(grid);/* deprecated */
  else  pair_subfaces_and_set_bc(grid);

  /* taking some precaution and adding more info at your whim */
  misc(grid);
  
  /* testing */
  test_subfaces(grid);
  
  /* printing boundary for test purposes */
  if(test_print(PRINT_INTERFACES))
    pr_interfaces(grid);
   
  /* freeing */
  free_2d_mem(point_flag,grid->np);
  FUNC_TOC
}


/* taking some precaution and adding more info at your whim */
static void misc(Grid_T *const grid)
{
  Interface_T **face;
  SubFace_T *subf;
  Uint pa,f,sf;
  
  /* for all patches */
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    /* for all faces */
    FOR_ALL(f,face)
    {
      /* for all subfaces */
      for (sf = 0; sf < face[f]->ns; ++sf)
      {
        subf = face[f]->subface[sf];
        
        if (subf->outerB)
        {
          /* these help to get segfault in case of using their value */
          subf->adjPatch = UINT_MAX;
          subf->adjFace  = UINT_MAX;
          if (grid->patch[pa]->outerB != 1)
            Error0("The patch->outerB has not been set correctly!\n");
        }
        if (subf->innerB)
        {
          /* these help to get segfault in case of using their value */
          subf->adjPatch = UINT_MAX;
          subf->adjFace  = UINT_MAX;
          if (grid->patch[pa]->innerB != 1)
            Error0("The patch->innerB has not been set correctly!\n"); 
        }
      }
    }
  }
  
}

/* testing various properties of subfaces */
static void test_subfaces(const Grid_T *const grid)
{
  Interface_T **face;
  SubFace_T *subf, *subf2;
  Flag_T flg;
  Uint pa,f,sf;
  
  /* check if all point have been found */
  FOR_ALL(pa,grid->patch)
    check_houseK(grid->patch[pa]);
  
  /* flag consistency */
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    FOR_ALL(f,face)
    {
      for (sf = 0; sf < face[f]->ns; ++sf)
      {
        subf = face[f]->subface[sf];
        
        if (subf->outerB && subf->innerB)
          Error0("Outerbound and innerboundary faces conflict.\n");
        if (subf->outerB && subf->touch)
          Error0("Outerbound and touch faces conflict.\n");
        if (subf->innerB && subf->touch)
          Error0("innerB and touch faces conflict.\n");
        if (subf->outerB && subf->copy)
          Error0("Outerbound and copy faces conflict.\n");
        if (subf->outerB && !subf->exterF)
          Error0("Outerbound and external face conflict.\n");
        if (subf->innerB && !subf->exterF)
          Error0("Innerbound and external face conflict.\n");
        if (subf->patch->pn == subf->adjPatch && 
            !subf->outerB && !subf->innerB && subf->exterF)
          Error0("patch and adjacent patch conflict.\n");
        if (subf->sameX && subf->sameY)
          Error0("Two constant for an interface.\n");
        if (subf->sameX && subf->sameZ)
          Error0("Two constant for an interface.\n");
        if (subf->sameZ && subf->sameY)
          Error0("Two constant for an interface.\n");
      }
    }
  }/* FOR_ALL(pa,grid->patch) */
  
  /* connection consistency */
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    FOR_ALL(f,face)
    {
      for (sf = 0; sf < face[f]->ns; ++sf)
      {
        subf = face[f]->subface[sf];
        
        if (subf->touch)
        {
          subf2 = get_paired_subface(subf);
          
          /* making sure two touching sufaces have different df_dn */
          if (subf->df_dn && subf2->df_dn)
            Error0("Wrong df_dn flags: both have df_dn = 1 (Neumann).\n");
          /* making sure two touching sufaces have different df_dn */
          if (!subf->df_dn && !subf2->df_dn)
            Error0("Wrong df_dn flags: both have df_dn = 0 (Dirichlet).\n");
          if (!subf2->touch)
            Error0("Wrong paired subface.\n");
        }

      }
    }
  }/* FOR_ALL(pa,grid->patch) */
  
  /* make sure not all of the subfaces have df_dn = 1
  // at least one of them must have df_dn = 0 or has outerB or innerB*/
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    flg = NO;
    FOR_ALL(f,face)
    {
      for (sf = 0; sf < face[f]->ns; ++sf)
      {
        subf = face[f]->subface[sf];
        if (subf->df_dn == 0 || subf->outerB == 1 || subf->innerB == 1)
        {
          flg = YES;
          break;
        }
      }
      if (flg == YES)
        break;
    }
    if (flg == NO)
      Error0("One of patches has all Neumann boundary condition.");
    
  }/* FOR_ALL(pa,grid->patch) */
  
}

/* setting df_dn flags inside sub-faces. this flag says if we have 
// interpolation or copy situation at an interface, 
// whether value of field is used or its derivative along normal to 
// that interface; in other words whether Drichlet or Nueman B.C applied.
// the point is no two of patches must have all Neumman B.C. otherwise you get
// inf condition number!
// moreover, it pairs the related subfaces. */
static void set_df_dn_and_pair(Grid_T *const grid)
{
  Interface_T **face;
  SubFace_T *subf;
  
  Uint pa,f;
  
  /* one can add a function here if needs specific arrangement of df_dn */

  /* go thru all patches and set one Dirichlet BC for each */
  FOR_ALL_PATCHES(pa,grid)
  {
    face = grid->patch[pa]->interface;
    
    /* set one of the subface->df_dn = 0 to rest assure 
    // there is at least one Dirichlet b.c. in each patch */
    set_one_Dirichlet_BC(face);
  
  }
  
  /* go thru all patches and setup df_dn */
  FOR_ALL_PATCHES(pa,grid)
  {
    face = grid->patch[pa]->interface;
    
    /* set the rest of df_dn's. favor Neumann BC */
    FOR_ALL(f,face)
    {
      Uint sf;
      Flag_T pristine_flg = YES;
      
      /* check if this face is pristine */
      for (sf = 0; sf < face[f]->ns; ++sf)
      {
        subf = face[f]->subface[sf];
        if (!subf->touch || strstr(subf->flags_str,"Dn:"))
        {
          pristine_flg = NO;
          break;
        }
      }
      if (pristine_flg == YES)
      {
        for (sf = 0; sf < face[f]->ns; ++sf)
        {
          subf = face[f]->subface[sf];
          /* make the whole chain and determine the last ring by null pointer and so for next pointer */
          Subf_T **chain = compose_the_chain(subf);
          Uint nc = countf(chain);
          if (nc == 1)
            Error0("There is no other subface matches to this subface.");
          set_df_dn(chain[0],1);
          free_2d(chain);
        }
      }
      else
      {
        for (sf = 0; sf < face[f]->ns; ++sf)
        {
          subf = face[f]->subface[sf];
          Subf_T **chain = 0, *ring = 0;
          Flag_T flg = NONE;
          Uint nc;
          
          if (!subf->touch || strstr(subf->flags_str,"Dn:"))
            continue;
        
          /* make the whole chain and determine the last ring by null pointer and so for next pointer */
          chain = compose_the_chain(subf);
          nc = countf(chain);
          
          if (nc == 1)
            Error0("There is no other subface matches to this subface.");
        
          /* first check if there is any clue for how to set this flags 
          // by looking to other memebers of this chain */
          for (ring = chain[1]; ring; ring = ring->next)
          {
            /* from the moment it finds a clue, it sets the afterward ones */
            if (ring->set)
            {
              if (strstr(ring->subf->flags_str,"Dn:1"))
              {
                set_df_dn(ring,1);
                flg = READY;
                break;
              }
              else if (strstr(ring->subf->flags_str,"Dn:0"))
              {
                set_df_dn(ring,0);
                flg = READY;
                break;
              }
              else
                Error0("Flags in subface are not set right!");
                
            }
          }
          
          /* if cluse didn't exist, favor Neumann bc  */
          if (flg != READY)
          {
            set_df_dn(chain[0],1);
            flg = READY;
          }
          
          if (chain != 0)
            free_2d(chain);
          
        }/* end of for (sf = 0; sf < face[f]->ns; ++sf) */
      }
    }/* end of FOR_ALL(f,face) */
  }/* end of FOR_ALL(pa,grid->patch) */
  
  /* go thru all of the patches and setup adjsn */
  FOR_ALL_PATCHES(pa,grid)
  {
    face = grid->patch[pa]->interface;
    
    /* set the rest of df_dn's. favor Neumann BC */
    FOR_ALL(f,face)
    {
      Uint sf;
      
      for (sf = 0; sf < face[f]->ns; ++sf)
      {
        subf = face[f]->subface[sf];
        
        if (!subf->touch)
          continue;
        
        SubFace_T *subf2 = find_subface(subf);/* juxtapose subface of subf */
        subf->adjsn = subf2->sn;
  
      }/* end of for (sf = 0; sf < face[f]->ns; ++sf) */
    }/* end of FOR_ALL(f,face) */
  }/* end of FOR_ALL(pa,grid->patch) */
  
}

/* set one of the subface->df_dn = 0 to rest assure 
// there is at least one Dirichlet b.c. in each patch. 
// it also sets the flags df_dn adn adjsn if there are trivial cases. */
static void set_one_Dirichlet_BC(Interface_T **const face)
{ 
  Subf_T **chain = 0;
  Uint nc;/* number of chain */
  SubFace_T *subf;
  Uint sf,f; 
  Flag_T flg = NONE,IsSet=NO;
  
  FOR_ALL(f,face)
  {
    flg   = NONE;
    IsSet = NO;
    /* check if there is any BC set for this face */
    for (sf = 0; sf < face[f]->ns; ++sf)
    {
      subf = face[f]->subface[sf];
      if (strstr(subf->flags_str,"Dn:") || !subf->touch )
      {
        flg = FOUND;
        break;
      }
    }
    if (flg == FOUND)/* if this face is used */
      continue;
    
    /* since it's pristine set Dn=0 for the whole face */
    for (sf = 0; sf < face[f]->ns; ++sf)
    {
      subf = face[f]->subface[sf];
      /* make the whole chain and determine the last ring 
      // by null pointer and so for next pointer */
      chain = compose_the_chain(subf);
      nc = countf(chain);
      if (nc == 1)
        Error0("There is no other subface matches to this subface.");
      set_df_dn(chain[0],0);
      free_2d(chain);
      IsSet = YES;
    }
    if(IsSet != YES)
    {
      char msg[999] = {'\0'};
      sprintf(msg,"Could not set Dirichlet for:\n(%s,%s)",
      face[f]->patch->name,FaceName[f]);
      Error0(msg);
    }
    break;/* => only need only one face */
  }/* end of FOR_ALL(f,face) */
    
}

/* given ring it sets df_dn flags in backward and forward direction 
// of given ring in a chain (INCLUDING the ring itself) start with df_dn.
// it takes df_dn as starting point. */
static void set_df_dn(Subf_T *const ring,const Uint df_dn)
{
  Subf_T *next = 0;
  Subf_T *prev = 0;
  const Uint AS = strlen(",Dn:?")+1;/* size of attached */
  char tail[10];
  Uint f,sl;
  
  f = df_dn;
  for (next = ring; next; next = next->next)
  {
    if (next->set == 1)
    { 
      if (next->subf->df_dn != f%2)
      {
        char msg[999] = {'\0'};
        sprintf(msg,
          "(%s, %s) <- patch\n(%s, %s) <- adjacent\nflags:\n%s\n%s\n",
        next->prev->subf->patch->name,FaceName[next->prev->subf->face],
        next->subf->patch->name,FaceName[next->subf->face],
        next->prev->subf->flags_str,next->subf->flags_str);
        Errors("Wrong df_dn flags was found!\n%s\n",msg);
      }
      else
      {
        ++f;
        continue;
      }
    }
    next->subf->df_dn = f%2;
    sl = (Uint) strlen(next->subf->flags_str)+AS;
    next->subf->flags_str = realloc(next->subf->flags_str,sl);
    IsNull(next->subf->flags_str);
    sprintf(tail,",Dn:%u",f%2);
    strcat(next->subf->flags_str,tail);
    next->set = 1;
    ++f;
  }
  
  f = df_dn+1;/* starting from ring->prev so f = df_dn+1 */
  for (prev = ring->prev; prev; prev = prev->prev)
  {
    if (prev->set == 1)
    { 
      if (prev->subf->df_dn != f%2)
        Error0("Wrong df_dn flags was found.");
      else
      {
        ++f;
        continue;
      }
    }
    
    prev->subf->df_dn = f%2;
    sl = (Uint) strlen(prev->subf->flags_str)+AS;
    prev->subf->flags_str = realloc(prev->subf->flags_str,sl);
    IsNull(prev->subf->flags_str);
    sprintf(tail,",Dn:%u",f%2);
    strcat(prev->subf->flags_str,tail);
    prev->set = 1;
    ++f;
  }
  
}

/* increase the size of chain by 2 and put the last element to null.
// ->return value: the last available ring of this chain. */
static Subf_T *add_ring(Subf_T ***chain_addrss)
{
  Subf_T *ring = 0;
  Subf_T **chain = *chain_addrss;
  Uint nc;/* number of chain */
  
  if (chain == 0)
  {
    chain = calloc(2,sizeof(*chain));
    IsNull(chain);
    chain[0] = calloc(1,sizeof(*chain[0]));
    IsNull(chain[0]);
    ring = chain[0];
    ring->n = 0;
  }
  else
  {
    nc = countf(chain);
    chain = realloc(chain,(nc+2)*sizeof(*chain));
    IsNull(chain);
    chain[nc] = calloc(1,sizeof(*chain[nc]));
    IsNull(chain[nc]);
    chain[nc+1] = 0;
    ring = chain[nc];
    ring->n = nc;
  }
  
  (*chain_addrss) = chain;
  return ring;
}

/* make the whole chain and determine the last ring of the chain by null pointer.
// furthermore, it determines the last ring by subf->next = 0 and 
// the first ring by subf->prev = 0.
// algorithm:
// it looks for subf's juxtapose and when it encounters repeated subface 
// it will stop put the next pointer to null.
// ->return value: chain composes of subfaces connected to each other. */
static Subf_T **compose_the_chain(SubFace_T *const subf1)
{
  Subf_T **chain = 0;
  Subf_T *ring1, *ring2;
  SubFace_T *subf2,*next;
  Flag_T flg;
  Uint nc = 0,r;
  
  ring1 = add_ring(&chain);
  ring1->subf = subf1;
  ring1->prev = 0;
  ++nc;
  
  if (strstr(subf1->flags_str,"Dn:"))
  {
    ring1->set = 1;
  }
    
  /* find the whole chain that determine the last ring by null pointer */
  flg = NO;
  next = subf1;
  while (flg == NO)
  {
    subf2 = find_subface(next);/* juxtapose subface of subf1 */
    
    /* check if there is any repetition */
    for (r = 0; r < nc; ++r)
      if (subf2 == chain[r]->subf)
      {
        flg = YES;
        break;
      }
      
    if (flg == YES)
      break;
      
    ring2 = add_ring(&chain);
    ring2->subf = subf2;
    ring2->prev = ring1;
    ring1->next = ring2;
    ++nc;
    
    if (strstr(subf2->flags_str,"Dn:"))
    {
      ring2->set = 1;
    }
    
    ring1 = ring2;
    next  = subf2;
  }
  
  return chain;
}      


/* getting a subface, find its correspondingly paired subface to it.
// algorithm: 
// 1. go to adjPatch and adjFace.
// 2. pick one of subfaces that are not outerboundary and are touch, sub2.
// 3. pick an adjacent point on sub and search this point in sub2.
// 4. if this point has been found return sub2, otherwise go to 2.
// note: searching of one point is enough, because you subfacess' points
// are mutually exclusive.
// note: this returning subface doesn't mean that sub and sub2 are
// necessarily paired.
// ->return value: the adjacent subface, if can't error. */
static SubFace_T *find_subface(const SubFace_T *const sub)
{
  const Interface_T *const face = 
    sub->patch->grid->patch[sub->adjPatch]->interface[sub->adjFace];
  SubFace_T *sub2 = 0;
  Uint po;
  Uint s;
  Flag_T flg = NONE;
  
  if (sub->copy)
  {
    po = sub->adjid[0];
    for (s = 0; s <face->ns; ++s)
    {
      sub2 = face->subface[s];
      if (sub2->copy == 1)
      {
        Uint i;
        for (i = 0; i < sub2->np; ++i)
        {
          if (sub2->id[i] == po)
          {
            flg = FOUND;
            return sub2;
          }
        }
      }
    }/* end of if (sub2->outerB == 0 && sub2->copy == 1) */
  }/* end of if (sub->copy == 1) */
  
  /* if the other interface has not any subface to be interpolated 
  // to subface, choose one of the subfaces which needs interpolating,
  // preferably the closest one.
  // I think that's fine, since sub after all will be interpolated 
  // to this interface 
  */
  else
  {
    Node_T **const node = sub->patch->grid->patch[sub->adjPatch]->node;
    double dis = DBL_MAX;
    SubFace_T *sub3;
    
    for (s = 0; s <face->ns; ++s)
    {
      sub3 = face->subface[s];
      if (sub3->touch)
      {
        Uint i;
        
        for (i = 0; i < sub3->np; ++i)
        {
          double *x = node[sub3->id[i]]->x;
          double r = root_square(3,x,0);
          
          if (LSS(r,dis))
          {
            dis = r;
            sub2 = sub3;
            flg = FOUND;
          }
        }
        
        break;
      }/* end of if (sub3->touch == 1) */
    }/* for (s = 0; s <face->ns; ++s) */
  }/* end of else */
  
  if (flg == NONE)
  {
    char msg[999] = {'\0'};
    sprintf(msg,"(%s,%s) <- patch\n(%s,%s) <- adjacent\n%s <- flags.",
    sub->patch->name,FaceName[sub->face],
    face->patch->name,FaceName[face->fn],sub->flags_str);
    Errors("The related subface could not be found!\n%s\n",msg);
  }
    
  return sub2;
}

/* having lead and pnt, see which subface best fit for this point */
static void add_to_subface(const Point_T *const pnt,const char *const lead)
{
  Interface_T *const face = pnt->patch->interface[pnt->face];
  SubFace_T *subface;
  Flag_T flg;
  Uint sf;
  
  /* compare lead with the previous sub-faces */
  flg = NONE;
  for (sf = 0; sf < face->ns; ++sf)
  {
    subface = face->subface[sf];
    if (strcmp_i(subface->flags_str,lead))
    {
      flg = FOUND;
      break;
    }
  }
  if (flg == NONE)
  {
    face->subface = 
      realloc(face->subface,(face->ns+1)*sizeof(*face->subface));
    IsNull(face->subface);
    face->subface[face->ns] = 
      calloc(1,sizeof(*face->subface[face->ns]));
    IsNull(face->subface[face->ns]);
    subface = face->subface[face->ns];
    face->ns++;
    
    /* setting flags of this new sub face */
    subface->patch     = pnt->patch;
    subface->flags_str = dup_s(lead);
    subface->sn        = face->ns-1;
    subface->face      = pnt->face;
    subface->adjFace   = pnt->adjFace;
    subface->adjPatch  = pnt->adjPatch;
    subface->sameX     = pnt->sameX;
    subface->sameY     = pnt->sameY;
    subface->sameZ     = pnt->sameZ;
    subface->touch     = pnt->touch;
    subface->copy      = pnt->copy;
    subface->exterF    = pnt->exterF;
    subface->outerB    = pnt->outerB;
    subface->innerB    = pnt->innerB;
    subface->precision_factor = PrecisionFactor;
    
  }
  
  /* add this point to the subface */
  add_point(subface,pnt);
}

/* add this point to the best fit of subfaces
// note: this also set pnt->houseK */
static void add_to_subface_scs(Point_T *const pnt)
{
  Interface_T *const face = pnt->patch->interface[pnt->face];
  SubFace_T *subface;
  char flags[999] = {'\0'};
  char tmp[999]   = {'\0'};
  Flag_T found;
  Uint sf;
  
  /* fill flags */
  sprintf(flags,"p%uf%u",pnt->patch->pn,pnt->face);
  /* touch */
  if (pnt->touch == 1) 	strcat(flags,"t1");
  else 			strcat(flags,"t0");
  /* copy */
  if (pnt->copy == 1)   strcat(flags,"c1");
  else 			strcat(flags,"c0");
  if (strstr(flags,"t1c0"))/* interpolation onto adjFace */
  {
    tmp[0] = '\0';
    sprintf(tmp,"ap%uaf%u",pnt->adjPatch,pnt->adjFace);
    strcat(flags,tmp);
    /* sameX */
    if (pnt->sameX == 1)  strcat(flags,"X1");
    else 		  strcat(flags,"X0");
    /* sameY */
    if (pnt->sameY == 1)  strcat(flags,"Y1");
    else		  strcat(flags,"Y0");
    /* sameZ */
    if (pnt->sameZ == 1)  strcat(flags,"Z1");
    else		  strcat(flags,"Z0");
  }
  else if (strstr(flags,"t1c1"))/* copy to adj */
  {
    tmp[0] = '\0';
    sprintf(tmp,"ap%uaf%u",pnt->adjPatch,pnt->adjFace);
    strcat(flags,tmp);
    strcat(flags,"X-Y-Z-");
  }
  else/* interpolation inside adjPatch.note: touch:0,copy:0 doen't exist*/
  {
    tmp[0] = '\0';
    sprintf(tmp,"ap%uaf-",pnt->adjPatch);
    strcat(flags,tmp);
    strcat(flags,"X-Y-Z-");
  }
  /* outerB */
  if (pnt->outerB == 1) strcat(flags,"oB1");
  else 			strcat(flags,"oB0");
  /* innerB */
  if (pnt->innerB == 1) strcat(flags,"iB1");
  else 			strcat(flags,"iB0");
  /* exterF */
  if (pnt->exterF == 1) strcat(flags,"eF1");
  else 			strcat(flags,"eF0");
  
  /* compare flags with the previous sub-faces */
  found = NONE;
  for (sf = 0; sf < face->ns; ++sf)
  {
    subface = face->subface[sf];
    if (!strcmp(subface->flags_str,flags))
    {
      found = FOUND;
      break;
    }
  }
  if (found == NONE)
  {
    face->subface = 
      realloc(face->subface,(face->ns+1)*sizeof(*face->subface));
    IsNull(face->subface);
    face->subface[face->ns] = 
      calloc(1,sizeof(*face->subface[face->ns]));
    IsNull(face->subface[face->ns]);
    subface = face->subface[face->ns];
    face->ns++;
    
    /* setting flags of this new sub face */
    subface->patch     = pnt->patch;
    subface->flags_str = dup_s(flags);
    subface->sn        = face->ns-1;
    subface->adjsn     = UINT_MAX;
    subface->face      = pnt->face;
    subface->adjFace   = pnt->adjFace;
    subface->adjPatch  = pnt->adjPatch;
    subface->sameX     = pnt->sameX;
    subface->sameY     = pnt->sameY;
    subface->sameZ     = pnt->sameZ;
    subface->touch     = pnt->touch;
    subface->copy      = pnt->copy;
    subface->exterF    = pnt->exterF;
    subface->outerB    = pnt->outerB;
    subface->innerB    = pnt->innerB;
    
    subface->precision_factor = PrecisionFactor;
  }
  
  /* add this point to the subface */
  add_point_scs(subface,pnt);
  pnt->houseK = 1;
}

/* adding one point adequately to the pertinent sub-face */
static void add_point(SubFace_T *const subface,const Point_T *const pnt)
{
  subface->id = 
    realloc(subface->id,(subface->np+1)*sizeof(*subface->id));
  IsNull(subface->id);
  subface->id[subface->np] = pnt->ind;
  
  if (pnt->copy == 1)
  {
    subface->adjid = 
      realloc(subface->adjid,(subface->np+1)*sizeof(*subface->adjid));
    IsNull(subface->adjid);
    subface->adjid[subface->np] = pnt->adjPoint->ind;
  }
  
  subface->np++;
}

/* adding one point adequately to the pertinent sub-face */
static void add_point_scs(SubFace_T *const subface,const Point_T *const pnt)
{
  subface->id = 
    realloc(subface->id,(subface->np+1)*sizeof(*subface->id));
  IsNull(subface->id);
  subface->id[subface->np] = pnt->ind;
  
  if (pnt->copy == 1)
  {
    subface->adjid = 
      realloc(subface->adjid,(subface->np+1)*sizeof(*subface->adjid));
    IsNull(subface->adjid);
    subface->adjid[subface->np] = pnt->adjInd;
  }
  
  subface->np++;
}

/* inspecting flags and put the result into a string 
// for purposes of convenient searching. note: the flags which are not
// defined, get '-' char.
// ->result value: string made of flags, null of it is unsuccessful.
*/
static char *inspect_flags(const Point_T *const pnt)
{
  char str[1000] = {'\0'}, *ret = 0;
  char tmp[100] = {'\0'};
  
  sprintf(str,"patch:%u,face:%u,",pnt->patch->pn,pnt->face);
  
  /* touch */
  if (pnt->touch == 1) 	strcat(str,"touch:1,");
  else 			strcat(str,"touch:0,");
  
  /* copy */
  if (pnt->copy == 1)   strcat(str,"copy:1,");
  else 			strcat(str,"copy:0,");
  
  if (strstr(str,"touch:1,copy:0"))/* interpolation onto adjFace */
  {
    tmp[0] = '\0';
    sprintf(tmp,"adjPatch:%u,adjFace:%u,",pnt->adjPatch,pnt->adjFace);
    strcat(str,tmp);
    /* sameX */
    if (pnt->sameX == 1)  strcat(str,"sameX:1,");
    else 		  strcat(str,"sameX:0,");
    
    /* sameY */
    if (pnt->sameY == 1)  strcat(str,"sameY:1,");
    else		  strcat(str,"sameY:0,");
    /* sameZ */
    if (pnt->sameZ == 1)  strcat(str,"sameZ:1,");
    else		  strcat(str,"sameZ:0,");
    
  }
  else if (strstr(str,"touch:1,copy:1,"))/* copy to adj */
  {
    tmp[0] = '\0';
    sprintf(tmp,"adjPatch:%u,adjFace:%u,",pnt->adjPatch,pnt->adjFace);
    strcat(str,tmp);
    strcat(str,"sameX:-,sameY:-,sameZ:-,");
  }
  else/* interpolation inside adjPatch.note: touch:0,copy:0 doen't exist*/
  {
    tmp[0] = '\0';
    sprintf(tmp,"adjPatch:%u,adjFace:-,",pnt->adjPatch);
    strcat(str,tmp);
    strcat(str,"sameX:-,sameY:-,sameZ:-,");
  }
  
  /* outerB */
  if (pnt->outerB == 1) strcat(str,"outerB:1,");
  else 			strcat(str,"outerB:0,");
  
  /* innerB */
  if (pnt->innerB == 1) strcat(str,"innerB:1,");
  else 			strcat(str,"innerB:0,");
  
  /* exterF */
  if (pnt->exterF == 1) strcat(str,"exterF:1,");
  else 			strcat(str,"exterF:0,");
  
  ret = dup_s(str);  
  return ret;
}

/* filling the geometry of point struct */
static void fill_geometry(Grid_T *const grid,Uint **const point_flag)
{
  sFunc_PtoV_T **func;
  Uint i;
  
  init_func_PtoV(&func);/* initialize struct */
  /* adding func to struct to be called. each coord must have 
  // its own func. note, also external face and inner face must be found
  // in first place; thus, for each new coord sys one must add below 
  // the related functions for these two purposes.
  */
  add_func_PtoV(&func,FindInnerB_Cartesian_coord,"FindInnerB",Cartesian);
  add_func_PtoV(&func,FindExterF_Cartesian_coord,"FindExterF",Cartesian);
  add_func_PtoV(&func,FindInnerB_CS_coord,"FindInnerB",CubedSpherical);
  add_func_PtoV(&func,FindExterF_CS_coord,"FindExterF",CubedSpherical);
  
  /* calling function to find external and internal faces*/
  FOR_ALL(i,grid->patch)
  {
    Patch_T *const patch = grid->patch[i];
    run_func_PtoV(func,"FindExterF",patch);/* find external faces */
    run_func_PtoV(func,"FindInnerB",patch);/* find inner boundary */
  }
  free_func_PtoV(func);/* freeing func struct */

  FOR_ALL(i,grid->patch)
  {
    realize_neighbor(grid->patch[i],point_flag);
  }
}

/* study neighbor of points to find quantities like 
// point->face, point->touch, point->adjPoint and others.
// algorithm:
// for each interface do
// 	a. seprate each interface to its inner and edge points
// 	b. realize the adjacency of inner points
// 	c. realize the adjacency of edge points
// the reason the inner and edge points are separated is they
// use different algorithm, since edge points need be treated specially
// ->return value-> EXIT_SUCCESS
*/
static int realize_neighbor(Patch_T *const patch,Uint **const point_flag)
{
  const int nf = (int)countf(patch->interface);
  int f;
  
  /* to maximize the outer/inner boundary subfaces 
  // ,good for elliptic solver, start with max f */
  for (f = nf-1; f >= 0 ; f--)
  {
    Interface_T *interface = patch->interface[f];
    PointSet_T **innerP    = interface->innerP;
    PointSet_T **edgeP     = interface->edgeP;
    
    realize_adj(innerP,point_flag);/* first, realizing the adjacency, INNER one */
    realize_adj(edgeP,point_flag);/* second, realizing the adjacency, EDGE one */
    
    free_PointSet(innerP);
    free_PointSet(edgeP);
    interface->innerP = 0;
    interface->edgeP  = 0;
  }
  
  return EXIT_SUCCESS;
}

/* realizing adjacency of given points
// ->return value: EXIT_SUCCESS.
*/
static int realize_adj(PointSet_T **const Pnt,Uint **const point_flag)
{
  Uint p;
  char *lead;
  
  FOR_ALL(p,Pnt)
  {
    Point_T *const pnt = Pnt[p]->Pnt;
    
    /* first take care of non external faces */
    if (Pnt[p]->Pnt->exterF == 0)
    {
      lead = inspect_flags(Pnt[p]->Pnt);
      add_to_subface(Pnt[p]->Pnt,lead);
      free(lead);
      Pnt[p]->Pnt->houseK = 1;
      point_flag[pnt->patch->pn][pnt->ind] = 1;
    }
    /* second take care of inner boundary */
    if (Pnt[p]->Pnt->innerB == 1)
    {
      lead = inspect_flags(Pnt[p]->Pnt);
      add_to_subface(Pnt[p]->Pnt,lead);
      free(lead);
      Pnt[p]->Pnt->houseK = 1;
      point_flag[pnt->patch->pn][pnt->ind] = 1;
    }
    /* if this point on this face is already considered. */
    if (pnt->houseK)
      continue;
    /* although this point on face has not been considered, it might
    // be already considered for other faces, so let's skip this
    // inorder to have only one BC for each edge points at the solver. */
    if (point_flag[pnt->patch->pn][pnt->ind]) 
      continue;
    
    analyze_adjPnt(Pnt[p],point_flag);
  }
  
  return EXIT_SUCCESS;
}

/* analyzing the adjacent point of point pnt; 
// see which one best describe the adjacent boundary of the interface
// and then fill the flags of Point_T. */
static void analyze_adjPnt(PointSet_T *const Pnt,Uint **const point_flag)
{
  Point_T *const p1 = Pnt->Pnt;
  Point_T *p2;
  AdjPoint_T *adjp1;
  char *lead = 0;
  
  /* NOTE: these following if's are order dependent, 
  // BE CAREFUL if you wanna change something. */
  
  /* this point is on outer boundary */
  if (IsOutBndry(Pnt))
  {
    p1->outerB = 1;
    p1->patch->outerB = 1;
    
    lead = inspect_flags(p1);
    add_to_subface(p1,lead);
  }
  /* this point is on inner boundary */
  else if (IsInnBndry(Pnt))
  {
    p1->innerB = 1;
    
    lead = inspect_flags(p1);
    add_to_subface(p1,lead);
  }
  /* if none of the points are on the interface */
  else if (IsOverlap(Pnt))
  {
    p1->touch = 0;
    p1->adjPatch = Pnt->adjPnt[0].p;
    p1->copy  = 0;
    
    lead = inspect_flags(p1);
    add_to_subface(p1,lead);
  }
  /* if this point can be matched with the susbface made by inner points.
  // this helps setting of B.C. between patches and low the condition number. */
  else if (IsMatchedOtherInnerSubface(Pnt))
  {
    return;
  }
  /* which point best meets the normal conditon (N2.N1 = -1)  */
  else if (IsNormalFit(Pnt))
  {
    adjp1 = &Pnt->adjPnt[Pnt->idFit];
    
    p1->touch = 1;
    p1->copy  = 1;
    p1->adjPatch = adjp1->p;
    p1->adjFace  = adjp1->CopyFace;
    p2 = get_p2(Pnt);
    assert(p2);
    p1->adjPoint = p2;
    
    if (p2->houseK == 0)
    {
      p2->touch = 1;
      p2->copy  = 1;
      p2->adjPatch = p1->patch->pn;
      p2->adjFace  = p1->face;
      p2->adjPoint = p1;
      p2->houseK = 1;
      point_flag[p2->patch->pn][p2->ind] = 1;
      lead = inspect_flags(p2);
      add_to_subface(p2,lead);
      free(lead);
    }
    
    lead = inspect_flags(p1);
    add_to_subface(p1,lead);
    
  }
  /* cases in which N1.N2 = 0 and reach outerbound */
  else if (IsOrthOutBndry(Pnt))
  {
    p1->outerB = 1;
    p1->patch->outerB = 1;
    lead = inspect_flags(p1);
    add_to_subface(p1,lead);
    
  }
  /* cases in which N1.N2 = 0 and reach innerbound */
  else if (IsOrthInnBndry(Pnt))
  {
    p1->innerB = 1;
    lead = inspect_flags(p1);
    add_to_subface(p1,lead);
  }
  /* cases in which the points needs interpolation */
  else if (IsInterpolation(Pnt))
  {
    p1->copy = 0;
    p1->adjPatch = Pnt->adjPnt[Pnt->idInterp].p;

    if (Pnt->overlap == 0) p1->touch = 1;
    else  		   p1->touch = 0;
    
    if (Pnt->adjPnt[Pnt->idInterp].FaceFlg == 1)
    {
      assert(Pnt->overlap == 0);
      p1->adjFace = Pnt->adjPnt[Pnt->idInterp].InterpFace;
      set_sameXYZ(p1,p1->adjFace);
    }
    
    lead = inspect_flags(p1);
    add_to_subface(p1,lead);
    
  }
  else
  {
    Uint ind = p1->ind;
    double *x = p1->patch->node[ind]->x;
    
    fprintf(stderr,"\nThis point(%f,%f,%f) at patch = '%s' "
      "and face = '%d' "
    " \nhas not been found.\n",
                    x[0],x[1],x[2],p1->patch->name,p1->face);
    Error0("Incomplete function.\n");
  }
  p1->houseK = 1;
  point_flag[p1->patch->pn][p1->ind] = 1;
  
  if (lead) free(lead);
}

/* if there is other subfaces in this face and this point Pnt
// can be matched with the same adjPatch then use the same
// adjPatch. it helps the B.C between the patches and lowers 
// the condition number. in a nut shell, this function helps all of 
// the possible subfaces on a same face use the same adjPatch.
// algorithm:
// 1. it first finds all of the subfaces on the same Pnt->Pnt->face 
// with the same adjPatch as the adjPnt has.
// 2. it then studies normal vector on this adjPatches to see if
// the match can be made.
// 3. if normal vector mathched, take this adjPatch as the matched one
// and return 1, else return 0.
// ->return value: 1 if yes; 0 otherwise.
*/
static int IsMatchedOtherInnerSubface(PointSet_T *const Pnt)
{
  /* this doesn't work for cubed spherical, investigate later. */
  return 0;
  
  if (!Pnt->NadjPnt) return 0;
  
  if (strcmp_i(PgetsEZ("Interface_BC_Maximum_Face_Match"),"no"))
    return 0;
  
  if (Pnt->type == INNER)
    return 0;
  
  Point_T *const p1 = Pnt->Pnt;
  Point_T *p2;
  AdjPoint_T *adjp1;
  Interface_T *const face = p1->patch->interface[p1->face];
  SubFace_T *subface;
  PointSet_T Pnt_copy = Pnt[0];
  Flag_T flg;
  char *lead = 0;
  Uint i,sf;
  
  /* go thru all adjPnt */
  for (i = 0; i < Pnt->NadjPnt; i++)
  {
    Uint adjPatch = Pnt->adjPnt[i].p;
    
    /* study a similar point but with only one adjPnt */
    Pnt_copy.NadjPnt = 1;
    Pnt_copy.adjPnt  = &Pnt->adjPnt[i];
    
    flg = NONE;
    /* study each subface on the face Pnt->Pnt->face with the same adjPatch */
    for (sf = 0; sf < face->ns; ++sf)
    {
      subface = face->subface[sf];
      
      /* we privilege the adjPnt with same adjPatch as the subface */
      if (subface->adjPatch != adjPatch)
        continue;
        
      /* which point best meets the normal conditon (N2.N1 = -1)  */
      if (IsNormalFit(&Pnt_copy))
      {
        Pnt->idFit    = i;
        Pnt->idOrth   = UINT_MAX;
        Pnt->idInterp = UINT_MAX;
        adjp1 = Pnt_copy.adjPnt;
        
        p1->touch = 1;
        p1->copy  = 1;
        p1->adjPatch = adjp1->p;
        p1->adjFace  = adjp1->CopyFace;
        p2 = get_p2(Pnt);
        assert(p2);
        p1->adjPoint = p2;
        
        if (p2->houseK == 0)
        {
          p2->touch = 1;
          p2->copy  = 1;
          p2->adjPatch = p1->patch->pn;
          p2->adjFace  = p1->face;
          p2->adjPoint = p1;
          p2->houseK = 1;
          
          lead = inspect_flags(p2);
          add_to_subface(p2,lead);
          free(lead);
        }
        
        lead = inspect_flags(p1);
        add_to_subface(p1,lead);
        
        flg = FOUND;
        break;
      }
      /* cases in which the points needs interpolation */
      else if (IsInterpolation(&Pnt_copy))
      {
        Pnt->idFit    = UINT_MAX;
        Pnt->idOrth   = UINT_MAX;
        Pnt->idInterp = i;
        Pnt->overlap = Pnt_copy.overlap;
        
        adjp1 = Pnt_copy.adjPnt;
        
        p1->copy = 0;
        p1->adjPatch = adjp1->p;

        if (Pnt_copy.overlap == 0)  p1->touch = 1;
        else  		      	    p1->touch = 0;
        
        if (adjp1->FaceFlg == 1)
        {
          assert(Pnt_copy.overlap == 0);
          p1->adjFace = adjp1->InterpFace;
          set_sameXYZ(p1,p1->adjFace);
        }
        
        lead = inspect_flags(p1);
        add_to_subface(p1,lead);
        
        flg = FOUND;
        break;
      }
      
    }/* end of for (sf = 0; sf < face->ns; ++sf) */
    if (flg == FOUND)
      break;
  }/* end of for (i = 0; i < Pnt->NadjPnt; i++) */
 
  if (lead) free(lead); 
  
  if (flg == FOUND)
  {
    return 1;
  }
  
  return 0;
}

/* does this need interpolation?
// algorithm:
// 1. if: adjPnt is not on an interface safely pick this patch as interpolation
// 2. else: pick all of adjPnts which are on an interface, chose the closest point
// on the same interface and get its normal, if dot product of this normal
// with Pnt normal is negative, save it and finally pick the one which
// has this dot product closest to -1.
// ->return value: 1 if yes; 0 otherwise.
*/
static int IsInterpolation(PointSet_T *const Pnt)
{
  struct Interp_S
  {
    Uint adjid;/* adjPnt id */
    Uint fid;/* adjPnt face id */
    double dot;/* N1 dot N2 */
  }*interp = 0, *interpEq = 0;
  double tmp;
  Uint i,j,id,f;
  
  /* check if Pnt reach inside any adj patch */
  for (i = 0; i < Pnt->NadjPnt; i++)
  {
    if (Pnt->adjPnt[i].FaceFlg == 0)
    {
      Pnt->idInterp = i;
      Pnt->overlap = 1;
      return 1;
    }
  }
  
  /* find the all related N1dotN2 */
  j = 0;
  for (i = 0; i < Pnt->NadjPnt; i++)
  {
    if (Pnt->adjPnt[i].FaceFlg == 1)
    {
      for (f = 0; f < TOT_FACE; f++)
      {
        struct Face_S *const fs = &Pnt->adjPnt[i].fs[f];
        if (fs->on_f == 1 && LSS(fs->N1dotN2,0))
        {
          interp = realloc(interp,(j+1)*sizeof(*interp));
          IsNull(interp);
          interp[j].adjid = i;
          interp[j].fid   = f;
          interp[j].dot   = fs->N1dotN2;
          j++;
        }
      }/* for (j = 0; j < TOT_FACE; j++) */
    }
  }/* end of for (i = 0; i < Pnt->NadjPnt; i++) */
  
  /* find the closest N1dotN2 to -1 */
  if (j > 0)
  {
    Uint e, id2 = UINT_MAX;
    struct Interp_S *int_S;
    Flag_T flg = NONE;
    
    id = 0;
    tmp = interp[id].dot;
    
    for (i = 1; i < j; i++)
      if (LSS(interp[i].dot,tmp))
        id = i;
        
    e = 0;
    interpEq = realloc(interpEq,(e+1)*sizeof(*interpEq));
    IsNull(interpEq);
    interpEq[e].adjid = interp[id].adjid;
    interpEq[e].fid   = interp[id].fid;
    interpEq[e].dot   = interp[id].dot;
    e++;
   
    for (i = 0; i < j; i++)
    {
      if (i == id) continue;
      
      if (EQL(interp[i].dot,interp[id].dot))
      {
        interpEq = realloc(interpEq,(e+1)*sizeof(*interpEq));
        IsNull(interpEq);
        interpEq[e].adjid = interp[i].adjid;
        interpEq[e].fid   = interp[i].fid;
        interpEq[e].dot   = interp[i].dot;
        e++;
      }
    }

    if (e > 1)
    {
      Point_T p1 = *Pnt->Pnt;
      char *lead = 0;
      
      p1.touch = 1;
      p1.copy  = 0;
      
      for (i = 0; i < e; i++)
      {
        int_S = &interpEq[i];
        p1.adjPatch = Pnt->adjPnt[int_S->adjid].p;
        p1.adjFace = int_S->fid;
        set_sameXYZ(&p1,p1.adjFace);
        
        lead = inspect_flags(&p1);
        if (IsOnSubface(&p1,lead))
        {
          free(lead);
          id2 = i;
          flg = FOUND;
          break;
        }
        free(lead);

      }/* end of for (i = 0; i < e; i++) */
    }/* end of if (e > 1)*/
    
    if (flg == FOUND)  int_S = &interpEq[id2];
    else	       int_S = &interp[id];
  
    Pnt->idInterp = int_S->adjid;
    Pnt->adjPnt[int_S->adjid].InterpFace = int_S->fid;
    Pnt->overlap = 0;
    
    free(interp);
    free(interpEq);
    return 1;
  }
  
  return 0;
}

/* is this an outerboundary point?
// yes if no points are in adjPnt or if we slightly move along the normal 
// reach nowhere.
// ->return value = 1 if outerbound; 0 otherwise. */
static int IsOutBndry(PointSet_T *const Pnt)
{
  /* confidently there is no point close by patch */
  if (Pnt->NadjPnt == 0 && !Pnt->Pnt->patch->innerB/* potentially it can be outerB. Note: innerB's are set a priory, but outerB's we don't know yet */)
    return 1;
  /* slightly move along normal and see if there is any other patches */
  else if (!Pnt->Pnt->patch->innerB)
  {
    Uint node = Pnt->Pnt->ind;
    Patch_T *const patch = Pnt->Pnt->patch;
    double *x = patch->node[node]->x;
    double *N1 = Pnt->Pnt->N;
    double eps = root_square(3,x,0)*ScaleFactor;
    double q[3];
    Uint ans;
    
    eps = GRT(eps,ScaleFactor) ? eps : ScaleFactor;
    
    if (Pnt->type == EDGE)
    {
      double Ntan[3];/* tanget vector for tilting */
      double Ntilt[3];/* tilting vector */
      const double T = 0.1;/* tilting angle is arctan(T) */
      double nrm;
      
      tangent(Pnt->Pnt,Ntan);
      nrm = root_square(3,Ntan,0);
      if (EQL(nrm,0))
        Error0("Normal vector is null!");
        
      /* make it unit */  
      Ntan[0] /= nrm;
      Ntan[1] /= nrm;
      Ntan[2] /= nrm;
      
      Ntilt[0] = N1[0]+T*Ntan[0];
      Ntilt[1] = N1[1]+T*Ntan[1];
      Ntilt[2] = N1[2]+T*Ntan[2];
      q[0] = x[0]+eps*Ntilt[0];
      q[1] = x[1]+eps*Ntilt[1];
      q[2] = x[2]+eps*Ntilt[2];
    }
    else if (Pnt->type == INNER)
    {
      q[0] = x[0]+eps*N1[0];
      q[1] = x[1]+eps*N1[1];
      q[2] = x[2]+eps*N1[2];
    }
    else
      Error0(NO_OPTION);
    
    Needle_T *needle = alloc_needle();
    needle->grid = patch->grid;
    needle_ex(needle,Pnt->Pnt->patch);
    needle->x = q;
    point_finder(needle);
    ans = needle->Nans;
    free_needle(needle);
    
    if (ans == 0)
      return 1;
  }
  
  return 0;
}

/* is this an inner boundary point?
// yes if no points are in adjPnt or if we slightly move along the normal 
// reach nowhere.
// ->return value = 1 if outerbound; 0 otherwise. */
static int IsInnBndry(PointSet_T *const Pnt)
{
  /* confidently there is no point close by patch */
  if (Pnt->NadjPnt == 0 && Pnt->Pnt->patch->innerB/* potentially it can be innerB */)
    return 1;
  /* slightly move along normal and see if there is any other patches */
  else if (Pnt->Pnt->patch->innerB)
  {
    Uint node = Pnt->Pnt->ind;
    Patch_T *const patch = Pnt->Pnt->patch;
    double *x = patch->node[node]->x;
    double *N1 = Pnt->Pnt->N;
    double eps = root_square(3,x,0)*ScaleFactor;
    double q[3];
    Uint ans;
    
    eps = GRT(eps,ScaleFactor) ? eps : ScaleFactor;
    
    if (Pnt->type == EDGE)
    {
      double Ntan[3];/* tanget vector for tilting */
      double Ntilt[3];/* tilting vector */
      const double T = 0.1;/* tilting angle is arctan(T) */
      double nrm;
      
      tangent(Pnt->Pnt,Ntan);
      nrm = root_square(3,Ntan,0);
      if (EQL(nrm,0))
        Error0("Normal vector is null!");
        
      /* make it unit */  
      Ntan[0] /= nrm;
      Ntan[1] /= nrm;
      Ntan[2] /= nrm;
      
      Ntilt[0] = N1[0]+T*Ntan[0];
      Ntilt[1] = N1[1]+T*Ntan[1];
      Ntilt[2] = N1[2]+T*Ntan[2];
      q[0] = x[0]+eps*Ntilt[0];
      q[1] = x[1]+eps*Ntilt[1];
      q[2] = x[2]+eps*Ntilt[2];
    }
    else if (Pnt->type == INNER)
    {
      q[0] = x[0]+eps*N1[0];
      q[1] = x[1]+eps*N1[1];
      q[2] = x[2]+eps*N1[2];
    }
    else
      Error0(NO_OPTION);
    
    Needle_T *needle = alloc_needle();
    needle->grid = patch->grid;
    needle_ex(needle,Pnt->Pnt->patch);
    needle->x = q;
    point_finder(needle);
    ans = needle->Nans;
    free_needle(needle);
    
    if (ans == 0)
      return 1;
  }
  
  return 0;
}

/* there are some cases -take place generally when interfaces of
// two patches reach outerbound - in which two normals are orthogonal 
// (or closely orthogonal) and those points can be flaged as outerbound.
// the algorithm goes like that:
// if each orthogonal normal be tilted toward the other one, and tip of
// this new vector couldn't be found in any other patches, flaged this
// as outerbound.
// ->return value: if their normals are orthogonal and surrounded by
// outerbound 1, otherwise 0. */
static int IsOrthOutBndry(PointSet_T *const Pnt)
{
  Uint i,j;
  
  if (Pnt->Pnt->patch->innerB)/* this patch has inner boundary, it's suspicious  */
    return 0;
  
  for (i = 0; i < Pnt->NadjPnt; i++)
  {
    if (Pnt->adjPnt[i].FaceFlg == 0) continue;
    
    for (j = 0; j < TOT_FACE; j++)
      if (Pnt->adjPnt[i].fs[j].OrthFlg == 1) 
        if (ReachBnd(Pnt,i,j))
        {
          Pnt->idOrth = i;
          return 1;
        }
  }  
  return 0;
}

/* there are some cases -take place generally when interfaces of
// two patches reach inner boundary - in which two normals are orthogonal 
// (or closely orthogonal) and those points can be flaged as inner boundary.
// the algorithm goes like that:
// if each orthogonal normal be tilted toward the other one, and tip of
// this new vector couldn't be found in any other patches, flaged this
// as innerbound.
// ->return value: if their normals are orthogonal and surrounded by
// innerbound 1, otherwise 0. */
static int IsOrthInnBndry(PointSet_T *const Pnt)
{
  Uint i,j;
  
  if (!Pnt->Pnt->patch->innerB)/* this patch dosen't have inner boundary, it's suspicious */
    return 0;
  
  for (i = 0; i < Pnt->NadjPnt; i++)
  {
    if (Pnt->adjPnt[i].FaceFlg == 0) continue;
    
    for (j = 0; j < TOT_FACE; j++)
      if (Pnt->adjPnt[i].fs[j].OrthFlg == 1) 
        if (ReachBnd(Pnt,i,j))
        {
          Pnt->idOrth = i;
          return 1;
        }
  }  
  return 0;
}

/* check if one tilts normal toward the adjPnt normal 
// there is no other patches, which means outerboundary or inner boundary.
// ->return value: 1 if boundary; 0 otherwise. */
static int ReachBnd(PointSet_T *const Pnt,const Uint p,const Uint f)
{
  Uint node = Pnt->Pnt->ind;
  Patch_T *const patch = Pnt->Pnt->patch;
  double *x = patch->node[node]->x;
  double *N1 = Pnt->Pnt->N;
  double *N2 = Pnt->adjPnt[p].fs[f].N2;
  double eps = root_square(3,x,0)*ScaleFactor;
  double q[3];
  Uint ans;
  
  eps = GRT(eps,ScaleFactor) ? eps : ScaleFactor;
  q[0] = x[0]+eps*N1[0]+ScaleFactor*N2[0];
  q[1] = x[1]+eps*N1[1]+ScaleFactor*N2[1];
  q[2] = x[2]+eps*N1[2]+ScaleFactor*N2[2];
  
  Needle_T *needle = alloc_needle();
  needle->grid = patch->grid;
  needle_ex(needle,Pnt->Pnt->patch);
  needle->x = q;
  point_finder(needle);
  ans = needle->Nans;
  free_needle(needle);
  
  if (ans == 0)
    return 1;
  
  return 0;
}

/* find the point p on the interface which has the fittest normal.
// ->return value: a pointer to point struct on interface
*/
static Point_T *get_p2(const PointSet_T *const Pnt)
{
  Uint f = Pnt->adjPnt[Pnt->idFit].CopyFace;
  Uint p = Pnt->adjPnt[Pnt->idFit].p;
  const Uint node = Pnt->adjPnt[Pnt->idFit].node;
  Patch_T *const patch = Pnt->Pnt->patch->grid->patch[p];
  Interface_T *const face = patch->interface[f];
  Uint i;
  
  FOR_ALL(i,face->point)
    if (face->point[i]->ind == node)
      return face->point[i];
  
  return 0;
}

/* if there is an adjacent point collocated with the point,
// and their interfaces are tangent correctly, i.e. their normals
// are in expected directions (N1dotN2 == -1). make sure
// that, among couple of possibilities, the one matches with the 
// previous subfaces, is chosen out.
// ->retrun value: 1 if fittest normal exists, 0 otherwise
*/
static int IsNormalFit(PointSet_T *const Pnt)
{
  Uint i,f;
  Point_T p1 = *Pnt->Pnt;/* copying */
  char *lead;
  Flag_T flg = NONE;
  
  p1.touch = 1;
  p1.copy  = 1;
  
  for (i = 0; i < Pnt->NadjPnt; i++)
  {
    if (Pnt->adjPnt[i].on_c == 1 && Pnt->adjPnt[i].FaceFlg == 1)
    {
      for (f = 0; f < TOT_FACE; f++)
        if (Pnt->adjPnt[i].fs[f].FitFlg == 1) 
        {
          flg = FOUND;
          Pnt->idFit = i;
          Pnt->adjPnt[i].CopyFace = f;
          
          p1.adjPatch = Pnt->adjPnt[i].p;
          p1.adjFace  = f;
          lead = inspect_flags(&p1);
          
          /* pick out the fittest one according to 
          // the previous made subfaces 
          */
          if (IsOnSubface(&p1,lead))
          {
            free(lead);
            return 1;
          }
          free(lead);
          
        }/* end of if (Pnt->adjPnt[i].fs[f].FitFlg == 1) */
    }
  }/* end of for (i = 0; i < Pnt->NadjPnt; i++) */
  
  /* if non of the found point has been on a subface, 
  // so it means it needs new subface. 
  */
  if (flg == FOUND)  return 1;
    
  return 0;
}
/* if none of the adjacent points are on an interface,
// thus, this point overlaps with other patch(es) 
// ->retrun value: 1 if it is overlap, 0 otherwise
*/
static int IsOverlap(PointSet_T *const Pnt)
{
  Uint i;
  
  for (i = 0; i < Pnt->NadjPnt; i++)
    if (Pnt->adjPnt[i].FaceFlg == 1) return 0;
  
  return 1;
}

/* finding the adjacent points of point pnt */
static void find_adjPnt(PointSet_T *const pnt)
{
  Uint ind = pnt->Pnt->ind;
  Uint *found;/* patches have been found */
  Uint nfp;/* number of found patches */
  double *x = pnt->Pnt->patch->node[ind]->x;
  Needle_T *needle = alloc_needle();
  needle->grid = pnt->Pnt->patch->grid;
  needle_ex(needle,pnt->Pnt->patch);
  needle->x = x;
  point_finder(needle); 
  found = needle->ans;  
  nfp = needle->Nans;
  add_adjPnt(pnt,found,nfp);
  
  /* freeing */    
  free_needle(needle);
}

/* initializing and finding inner points*/
static void init_Points(const Interface_T *const interface,PointSet_T ***const innP,PointSet_T ***const edgP)
{
  PointSet_T **pnt_in,**pnt_ed;/* inner and edge points */
  Uint N_in, N_ed;/* number of inner and edge points */
  const Uint f = interface->fn;
  const Uint *const n = interface->patch->n;
  Uint i,j,k,im,iM,jm,jM,km,kM;/* m for min and M for max */
  Uint sum,ed,in;
  
  N_in = NumPoint(interface,INNER);
  N_ed = NumPoint(interface,EDGE);
  
  set_min_max_sum(n,f,&im,&iM,&jm,&jM,&km,&kM,&sum);
  assert(N_in+N_ed == sum);
  
  ed = in = 0;
  pnt_in = 0;
  pnt_ed = 0;
  FOR_ijk(i,j,k,im,iM,jm,jM,km,kM)
  {
    Uint l = i_j_k_to_ijk(n,i,j,k);
    Uint p = L2(n,f,i,j,k);
    
    if (IsOnEdge(n,l))
    {
      alloc_PointSet(ed+1,&pnt_ed);
      pnt_ed[ed]->Pnt = interface->point[p];
      pnt_ed[ed]->type = EDGE;
      ed++;
    }
    else
    {
      alloc_PointSet(in+1,&pnt_in);
      pnt_in[in]->Pnt = interface->point[p];
      pnt_in[in]->type = INNER;
      in++;
    }  
  }
  
  (*innP) = pnt_in;
  (*edgP) = pnt_ed;
}

/* number of inner points on an interface
// ->return value: number of points on the given interface.
*/
static Uint NumPoint(const Interface_T *const interface,const enum Type type)
{
  Uint v = UINT_MAX;
  Uint *n = interface->patch->n;
  Uint f = interface->fn;
  
  if (type == INNER)
  {
    switch(f)
    {
      case I_0:
        v = (n[1]-2)*(n[2]-2);
        break;
      case I_n0:
        v = (n[1]-2)*(n[2]-2);
        break; 
      case J_0:
        v = (n[0]-2)*(n[2]-2);
        break; 
      case J_n1:
        v = (n[0]-2)*(n[2]-2);
        break; 
      case K_0:
        v = (n[0]-2)*(n[1]-2);
        break; 
      case K_n2:
        v = (n[0]-2)*(n[1]-2);
        break;
      default:
        Error0("There is not such interface.\n");
    }
  }
  else if (type == EDGE)
  {
    switch(f)
    {
      case I_0:
        v = 2*(n[1]+n[2]-2);
        break;
      case I_n0:
        v = 2*(n[1]+n[2]-2);
        break; 
      case J_0:
        v = 2*(n[0]+n[2]-2);
        break; 
      case J_n1:
        v = 2*(n[0]+n[2]-2);
        break; 
      case K_0:
        v = 2*(n[0]+n[1]-2);
        break; 
      case K_n2:
        v = 2*(n[0]+n[1]-2);
        break;
      default:
        Error0("There is not such interface.\n");
    }
  }
  else
    Error0("There is no such type.\n");
  
  return v;
}

/* find inner boundary for Cartesian type */
static void FindInnerB_Cartesian_coord(Patch_T *const patch)
{
  Interface_T **interface = patch->interface;
  Uint i,f;
  
  FOR_ALL(f,interface)
  {
    Point_T **point = interface[f]->point;
    FOR_ALL(i,point)
    {
      point[i]->innerB = 0;
    }
  }

}

/* find external faces for Cartesian type */
static void FindExterF_Cartesian_coord(Patch_T *const patch)
{
  Interface_T **interface = patch->interface;
  Uint i,f;
  
  FOR_ALL(f,interface)
  {
    Point_T **point = interface[f]->point;
    FOR_ALL(i,point)
    {
      point[i]->exterF = 1;
    }
  }
  
}

/* find inner boundary for cubed spherical type */
static void FindInnerB_CS_coord(Patch_T *const patch)
{
  Interface_T **interface = patch->interface;
  Uint i,f;
  
  if (patch->grid->kind == Grid_CubedSpherical_NSNS)
  {
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      FOR_ALL(i,point)
      {
        point[i]->innerB = 0;
      }
    }
  }
  else if (patch->grid->kind == Grid_CubedSpherical_BHNS)
  {
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      
      if (patch->innerB && f == K_0)
        FOR_ALL(i,point)
          point[i]->innerB = 1;
      else
        FOR_ALL(i,point)
          point[i]->innerB = 0;
    }
  }
  else if (patch->grid->kind == Grid_SplitCubedSpherical_BHNS)
  {
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      
      if (patch->innerB && f == K_0)
        FOR_ALL(i,point)
          point[i]->innerB = 1;
      else
        FOR_ALL(i,point)
          point[i]->innerB = 0;
    }
  }
  /* else if (strcmp_i(patch->grid->kind,"SNS_CubedSpherical+Box_grid"))
  {
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      FOR_ALL(i,point)
      {
        point[i]->innerB = 0;
      }
    }
  }
  else if (strcmp_i(patch->grid->kind,"SNS_CubedSpherical_grid"))
  {
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      FOR_ALL(i,point)
      {
        point[i]->innerB = 0;
      }
    }
  }
  else if (strcmp_i(patch->grid->kind,"SBH_CubedSpherical_grid"))
  {  
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      
      if (patch->innerB && f == K_0)
        FOR_ALL(i,point)
          point[i]->innerB = 1;
      else
        FOR_ALL(i,point)
          point[i]->innerB = 0;
    }
  }*/
  else
    Error0(INCOMPLETE_FUNC);
}

/* find external faces for cubed spherical type */
static void FindExterF_CS_coord(Patch_T *const patch)
{
  Interface_T **interface = patch->interface;
  Uint i,f;
  
  if (patch->grid->kind == Grid_CubedSpherical_NSNS)
  {
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      FOR_ALL(i,point)
      {
        point[i]->exterF = 1;
      }
    }
  }
  else if (patch->grid->kind == Grid_CubedSpherical_BHNS)
  {
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      FOR_ALL(i,point)
      {
        point[i]->exterF = 1;
      }
    }
  }
  else if (patch->grid->kind == Grid_SplitCubedSpherical_BHNS)
  {
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      FOR_ALL(i,point)
      {
        point[i]->exterF = 1;
      }
    }
  }
  /*else if (strcmp_i(patch->grid->kind,"SNS_CubedSpherical+Box_grid"))
  {
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      FOR_ALL(i,point)
      {
        point[i]->exterF = 1;
      }
    }
  }
  else if (strcmp_i(patch->grid->kind,"SNS_CubedSpherical_grid"))
  {
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      FOR_ALL(i,point)
      {
        point[i]->exterF = 1;
      }
    }
  }
  else if (strcmp_i(patch->grid->kind,"SBH_CubedSpherical_grid"))
  {
    FOR_ALL(f,interface)
    {
      Point_T **point = interface[f]->point;
      FOR_ALL(i,point)
      {
        point[i]->exterF = 1;
      }
    }
  } */
  else
    Error0(NO_OPTION);
}

/* filling point[?]->N */
static void 
fill_N
  (
    Patch_T *const patch
  )
{
  Interface_T **interface = patch->interface;
  Uint i,f;
  
  /* set normal vector for all face points */
  FOR_ALL(f,interface)
  {
    Point_T **point = interface[f]->point;
    
    FOR_ALL(i,point)
    {
      normal_vec(point[i]);
    }
  }
  /* set normal vector for only center (roughly) of each face */
  FOR_ALL(f,interface)
  {
    const Uint *n = interface[f]->patch->n;
    Uint ih = (n[0]-1)/2;
    Uint jh = (n[1]-1)/2;
    Uint kh = (n[2]-1)/2;
    Uint ijk= UINT_MAX;
    
    switch(f)
    {
      case I_0:
        ijk = i_j_k_to_ijk(n,0,jh,kh);
        break;
      case I_n0:
        ijk = i_j_k_to_ijk(n,n[0]-1,jh,kh);
        break; 
      case J_0:
        ijk = i_j_k_to_ijk(n,ih,0,kh);
        break; 
      case J_n1:
        ijk = i_j_k_to_ijk(n,ih,n[1]-1,kh);
        break; 
      case K_0:
        ijk = i_j_k_to_ijk(n,ih,jh,0);
        break; 
      case K_n2:
        ijk = i_j_k_to_ijk(n,ih,jh,n[2]-1);
        break;
      default:
        Error0("There is not such interface.\n");
    }
    Point_T point[1] = {0};
    point->ind   = ijk;
    point->patch = interface[f]->patch;
    point->face  = f;
    normal_vec(point);
    /* center normal */
    interface[f]->centerN[0] = point->N[0];
    interface[f]->centerN[1] = point->N[1];
    interface[f]->centerN[2] = point->N[2];
    /* center x coords */
    interface[f]->centerx[0] = patch->node[ijk]->x[0];
    interface[f]->centerx[1] = patch->node[ijk]->x[1];
    interface[f]->centerx[2] = patch->node[ijk]->x[2];
  }
}

/* 
// filling  basic elements:
// point[?]->ind, point[?]->face and point[?]->patch and
// also realize those nodes reach boundary of patch
*/
static void fill_basics(Patch_T *const patch)
{
  Uint *n = patch->n;
  Uint f;
  
  for (f = I_0; f < TOT_FACE; f++)
  {
    Point_T **point;
    Uint i,j,k,l,im,iM,jm,jM,km,kM,sum;/* m for min and M for max */
    Uint t = 0;/* test purposes */
    
    set_min_max_sum(n,f,&im,&iM,&jm,&jM,&km,&kM,&sum);
    patch->interface[f]->point = alloc_point(sum);
    patch->interface[f]->patch = patch;
    patch->interface[f]->fn = f;
    patch->interface[f]->np = sum;
    point = patch->interface[f]->point;
    
    t = 0;
    FOR_ijk(i,j,k,im,iM,jm,jM,km,kM)
    {
      Uint p = L2(n,f,i,j,k);
      point[p]->ind = l = i_j_k_to_ijk(n,i,j,k);
      point[p]->face = f;
      point[p]->patch = patch;
      point[p]->exterF= 1;
      if (IsOnEdge(n,l))
        point[p]->IsOnEdge = 1;
      t++;
    }
    assert(t == sum);
  }
}

/* 2d linear index format for interface points.
// L = k+n2*(j+n1*i) so depends on interface some of i,j,j,n2 and n1
// are zero or supremum as follows:
// ->return value: 2d index.
*/
static Uint L2(const Uint *const n,const Uint f, const Uint i, const Uint j, const Uint k)
{
  Uint v = UINT_MAX;
  
  switch(f)
  {
    case I_0:/* i = 0 */
      v = k+n[2]*j;
      break;
    case I_n0:/* i = n[0]-1 */
      v = k+n[2]*j;
      break;
    case J_0:/* j = 0 */
      v = k+n[2]*i;
      break;
    case J_n1:/* j = n[1]-1 */
      v = k+n[2]*i;
      break;
    case K_0:/* k = 0 */
      v = j+n[1]*i;
      break;
    case K_n2:/* k = n[2]-1 */
      v = j+n[1]*i;
      break;
    default:
      Error0("Exceeding from total number of face.\n");
      break; 
  }
  
  return v;
}

/* setting min and max of total number of point on an interface.
// the purpose it serves is making loop over interface point easier.
// note: since in FOR_ijk the condition for loop is less than sign <
// then the maximum value must be the supremum; but the slice 
// should be at the exact value.
*/
static void set_min_max_sum(const Uint *const n,const Uint f,Uint *const im,Uint *const iM,Uint *const jm,Uint *const jM,Uint *const km,Uint *const kM,Uint *const sum)
{
  switch(f)
  {
    case I_0:
      *im = 0;
      *iM = 1;
      *jm = 0;
      *jM = n[1];
      *km = 0;
      *kM = n[2];
      *sum = n[1]*n[2];
      break;
    case I_n0:
      *im = n[0]-1;
      *iM = n[0];
      *jm = 0;
      *jM = n[1];
      *km = 0;
      *kM = n[2];
      *sum = n[1]*n[2];
      break;
    case J_0:
      *im = 0;
      *iM = n[0];
      *jm = 0;
      *jM = 1;
      *km = 0;
      *kM = n[2];
      *sum = n[0]*n[2];
      break;    
    case J_n1:
      *im = 0;
      *iM = n[0];
      *jm = n[1]-1;
      *jM = n[1];
      *km = 0;
      *kM = n[2];
      *sum = n[0]*n[2];
      break;    
    case K_0:
      *im = 0;
      *iM = n[0];
      *jm = 0;
      *jM = n[1];
      *km = 0;
      *kM = 1;
      *sum = n[0]*n[1]; 
      break;    
    case K_n2:
      *im = 0;
      *iM = n[0];
      *jm = 0;
      *jM = n[1];
      *km = n[2]-1;
      *kM = n[2];
      *sum = n[0]*n[1];
      break;
    default:
      Error0("Exceeding from total number of face.\n");
      break;
  }
}

/* normal vector at interface of patch, points "OUTWARD" and "NORMALIZED" and
// in calculated in Cartesian coords.
// the normal vector is written in N at point structure;
// note: the members in point struct that must be filled before
// passing to this function are: "ind","patch","face".
// -> return value: a pointer to N.
*/
double *normal_vec(Point_T *const point)
{
  if (point->patch->coordsys == Cartesian)
  {
    normal_vec_Cartesian_coord(point);
  }
  else if (point->patch->coordsys == CubedSpherical)
  {
    normal_vec_CS_coord(point);
  }
  else
    Error0("No Normal defined for such coordinate yet!\n");
  
  return point->N;
}

/* make sure the normal is outward. 
// it checks if q = x+eps*N won't be found in the patch 
// where the normal on its interface is being found.
// if it finds then it will multiplies N by -1. */
void make_normal_outward(Point_T *const point)
{
  const double *const x = point->patch->node[point->ind]->x;
  
  double eps = root_square(3,x,0)*ScaleFactor;
  Needle_T *needle = alloc_needle();
  const double *N = point->N;
  double q[3];/* q = pnt+eps*N */
  
  eps = GRT(eps,ScaleFactor) ? eps : ScaleFactor;
  q[0] = x[0]+eps*N[0];
  q[1] = x[1]+eps*N[1];
  q[2] = x[2]+eps*N[2];
  needle->grid = point->patch->grid;
  needle->x = q;
  needle_in(needle,point->patch);
  point_finder(needle);
  
  if (needle->Nans > 0)
  {
    point->N[0] *= -1;
    point->N[1] *= -1;
    point->N[2] *= -1;
  }
  
  free_needle(needle);
}

/* finding normal for cubed spherical coord */
static void normal_vec_CS_coord(Point_T *const point)
{
  const Uint p = point->ind;
  Patch_T *const patch = point->patch;
  double N;
  
  switch(point->face)
  {
    case I_0:
      point->N[0] = -dq2_dq1(patch,_a_,_x_,p);
      point->N[1] = -dq2_dq1(patch,_a_,_y_,p);
      point->N[2] = -dq2_dq1(patch,_a_,_z_,p);
    break;
    case I_n0:
      point->N[0] = dq2_dq1(patch,_a_,_x_,p);
      point->N[1] = dq2_dq1(patch,_a_,_y_,p);
      point->N[2] = dq2_dq1(patch,_a_,_z_,p);
    break;
    case J_0:
      point->N[0] = -dq2_dq1(patch,_b_,_x_,p);
      point->N[1] = -dq2_dq1(patch,_b_,_y_,p);
      point->N[2] = -dq2_dq1(patch,_b_,_z_,p);
    break;
    case J_n1:
      point->N[0] = dq2_dq1(patch,_b_,_x_,p);
      point->N[1] = dq2_dq1(patch,_b_,_y_,p);
      point->N[2] = dq2_dq1(patch,_b_,_z_,p);
    break;
    case K_0:
      point->N[0] = -dq2_dq1(patch,_c_,_x_,p);
      point->N[1] = -dq2_dq1(patch,_c_,_y_,p);
      point->N[2] = -dq2_dq1(patch,_c_,_z_,p);
    break;
    case K_n2:
      point->N[0] = dq2_dq1(patch,_c_,_x_,p);
      point->N[1] = dq2_dq1(patch,_c_,_y_,p);
      point->N[2] = dq2_dq1(patch,_c_,_z_,p);
    break;
    default:
      Error0("There is no such face.\n");
  }
  
  N = root_square(3,point->N,0);
  if (EQL(N,0))
    Error0("Normal vector is null!");
  
  /* make it unit */  
  point->N[0] /= N;
  point->N[1] /= N;
  point->N[2] /= N;
}

/* finding normal for Cartesian coord */
static void normal_vec_Cartesian_coord(Point_T *const point)
{
  switch(point->face)
  {
    case I_0:
      point->N[0] = -1;
      point->N[1] = 0;
      point->N[2] = 0;
      break;
    case I_n0:
      point->N[0] = 1;
      point->N[1] = 0;
      point->N[2] = 0;
      break;
    case J_0:
      point->N[0] = 0;
      point->N[1] = -1;
      point->N[2] = 0;
      break;
    case J_n1:
      point->N[0] = 0;
      point->N[1] = 1;
      point->N[2] = 0;
      break;
    case K_0:
      point->N[0] = 0;
      point->N[1] = 0;
      point->N[2] = -1;
      break;
    case K_n2:
      point->N[0] = 0;
      point->N[1] = 0;
      point->N[2] = 1;
      break;
    default:
      Error0("There is no such face.\n");
  }
  
  
}

/* reallocation mem for PointSet_T strcut with one extera block
// to put the last pointer to null.
*/
static void alloc_PointSet(const Uint N,PointSet_T ***const pnt)
{
  assert(pnt);
  
  (*pnt) = realloc((*pnt),(N+1)*sizeof(*(*pnt)));
  IsNull(*pnt);
  
  (*pnt)[N-1] = calloc(1,sizeof(*(*pnt)[N-1]));
  IsNull((*pnt)[N-1]);
  
  (*pnt)[N] = 0;
}

/* free mem for PointSet_T strcut */
static void free_PointSet(PointSet_T **const pnt)
{
  if(!pnt) return;
  
  Uint i;
  FOR_ALL (i,pnt)
  {
    if (pnt[i]->adjPnt)
    {
      free(pnt[i]->adjPnt);
    }
    
    free(pnt[i]);
  }
  
  free(pnt);
}

/* adding all of patches found in to point->adjPnt;
// and find if there is any node in adjacent patch correspond to
// this point pnt. moreover, study if this pnt located on some interfaces
// of adjacent patch, if yes put adjPnt->fs[face_nun].on = 1 and if not
// it is 0.
// note: it won't add any thing if the patch is already exists.
*/
static void add_adjPnt(PointSet_T *const pnt,const Uint *const p, const Uint np)
{
  Uint i,j;
  for (i = 0; i < np; i++)
  {
    for(j = 0; j < pnt->NadjPnt; j++)
      if(p[i] == pnt->adjPnt[j].p) return;
        
    /* realloc memory for new found point */
    pnt->adjPnt = 
      realloc(pnt->adjPnt,(pnt->NadjPnt+1)*sizeof(*pnt->adjPnt));
    IsNull(pnt->adjPnt);
    
    /* initializing */
    AdjPoint_T *const adjPnt = &pnt->adjPnt[pnt->NadjPnt];
    adjPnt->p = p[i];/* fill p */
    adjPnt->FaceFlg = 0;
    adjPnt->on_c = 0;
    
    fill_adjPnt(pnt,pnt->NadjPnt);/* filling elemetns of adjPnt struct */
    
    pnt->NadjPnt++;
  }
}

/* study:
// 1. if this pnt collocated with some node of adjacent patch
// 2. if this pnt located on some interfaces of adjacent patch
// 3. find adjPnt normal vector and its dot product with the pnt
// the last one gets used for realizing the geometry of neighbor of pnt
*/
static void fill_adjPnt(PointSet_T *const pnt,const Uint N)
{
  AdjPoint_T *const adjPnt = &pnt->adjPnt[N];
  Uint ind = pnt->Pnt->ind;
  const double *const x = pnt->Pnt->patch->node[ind]->x;
  Uint f[TOT_FACE];
  const Grid_T *const grid = pnt->Pnt->patch->grid;
  Flag_T flg;
  Uint ind2;
  
  /* find if it is on a node (collocation) */
  ind2 = find_node(x,grid->patch[adjPnt->p],&flg);
  if (flg == FOUND)
  {
    adjPnt->node = ind2;
    adjPnt->on_c = 1;
  }
  
  /* find if it is on a face or not */  
  if (IsOnFace(x,grid->patch[adjPnt->p],f,PrecisionFactor))
  {
    Point_T po;
    double *N2;
    Uint i;
    
    po.patch = grid->patch[adjPnt->p];
    adjPnt->FaceFlg = 1;
    
    for (i = 0; i < TOT_FACE; i++)
    {
      adjPnt->fs[i].FitFlg  = 0;
      adjPnt->fs[i].OrthFlg = 0;
      
      if (f[i])
      {
        po.ind = node_onFace(x,i,po.patch);
        po.face = i;
        N2 = normal_vec(&po);
        
        adjPnt->fs[i].on_f = 1;
        adjPnt->fs[i].N2[0] = N2[0];
        adjPnt->fs[i].N2[1] = N2[1];
        adjPnt->fs[i].N2[2] = N2[2];
        adjPnt->fs[i].N1dotN2 = dot(3,pnt->Pnt->N,N2);
        
        if (EQL(adjPnt->fs[i].N1dotN2,-1)) 
          adjPnt->fs[i].FitFlg  = 1;
        else if (EQL(adjPnt->fs[i].N1dotN2,0)) 
          adjPnt->fs[i].OrthFlg = 1;
      }
      else
        adjPnt->fs[i].on_f = 0;
    }
  }/* if (IsOnFace(...)) */
  
}

/* find the closest inner point (not on edge) to this point and then 
// put the vector made by the difference of these two points in N.
// this is the approximate tangent vector. note, the tangent
// vector won't be "NORMALIZED". I want to keep the magnetitue to be of
// the order of grid size. 
// NOTE: the given pnt must be on edge otherwise N is 0. */
void tangent(const Point_T *const pnt,double *const N)
{
  const Uint *const n = pnt->patch->n;
  const Uint ind = pnt->ind;
  const Uint f = pnt->face;
  double *x = pnt->patch->node[ind]->x;
  double *y;
  double s_ds = DBL_MAX;/* smallest distance */
  Uint s_in = UINT_MAX;/* index referring to point with smallest distance */
  Uint i,j,k,im,iM,jm,jM,km,kM,sum;/* m for min and M for max */
  Flag_T flg = NONE;
  
  set_min_max_sum(n,f,&im,&iM,&jm,&jM,&km,&kM,&sum);
  
  FOR_ijk(i,j,k,im,iM,jm,jM,km,kM)
  {
    Uint l = i_j_k_to_ijk(n,i,j,k);
    double nrm;
    
    if (!IsOnEdge(n,l))
    {
      y = pnt->patch->node[l]->x;
      nrm = root_square(3,x,y);
      
      if (LSS(nrm,s_ds)) 
      {
        flg = FOUND;
        s_in = l;
      }
    }
  }
  
  if (flg == NONE)
    Error0("Tangent vector could not be found.\n");
  
  /* note: it must be y-x to tilt toward interface not out of it */  
  y = pnt->patch->node[s_in]->x;
  N[0] = y[0]-x[0];
  N[1] = y[1]-x[1];
  N[2] = y[2]-x[2];
}

/* getting interface, figures what is the constant coordinates
// on this interface.
*/
static void set_sameXYZ(Point_T *const p,const Uint f)
{
  switch(f)
  {
    case I_0:
      p->sameX = 1;
      p->sameY = 0;
      p->sameZ = 0;
      break;
    case I_n0:
      p->sameX = 1;
      p->sameY = 0;
      p->sameZ = 0;
      break; 
    case J_0:
      p->sameX = 0;
      p->sameY = 1;
      p->sameZ = 0;
      break; 
    case J_n1:
      p->sameX = 0;
      p->sameY = 1;
      p->sameZ = 0;
      break; 
    case K_0:
      p->sameX = 0;
      p->sameY = 0;
      p->sameZ = 1;
      break; 
    case K_n2:
      p->sameX = 0;
      p->sameY = 0;
      p->sameZ = 1;
      break;
    default:
      Error0("There is not such interface.\n");
  }
}

/* see if there is a subface with these lead as its flags_str 
// ->return value: 1 if found, 0 otherwise.
*/
static int IsOnSubface(const Point_T *const pnt, const char *const lead)
{
  Interface_T *const face = pnt->patch->interface[pnt->face];
  Uint sf;
  
  /* compare lead with the previous sub-faces */
  for (sf = 0; sf < face->ns; ++sf)
  {
    if (strcmp_i(face->subface[sf]->flags_str,lead))
      return 1;
  }

  return 0;
}

/* memory alloc nodes */
void alloc_nodes(Grid_T *const grid)
{
  Uint *n;
  Uint i;
  
  FOR_ALL(i,grid->patch)
  {
    Uint j,U;
    Node_T **node;
    
    n = grid->patch[i]->n;
    grid->patch[i]->node = 
      calloc((n[0]*n[1]*n[2]+1),sizeof(*grid->patch[i]->node));
    IsNull(grid->patch[i]->node);
    
    node = grid->patch[i]->node;
    node[n[0]*n[1]*n[2]] = 0;
    
    U = n[0]*n[1]*n[2];
    for (j = 0; j < U; j++)
    {
      node[j] = calloc(1,sizeof(*node[j]));
      IsNull(node[j]);
    }
    
  }
}

/* memory allocation for interface struct */
void alloc_interface(Patch_T *const patch)
{
  Uint i;
  assert(patch);
  
  patch->interface = calloc(NFaces+1,sizeof(*patch->interface));
  IsNull(patch->interface);
  
  for (i = 0; i < NFaces; i++)
  {
    patch->interface[i] = calloc(1,sizeof(*patch->interface[i]));
    IsNull(patch->interface[i]);
  }
}

/*
// memory allocation for point struct;
// s is the number of point which is demanded 
// ->return value: pointer to new allocated memory
*/
void *alloc_point(const Uint s)
{
  Point_T **point;
  Uint i;
  
  point = calloc(s+1,sizeof(*point));
  IsNull(point);
  
  for (i = 0; i < s; i++)
  {
    point[i] = calloc(1,sizeof(*point[i]));
    IsNull(point[i]);
  }
  
  return point;
}

/* feeing memory of Point_T inside grid */
void free_points(Grid_T *const grid)
{
  Uint pa;
  
  FOR_ALL(pa,grid->patch)
  {
    Interface_T **face = grid->patch[pa]->interface;
    Uint f;
    
    FOR_ALL(f,face)
    {
      free_2d_mem(face[f]->point,face[f]->np);
      face[f]->np = 0;
      face[f]->point = 0;
    }

  }
}

/* free thoroughly patch->interface */
void free_patch_interface(Patch_T *const patch)
{
  Interface_T **face = patch->interface;
  SubFace_T **subface = 0;
  Uint f,i;
  
  if (face)
  FOR_ALL(f,face)
  {
    /* free point */
    free_2d_mem(face[f]->point,face[f]->np);
    face[f]->np = 0;
    
    /* free subface */
    subface = face[f]->subface;
    for (i = 0; i < face[f]->ns; ++i)
    {
      Free(subface[i]->flags_str);
      Free(subface[i]->id);
      Free(subface[i]->adjid);
      Free(subface[i]);
    }
    Free(subface);
    
    /* free point set */
    Free(face[f]->innerP);
    Free(face[f]->edgeP);
    
    /* free face */
    Free(face[f]);
  }
  Free(face);

  patch->interface = 0;
}

/* check if all of houseKs have been marked */
static void check_houseK(Patch_T *const patch)
{
  Interface_T **const interface = patch->interface;
  const Uint nf = countf(interface);
  Node_T *node;
  Uint i,f;
  
  for (f = 0; f < nf; f++)
    for (i = 0; i < interface[f]->np; i++)
      if (interface[f]->point[i]->houseK == 0)
      {
        node = patch->node[interface[f]->point[i]->ind];
        double *x = node->x;
        fprintf(stderr,"This point(%f,%f,%f) has not been found.\n",
                        x[0],x[1],x[2]);
        Error0("Incomplete function.\n");
      }
}

/* setting all of houseK (house keeping) flags 
// in Point_T to zero for given patch */
static void flush_houseK(Patch_T *const patch)
{
  Interface_T **const interface = patch->interface;
  const Uint nf = countf(interface);
  Uint i,f;
  
  for (f = 0; f < nf; f++)
    for (i = 0; i < interface[f]->np; i++)
      interface[f]->point[i]->houseK = 0;
}


/* displace the center point along the normal and find where the new
// point takes place to find the adjacent patch and adjacent face.
// also it finds innerB and outerB and set their flags.
// this function has many assumptions which is only for 
// split cubed spherical. */
static void 
find_adjacent_scs
  (
  Grid_T *const grid,
  Patch_T *const patch,
  Uint *const point_flag
  )
{
  Interface_T **const interface = patch->interface;
  const double idealN1N2 = -1;
  double q[3] = {0};/* q = x+eps*N */
  double eps,N1dotN2;
  double min = DBL_MAX;
  char s[999] = {'\0'};
  int f;
  Uint i,p;
  
  for (f = (int)NFaces-1; f >= 0; --f)
  {
    /* FIRST: find the adjacent patch to the center. */
    Needle_T *needle0 = 0;
    const double *centerN, *centerx;
    Uint Nfound0;
    
    centerN = interface[f]->centerN;
    centerx = interface[f]->centerx;
    
    /* try without normal */
    q[0] = centerx[0];
    q[1] = centerx[1];
    q[2] = centerx[2];
    needle0 = alloc_needle();
    needle0->precision_factor = PrecisionFactor;
    needle0->grid = grid;
    needle_ex(needle0,patch);
    needle0->x = q;
    point_finder(needle0);
    Nfound0 = needle0->Nans;
    /* try again with normal */
    if (!Nfound0)
    {
      eps = root_square(3,centerx,0)*ScaleFactor;
      eps = GRT(eps,ScaleFactor) ? eps : ScaleFactor;
      q[0] = centerx[0]+eps*centerN[0];
      q[1] = centerx[1]+eps*centerN[1];
      q[2] = centerx[2]+eps*centerN[2];
      free_needle(needle0);
      needle0 = alloc_needle();
      needle0->precision_factor = PrecisionFactor;
      needle0->grid = grid;
      needle_ex(needle0,patch);
      needle0->x = q;
      point_finder(needle0);
      Nfound0 = needle0->Nans;
    }
    
    /* if still no patch found 3 possibilities = {innerB,outerB,gap} */
    if (!Nfound0)
    {
      /* make sure there is really no other patches */
      Uint Nfound1;
        
      eps = root_square(3,centerx,0)*ScaleFactor;
      eps = GRT(eps,ScaleFactor) ? eps : ScaleFactor;
      eps += Huge_Eps;
      q[0] = centerx[0]+eps*centerN[0];
      q[1] = centerx[1]+eps*centerN[1];
      q[2] = centerx[2]+eps*centerN[2];
  
      Needle_T *needle1 = alloc_needle();
      needle1->precision_factor = PrecisionFactor;
      needle1->grid = grid;
      needle_ex(needle1,patch);
      needle1->x = q;
      point_finder(needle1);
      Nfound1 = needle1->Nans;
      free_needle(needle1);
      
      if (!Nfound1)/* there is really no other patch => outer/innter-boundary */
      {
        Point_T **const pnt = interface[f]->point;
        
        if (patch->innerB)/* we expect this */
        {
          /* make all of this interface point in one innerB subface */
          for (p = 0; p < interface[f]->np; ++p)
          {
            pnt[p]->innerB = 1;
            add_to_subface_scs(pnt[p]);
            point_flag[pnt[p]->ind] = 1;
          }
          interface[f]->innerB = 1;
        } 
        else/* it must be outerB */
        {
          patch->outerB  = 1;
          /* make all of this interface point in one innerB subface */
          for (p = 0; p < interface[f]->np; ++p)
          {
            pnt[p]->outerB = 1;
            add_to_subface_scs(pnt[p]);
            point_flag[pnt[p]->ind] = 1;
          }
          interface[f]->outerB = 1;
        }
      }
      else
      {
        err_spr(s,centerx);
        Errors(" Unexpected gap between patches!\n %s\n",s);
      }
    }/* if (!Nfound0) */
    else/* find adjacent subfaces for this interface */
    {
      Point_T **const pnt = interface[f]->point;
      const Patch_T *center_adjpatch = 0;
      Uint center_adjface = UINT_MAX;
      
      min = DBL_MAX;
      
      /* find the best adjacet adjPatch for center */
      for (i = 0; i < Nfound0; ++i)
      {
        const Patch_T *adjpatch = grid->patch[needle0->ans[i]];
        Uint adjfn;
        
        /* find the best adjface */
        for (adjfn = 0; adjfn < NFaces; ++adjfn)
        {
          N1dotN2 = 
            dot(3,adjpatch->interface[adjfn]->centerN,centerN);
          double fabsdif = fabs(N1dotN2-idealN1N2);
          if (min > fabsdif)
          {
            center_adjface  = adjfn;
            center_adjpatch = adjpatch;
            min             = fabsdif;
          }
        }
      }
      if (min > ScaleFactor)
      {
        warn_normal_spr(s,patch,center_adjpatch,center_adjface);
        printf("%s\n",s);
        fflush(stdout);
      }
      
      /* havining found the potential face to this patch 
      // find point to the pertinent adjacent. */
      for (p = 0; p < interface[f]->np; ++p)
      {
        if (!point_flag[pnt[p]->ind])
        {
          Uint ind = pnt[p]->ind;
          Uint adjind;
          const double *pnt_x = patch->node[ind]->x;
          double adjX[3] = {0};
          Flag_T iscp = NONE;
          
          /* first try to find this point on center adjpatch */
          if (X_of_x_precision(adjX,pnt_x,center_adjpatch,PrecisionFactor))
          {
            pnt[p]->touch    = 1;
            pnt[p]->adjPatch = center_adjpatch->pn;
            pnt[p]->adjFace  = center_adjface;
            
            /* if this is a copy? */
            adjind = find_node(pnt_x,center_adjpatch,&iscp);
            if (iscp == FOUND)
            {
              Uint adjp = pnt[p]->adjPatch;
              Uint adjf = pnt[p]->adjFace;
              Uint nadjp = grid->patch[adjp]->interface[adjf]->np;
              Point_T **adjpnt = 
                grid->patch[adjp]->interface[adjf]->point;
              Uint ip;
              
              pnt[p]->copy = 1;
              pnt[p]->adjInd = adjind;
              pnt[p]->adjIndF= UINT_MAX;/* catch error */
              /* find index of adjacent point on face. */
              for (ip = 0; ip < nadjp; ++ip)
              {
                 if (adjpnt[ip]->ind == adjind)
                 {
                   pnt[p]->adjIndF = ip;
                   break;
                 }
              }
              assert(pnt[p]->adjIndF != UINT_MAX);
            }
            else
            {
              /* if this is a interpolation => copy = 0*/
              set_sameXYZ(pnt[p],center_adjface);
            }
          }
          /* it means not all points on this surface are on a same
          // adjacent patch, for example ,filling_box reaches few patches.
          // thus, we should find their adjacent patch. */
          else
          {
            Needle_T *needle2 = 0;
            Uint Nfound2 = 0;
            
            eps = root_square(3,pnt_x,0)*ScaleFactor;
            eps = GRT(eps,ScaleFactor) ? eps : ScaleFactor;
            
            if (pnt[p]->IsOnEdge)
            {
              q[0] = pnt_x[0];
              q[1] = pnt_x[1];
              q[2] = pnt_x[2];
              needle2 = alloc_needle();
              needle2->precision_factor = PrecisionFactor;
              needle2->grid = grid;
              needle_ex(needle2,patch);
              needle_ex(needle2,center_adjpatch);
              needle2->x = q;
              point_finder(needle2);
              Nfound2 = needle2->Nans;
            }
            else
            {
              q[0] = pnt_x[0];
              q[1] = pnt_x[1];
              q[2] = pnt_x[2];
              needle2 = alloc_needle();
              needle2->precision_factor = PrecisionFactor;
              needle2->grid = grid;
              needle_ex(needle2,patch);
              needle_ex(needle2,center_adjpatch);
              needle2->x = q;
              point_finder(needle2);
              Nfound2 = needle2->Nans;
              /* try with normal */
              if (!Nfound2)
              {
                free_needle(needle2);
                q[0] = pnt_x[0]+eps*pnt[p]->N[0];
                q[1] = pnt_x[1]+eps*pnt[p]->N[1];
                q[2] = pnt_x[2]+eps*pnt[p]->N[2];
                needle2 = alloc_needle();
                needle2->precision_factor = PrecisionFactor;
                needle2->grid = grid;
                needle_ex(needle2,patch);
                needle_ex(needle2,center_adjpatch);
                needle2->x = q;
                point_finder(needle2);
                Nfound2 = needle2->Nans;
              }
            }
            
            if (!Nfound2)
            {
            
              /* test why could not find the adjacent patch */
              double tt_X[3];
              int tt_ret = X_of_x_precision(tt_X,pnt_x,center_adjpatch,PrecisionFactor);
              printf("%s: [%g,%g]x[%g,%g]x[%g,%g]\n"
                     "%d <- X_of_x((%0.15f,%0.15f,%0.15f))\n",
              center_adjpatch->name,
              center_adjpatch->min[0],center_adjpatch->max[0],
              center_adjpatch->min[1],center_adjpatch->max[1],
              center_adjpatch->min[2],center_adjpatch->max[2],
              tt_ret,tt_X[0],tt_X[1],tt_X[2]);
              fflush(stdout);
              err_spr_adj(s,pnt_x,center_adjpatch,center_adjface);
              Errors(" Point not found!\n %s\n",s);
            }
            
            /* first find the best adjface */
            const Patch_T *pnt_adjpatch = 0;
            Uint pnt_adjface = UINT_MAX;
            
            min = DBL_MAX;
            for (i = 0; i < Nfound2; ++i)
            {
              const Patch_T *adjpatch = grid->patch[needle2->ans[i]];
              Uint onface[NFaces];
              Uint adjfn;
              
              /* find the best adjface if it is on face at all */
              if(IsOnFace(pnt_x,adjpatch,onface,PrecisionFactor))
              for (adjfn = 0; adjfn < NFaces; ++adjfn)
              {
                if (!onface[adjfn])
                  continue;
                  
                N1dotN2 = 
                  dot(3,centerN,adjpatch->interface[adjfn]->centerN);
                double fabsdif = fabs(N1dotN2-idealN1N2);
                
                if (min > fabsdif)
                {
                  pnt_adjface  = adjfn;
                  pnt_adjpatch = adjpatch;
                  min          = fabsdif;
                }
              }
            }
            if (min > ScaleFactor)
            {
              warn_normal_spr(s,patch,pnt_adjpatch,pnt_adjface);
              printf("%s\n",s);
              fflush(stdout);
            }
            /* find this point on adjpatch */
            if (X_of_x_precision(adjX,pnt_x,pnt_adjpatch,PrecisionFactor))
            {
              pnt[p]->touch    = 1;
              pnt[p]->adjPatch = pnt_adjpatch->pn;
              pnt[p]->adjFace  = pnt_adjface;
              
              /* if this is a copy? */
              adjind = find_node(pnt_x,pnt_adjpatch,&iscp);
              if (iscp == FOUND)
              {
                Uint adjp = pnt[p]->adjPatch;
                Uint adjf = pnt[p]->adjFace;
                Uint nadjp = grid->patch[adjp]->interface[adjf]->np;
                Point_T **adjpnt = 
                  grid->patch[adjp]->interface[adjf]->point;
                Uint ip;
                
                pnt[p]->copy = 1;
                pnt[p]->adjInd = adjind;
                pnt[p]->adjIndF= UINT_MAX;/* catch error */
                /* find index of adjacent point on face. */
                for (ip = 0; ip < nadjp; ++ip)
                {
                   if (adjpnt[ip]->ind == adjind)
                   {
                     pnt[p]->adjIndF = ip;
                     break;
                   }
                }
                
                /* could not find the point! why? */
                assert(pnt[p]->adjIndF != UINT_MAX);
              }
              else
              {
                /* if this is a interpolation => copy = 0*/
                set_sameXYZ(pnt[p],pnt_adjface);
              }
            }
            else
            {
              /* test why could not find the adjacent patch */
              double tt_X[3];
              int tt_ret = X_of_x_precision(tt_X,pnt_x,center_adjpatch,PrecisionFactor);
              printf("%s: [%g,%g]x[%g,%g]x[%g,%g]\n"
                     "%d <- X_of_x((%0.15f,%0.15f,%0.15f))\n",
              center_adjpatch->name,
              center_adjpatch->min[0],center_adjpatch->max[0],
              center_adjpatch->min[1],center_adjpatch->max[1],
              center_adjpatch->min[2],center_adjpatch->max[2],
              tt_ret,tt_X[0],tt_X[1],tt_X[2]);
              fflush(stdout);
              err_spr_adj(s,pnt_x,center_adjpatch,center_adjface);
              Errors(" Point not found!\n %s\n",s);
            }
            free_needle(needle2);
        
          }
          /* NOTE: for interpolation and copy points the houseK is not
          // set. */
          point_flag[pnt[p]->ind] = 1;
          
        }/* if (!point_flag[pnt[p]->ind]) */
        
        else/* means point_flag[pnt[p]->ind] is set but houseK might no
            // set, for instance an edge point might set from one side
            // but not from the other side. */
        {
          pnt[p]->houseK = 1;
        }
      }/* for (p = 0; p < interface[f]->np; ++p) */
    }
    free_needle(needle0);
  }
}

/* having found the adjacent points now set the correspondence
// subfaces consistently. */
static void set_subfaces_scs(Grid_T *const grid,Patch_T *const patch)
{
  Interface_T **const interface = patch->interface;
  Uint f;
  
  for (f = 0; f < NFaces; ++f)
  {
    Uint p;
    for (p = 0; p < interface[f]->np; ++p)
    {
      Point_T **const pnt = interface[f]->point;
      
      if (!pnt[p]->houseK)
      {
        add_to_subface_scs(pnt[p]);
        if (pnt[p]->copy)/* copy */
        {
          Uint adjf = pnt[p]->adjFace;
          Uint adjp = pnt[p]->adjPatch;
          Uint adji = pnt[p]->adjIndF;
          Point_T *adjpnt = 
            grid->patch[adjp]->interface[adjf]->point[adji];
          
          /* set flags and add to subface */
          adjpnt->touch = 1;
          adjpnt->copy  = 1;
          adjpnt->adjFace  = pnt[p]->face;
          adjpnt->adjPatch = pnt[p]->patch->pn;
          adjpnt->adjInd   = pnt[p]->ind;
          add_to_subface_scs(adjpnt);
        }
      }
    }
  }
  
}

/* pairing subfaces and setting BC Drichlet or Neumann BC consistently */
static void pair_subfaces_and_set_bc(Grid_T *const grid)
{
  /* struct for ranking of faces to more Nof the higher rank */
  struct face_rank_S
  {
    Uint pn;/* patch number */
    Uint fn;/* face number */
    Uint Nof;/* number of other faces */
  };
  
  const Uint Nfrank = grid->np*NFaces;
  /* face rank */
  struct face_rank_S *frank = calloc(Nfrank,sizeof(*frank));IsNull(frank);
  Uint nfr = 0;
  /* if this patch has Drichlet 1 otherwise 0 */
  Uint *isD = calloc(grid->np,sizeof(*isD));IsNull(isD);
  /* notation: o prefix letter means other */
  Interface_T **faces;
  Interface_T *face;
  Interface_T *oface;
  Patch_T *patch;
  Patch_T *opatch;
  SubFace_T *ssubf = 0;
  SubFace_T *osubf = 0;
  Uint bc;/* 0 = Dirichlet, 1 = Neumann */
  char sflgs[999] = {'\0'};
  Uint p,f,sf,i,j;
  
  /* pair all of the subfaces and count number of adjacent faces. */
  nfr = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    faces = grid->patch[p]->interface;
    
    FOR_ALL(f,faces)
    {
      frank[nfr].pn = grid->patch[p]->pn;
      frank[nfr].fn = f;
      frank[nfr].Nof = counter_n_adjacent_faces(faces[f]);
      nfr++;
      
      for (sf = 0; sf < faces[f]->ns; ++sf)
      {
        ssubf = faces[f]->subface[sf];
        
        if (!ssubf->touch)
          continue;
        
        osubf = find_subface(ssubf);/* juxtapose subface */
        ssubf->adjsn = osubf->sn;
      }
    }
  }
  assert(nfr == Nfrank);
  
  /* rank faces */
  for (i = 0; i < Nfrank; ++i)
  {
    for (j = i+1; j < Nfrank; ++j)
    {
      /* swap */
      if (frank[j].Nof > frank[i].Nof)
      {
        struct face_rank_S temp;
        temp     = frank[j];
        frank[j] = frank[i];
        frank[i] = temp;
      }
    }
  }
  
  /* preconditioning BC to balance Neumann and Dirichlet.
  // this is important to get converge elliptic solve. */
  if (grid->kind == Grid_SplitCubedSpherical_BHNS)
  FOR_ALL_PATCHES(p,grid)
  {
    patch = grid->patch[p];
    
    /* if condition number is high or elliptic solver fails
    // you can remove comments and try them. 
    // the ideal case that found out in experiments is equality 
    // in the number of Neuman and Drichlet bc for each "patch";
    // thus in future if elliptic solver fails you might try this. */
    //if (!IsItCovering(patch,"NS_around"))
    //if (!IsItCovering(patch,"BH_around"))
    if (IsItCovering(patch,"BH_around,NS_around"))
    {
      /* try to first set face Z = 1 to Dirichlet */
      i = 5;
    }
    else if (IsItCovering(patch,"NS"))
    {
      /* try to first set face Y = 1 to Dirichlet */
      i = 3;
    }
    else
    {
      continue;
    }
    
    /* set Dirichlet */
    {
      face  = patch->interface[i];
      int Favor_Dirichlet = 1;
      
      /* NOTE: an inner or outerB face is assumed to have ONE subface */
      if (face->outerB || face->innerB)
      {
        assert(face->ns == 1);
        continue;
      }
      if (face->df_dn_set)
        continue;
    
        
      /* test if we can favor Dirichlet: */
      
      /* if it has already Dirichlet */
      if (isD[patch->pn]) bc = 1;
      else                bc = 0;
      
      for (sf = 0; sf < face->ns; ++sf)
      {
        ssubf = face->subface[sf];
        if (!ssubf->touch) continue;
        
        opatch = grid->patch[ssubf->adjPatch];
        oface  = opatch->interface[ssubf->adjFace];
        
        /* if confilict => can't do Dirichlet */
        if(oface->df_dn_set && oface->df_dn != (bc+1)%2)
        {
          Favor_Dirichlet = 0;
          break;
        }
      }
      
      if (Favor_Dirichlet)
      {
        if (isD[patch->pn]) bc = 1;
        else                bc = 0;
      }
      else/* favor Nuemann */
      {
        if (isD[patch->pn]) bc = 0;
        else                bc = 1;
      }
      
      face->df_dn = bc%2;
      face->df_dn_set = 1;
      
      if (!face->df_dn) isD[patch->pn] = 1;
      
      for (sf = 0; sf < face->ns; ++sf)
      {
        ssubf = face->subface[sf];
        if (!ssubf->touch) continue;
        
        opatch = grid->patch[ssubf->adjPatch];
        oface  = opatch->interface[ssubf->adjFace];
        
        if(oface->df_dn_set && oface->df_dn != (bc+1)%2)
          Error0("Wrong BC arrangement!");
        
        oface->df_dn = (bc+1)%2;
        oface->df_dn_set = 1;
        
        if (!oface->df_dn) isD[opatch->pn] = 1;
      }
    }
  }
  
  
  /* set the other BC start with the largest number of neighbors.
  // start with Dirichlet since derivative might be not very accurate
  // in outermost which generally have the highest rank */
  for (i = 0; i < Nfrank; ++i)
  {
    patch = grid->patch[frank[i].pn];
    face  = patch->interface[frank[i].fn];
    int Favor_Dirichlet = 1;
    
    /* NOTE: an inner or outerB face is assumed to have ONE subface */
    if (face->outerB || face->innerB)
    {
      assert(face->ns == 1);
      continue;
    }
    if (face->df_dn_set)
      continue;
    
    /* test if we can favor Dirichlet: */
    
    /* if it has already Dirichlet */
    if (isD[patch->pn]) bc = 1;
    else                bc = 0;
    
    for (sf = 0; sf < face->ns; ++sf)
    {
      ssubf = face->subface[sf];
      if (!ssubf->touch) continue;
      
      opatch = grid->patch[ssubf->adjPatch];
      oface  = opatch->interface[ssubf->adjFace];
      
      /* if confilict => can't do Dirichlet */
      if(oface->df_dn_set && oface->df_dn != (bc+1)%2)
      {
        Favor_Dirichlet = 0;
        break;
      }
    }
    
    if (Favor_Dirichlet)
    {
      if (isD[patch->pn]) bc = 1;
      else                bc = 0;
    }
    else/* favor Nuemann */
    {
      if (isD[patch->pn]) bc = 0;
      else                bc = 1;
    }
    
    /* NOTE: innerB and outerB might considered as Neumann or Dirichlet. */
    face->df_dn = bc%2;
    face->df_dn_set = 1;
    
    if (!face->df_dn) isD[patch->pn] = 1;
    
    for (sf = 0; sf < face->ns; ++sf)
    {
      ssubf = face->subface[sf];
      if (!ssubf->touch) continue;
      
      opatch = grid->patch[ssubf->adjPatch];
      oface  = opatch->interface[ssubf->adjFace];
      
      if(oface->df_dn_set && oface->df_dn != (bc+1)%2)
        Error0("Wrong BC arrangement!");
      
      oface->df_dn = (bc+1)%2;
      oface->df_dn_set = 1;
      
      if (!oface->df_dn) isD[opatch->pn] = 1;
    }
  }
  
  /* set the subface flags correspondingly */
  FOR_ALL_PATCHES(p,grid)
  {
    faces = grid->patch[p]->interface;
    
    FOR_ALL(f,faces)
    {
      for (sf = 0; sf < faces[f]->ns; ++sf)
      {
        ssubf = faces[f]->subface[sf];
        
        if (!ssubf->touch)
          continue;
        
        assert(faces[f]->df_dn_set);
        ssubf->df_dn = faces[f]->df_dn;
        sprintf(sflgs,"%s,Dn:%u",ssubf->flags_str,ssubf->df_dn);
        Free(ssubf->flags_str);
        ssubf->flags_str = dup_s(sflgs);
      }
    }
  }
  
  /* free */
  Free(frank);
  Free(isD);
}

/* count the number of adjacent faces the give face has */
static Uint counter_n_adjacent_faces(const Interface_T *const face)
{
  Uint count = 0;
  Uint *pn = 0,*fn = 0;/* adjacent patch and face number */
  Uint sf;
  
  for (sf = 0; sf < face->ns; ++sf)
  {
    const SubFace_T *subf = face->subface[sf];
    
    if (!subf->touch)
      continue;
    
    /* check repetition */
    Flag_T rep = NO;
    Uint r;
    for (r = 0; r < count; ++r)
    {
      if (pn[r] == subf->adjPatch && fn[r] == subf->adjFace)
      {
        rep = YES;
        break;
      }
    }
    if (rep == YES)
      continue;
    
    /* no repetition => add them */
    pn = realloc(pn,(count+1)*sizeof(*pn)); 
    IsNull(pn);
    fn = realloc(fn,(count+1)*sizeof(*fn)); 
    IsNull(fn);
    
    pn[count] = subf->adjPatch;
    fn[count] = subf->adjFace;
    
    count++;
  }
  
  Free(pn);
  Free(fn); 
  
  return count;   
}

/* carry over the interfaces from old grid to the new one.
// NOTE: it empties old interface, i.e. puts its interface pointer to 0. */
void carryover_interfaces(Grid_T *const new_grid,Grid_T *const old_grid)
{
  FUNC_TIC
  
  /* some checks */
  assert(new_grid && old_grid);
  assert(new_grid->kind == Grid_SplitCubedSpherical_BHNS ||
         new_grid->kind == Grid_SplitCubedSpherical_NSNS ||
         new_grid->kind == Grid_SplitCubedSpherical_BHBH ||
         new_grid->kind == Grid_SplitCubedSpherical_SNS  ||
         new_grid->kind == Grid_SplitCubedSpherical_SBH);
    
  Uint p2,p1;
  Uint f,sf;
  
  FOR_ALL_PATCHES(p2,new_grid)
  {
    Patch_T *patch1   = 0;/* old patch */
    Patch_T *patch2   = new_grid->patch[p2];
    const char *name2 = strstr(patch2->name,"_");
    
    /* find the corresponding patch */
    FOR_ALL_PATCHES(p1,old_grid)
    {
      patch1             = old_grid->patch[p1];
      const char *name1  = strstr(patch1->name,"_");
  
      if (!strcmp(name2,name1))
        break;
    }/* FOR_ALL_PATCHES(p2,old_grid) */
    assert(patch1);
    
    /* move geometry */
    patch2->interface   = patch1->interface;
    patch2->innerB      = patch1->innerB;
    patch2->outerB      = patch1->outerB;
    patch2->is_a_closed = patch1->is_a_closed;
    patch2->is_b_closed = patch1->is_b_closed;
    patch2->is_c_closed = patch1->is_c_closed;
    patch1->interface   = 0;
                          
    /* update the internal pointers */
    Interface_T **face = patch2->interface;
    
    /* for all interfaces */
    FOR_ALL(f,face)
    {
      face[f]->patch = patch2;
      
      /* why is it not empty??? */
      assert(!face[f]->point);
      assert(!face[f]->innerP);
      assert(!face[f]->edgeP);
      
      /* for all subfaces */
      for (sf = 0; sf < face[f]->ns; ++sf)
      {
        SubFace_T *subf = face[f]->subface[sf];
        subf->patch = patch2;
      }
    }

  }/* end of FOR_ALL_PATCHES(p2,new_grid) */
  
  FUNC_TOC
}
