/*
// Alireza Rashti
// January 2019
*/

#include "solve_eqs_ddm_schur_complement.h"

/* 
   |B E||x|   |f|
   |F C||y| = |g|
   
   => Bx+Ey=f
      Fx+Cy=g
      
   => x = B^-1*(f-Ey)
   => (C-F*B^-1*E)y = g-F*B^-1*f
   
   Algorithm:
   
   1. solve BE'= E and Bf'= f
   2. compute g' = g - Ff'
   3. compute S  = C - FE'
   4. solve Sy = g'
   5. compute x = f'-E'y
  
// some definitions used here:
// subdomain = inner nodes + outer boundary points
// interface = all nodes - subdomains.
// note: I interface struct includes all boundary points.
*/

/* using Schur Complement domain decomposition method
// to solve equation. This method is capable of using direct solver
// like UMFPACK and also it is parallelizable.
*/
int ddm_schur_complement(Grid_T *const grid)
{
  char **field_name = 0;
  unsigned nf = 0;
  unsigned f;
  
  /* picking up labeling, mapping etc. */
  preparing_ingredients(grid);
  
  /* read order of fields to be solved from input */
  field_name = read_fields_in_order(&nf);
  
  /* solving fields in order */
  for (f = 0; f < nf; ++f)
  {
    printf("Solving Equation for %s:\n",field_name[f]);
    set_cf(grid,field_name[f]);/* solving_man->cf */
    solve_field(grid);/* solve field[f] */
  }
  
  /* free names */
  free_2d_mem(field_name,nf);
  
  return EXIT_SUCCESS;
}  

static int solve_field(Grid_T *const grid)
{
  /* residual determined in the input file */
  const double res_input = fabs(GetParameterD_E("Solving_Residual"));
  const int NumIter = GetParameterI_E("Linear_Solver_Number_of_Iteration");
  Flag_T IsItSolved = NO;
  int iter = 0;
  
  while (IsItSolved == NO && iter < NumIter)
  {
    unsigned p;
    
    DDM_SCHUR_COMPLEMENT_OpenMP(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *patch = grid->patch[p];
      make_f(patch);
      make_g_partial(patch);
    }
    //make_g(grid);
    
    IsItSolved = check_residual(grid,res_input);
    if (IsItSolved == YES)
      break;
   /*   
    DDM_SCHUR_COMPLEMENT_OpenMP(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *patch = grid->patch[p];
      making_B_and_E(patch);
      making_E_prime_and_f_prime(patch);
      making_F_and_C(patch);
      compute_g_prime(patch);
      compute_S(patch);
    }
    */
    /* solve Sy = g' */
    //solve_Sy_g_prime(grid);
    /*
    DDM_SCHUR_COMPLEMENT_OpenMP(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *patch = grid->patch[p];*/
      /* x = f'-E'y */
      //compute_x(patch);
      
    //}
    
  }/* end of while (IsItSolved == NO && iter < NumIter) */
  
  return EXIT_SUCCESS;
}

/* since Shur Complement needs:
// 1. specific labeling for grid 
// 2. boundary points counted only once
// 3. adjacency patch of each patch and how other 
// patches use this patch to impose their boundary conditions.
// this function provides this ingredients
*/
static void preparing_ingredients(Grid_T *const grid)
{
  unsigned p;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    /* if this path already has Schur stucture */
    if (patch->solving_man->method->Schur_Complement)
      continue;
      
    DDM_Schur_Complement_T *SchurC = calloc(1,sizeof(*SchurC));
    pointerEr(SchurC);
    SchurC->map = malloc(patch->nn*sizeof(*SchurC->map));
    pointerEr(SchurC->map);
    SchurC->inv = malloc(patch->nn*sizeof(*SchurC->inv));
    pointerEr(SchurC->inv);  
    patch->solving_man->method->Schur_Complement = 1;
    patch->solving_man->method->SchurC =  SchurC;
    
    /* making map and inv */
    make_map_and_inv(patch);
  }
  
}

/* making map and inv map for Schur complement method.
// inner mesh nodes come first, then outer boundary nodes 
// and then finally inner boundar nodes.
*/
static void make_map_and_inv(Patch_T *const patch)
{
  DDM_Schur_Complement_T *const SchurC = 
    patch->solving_man->method->SchurC;
  unsigned *const map = SchurC->map;
  unsigned *const inv = SchurC->inv;
  unsigned *Imap = 0;
  unsigned *Iinv = 0;
  const unsigned nn = patch->nn;
  const unsigned *n = patch->n;
  const unsigned nintfc = countf(patch->interface);
  unsigned i,j,j2, intfc;/* dummy indices and counters */
  
  /* keep tracking of points, 1 means counted, 0 means not */
  unsigned *flag_point = calloc(nn,sizeof(*flag_point));
  pointerEr(flag_point);
  
  /* filling inner points */
  j = 0;
  for (i = 0; i < nn; ++i)
  {
    if (!IsOnEdge(n,i))
    {
      map[i] = j;
      inv[j] = i;
      flag_point[i] = 1;
      j++;
    }
  }
  
  /* filling outer boundary points */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    unsigned nsfc = interface->ns;
    unsigned sfc;
    
    /* loop over all subfaces to filling outerbound */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      const unsigned *const B = subface->id;
      
      if (subface->outerB)/* if it reaches outer boundary */
      {
        for (i = 0; i < subface->np; ++i)
          if (!flag_point[B[i]])
          {
            map[B[i]] = j;
            inv[j] = B[i];
            flag_point[B[i]] = 1;
            j++;
          }
      }
    }
  }
  
  /* all of subdomain points already considered so : */
  patch->solving_man->method->SchurC->NS = j;
  j2 = 0;
  Imap = calloc(nn,sizeof(*Imap));
  pointerEr(Imap);
  Iinv = calloc(nn-j,sizeof(*Iinv));
  pointerEr(Iinv);
  
  /* filling the other remaining points */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    unsigned nsfc = interface->ns;
    unsigned sfc;
    
    /* loop over all subfaces to filling outerbound */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      const unsigned *const T = subface->id;
      
      if (!subface->innerB && subface->exterF)/* not boundary and internal faces */
      {
        for (i = 0; i < subface->np; ++i)
          if (!flag_point[T[i]])
          {
            map[T[i]] = j;
            inv[j] = T[i];
            Imap[T[i]] = j2;
            Iinv[j2] = T[i];
            flag_point[T[i]] = 1;
            ++j;
            ++j2;
          }
      }
    }
  }
  
  if (j != nn)
  {
    abortEr("Not all points are mapped.\n");
  }
  
  SchurC->Imap = Imap;
  SchurC->Iinv = Iinv;
  
  free(flag_point);
}

/* read order of fields to be solved from inout file.
// ->return value: name and number of fields.
*/
static char **read_fields_in_order(unsigned *const nf)
{
  char COMMA = ',';
  const char *par_f = GetParameterS_E("Solving_Order");
  char *par;
  char **field_name = 0;
  char *tok,*save = 0;
  
  /* finding fields's name */
  par = dup_s(par_f);/* par = f1,f2,... */
  tok = tok_s(par,COMMA,&save);/* tok = f1 */
  while(tok)
  {
    field_name = realloc(field_name,(*nf+1)*sizeof(*field_name));
    pointerEr(field_name);
    field_name[*nf] = dup_s(tok);
    tok = tok_s(0,COMMA,&save);/* tok = f2 */
    (*nf)++;
  }
  free(par);
  
  return field_name;
}

static void set_cf(Grid_T *const grid,const char *const field_name)
{
  unsigned p;
    
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    char **fields = patch->solving_man->field_name;
    unsigned nf = patch->solving_man->nf;
    patch->solving_man->cf = find_index_string(fields,nf,field_name);
  }
}

/* making f column.
// THEARD SAFA.
*/
static void make_f(Patch_T *const patch)
{
  DDM_Schur_Complement_T *Schur = patch->solving_man->method->SchurC;
  assert(!Schur->f);
  double *f = alloc_double(patch->nn);
  Schur->f = f;
  
  f_in_equation_part(patch);
  f_in_outerboundary_part(patch);
  
}

/* making partial of g which coming from patch.
// for example if g = y1-y2, g partly coming from patch1
// and partly coming from patch2. when all of partial g's are made
// the g coloumn itself will be populated.
*/
static void make_g_partial(Patch_T *const patch)
{
  DDM_Schur_Complement_T *Schur = patch->solving_man->method->SchurC;
  assert(!Schur->pg);
  double *pg = alloc_double(patch->nn-Schur->NS);
  Schur->pg = pg;
  const unsigned nintfc = countf(patch->interface);
  unsigned intfc;
  
  /* loop over all interfaces */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    unsigned nsfc = interface->ns;
    unsigned sfc;
    
    /* loop over all subfaces */
    for (sfc = 0; sfc < nsfc; ++sfc)
{
      SubFace_T *subface = interface->subface[sfc];
      
      if (!subface->exterF)/* if subface is internal */
      {
        abortEr(INCOMPLETE_FUNC);
      }
      else if (subface->innerB)/* if there is inner boundary */
      {
        abortEr(INCOMPLETE_FUNC);
      }
      else if (subface->outerB)/* if it reaches outer boundary */
      {
        continue;
      }
      else if (subface->touch)/* if two patches are in touch */
      {
        if (subface->copy)/* if the collocated point */
        {
          /* copy values */
          b_bndry_copy_ppm(&BC);
        }
        else
        {
          /* interpolate values */
          b_bndry_interpolate_ppm(&BC);
        }
      }
      else /* if there is an overlap case */
      {
        /* interpolate values */
        b_bndry_interpolate_ppm(&BC);
      }
      
    }/* end of for (sfc = 0; sfc < nsfc; ++sfc) */
    
  }/* end of for (intfc = 0; intfc < nintfc; ++intfc) */

}

/* calculating the part of f coming from equation */
static void f_in_equation_part(Patch_T *const patch)
{
  Solving_Man_T *const S      = patch->solving_man;
  const unsigned cf           = S->cf;
  const char *const fname     = S->field_name[cf];
  fEquation_T *const field_eq = S->field_eq[cf];
  Field_T *const field        = patch->pool[Ind(fname)];
  
  field_eq(field,S->method->SchurC);
}

/* calculating the part of f coming from outerboundary points */
static void f_in_outerboundary_part(Patch_T *const patch)
{
  const unsigned nintfc    = countf(patch->interface);
  Solving_Man_T *const S   = patch->solving_man;
  const unsigned cf        = S->cf;
  const char *const fname  = S->field_name[cf];
  fEquation_T *const bc_eq = S->bc_eq[cf];
  Field_T *const field     = patch->pool[Ind(fname)];
  unsigned intfc;
  
  /* loop over all interfaces */
  for (intfc = 0; intfc < nintfc; ++intfc)
  {
    Interface_T *interface = patch->interface[intfc];
    unsigned nsfc = interface->ns;
    unsigned sfc;
    
    /* loop over all subfaces and look for outerboundary */
    for (sfc = 0; sfc < nsfc; ++sfc)
    {
      SubFace_T *subface = interface->subface[sfc];
      Boundary_Condition_T bc;
      bc.subface = subface;
      bc.node    = subface->id;
      bc.nn      = subface->np;
      bc.field   = field;
      bc.patch   = patch;
      
      if (subface->outerB)/* if it reaches outer boundary */
        bc_eq(&bc,S->method->SchurC);
    }
  }
}



/* find out the residual of each patch and decide weather the equations
// are already solved or not.
// ->return value: YES if EQs are solved, NO otherwise.
*/
static Flag_T check_residual(const Grid_T *const grid,const double res_input)
{
  Flag_T flg = YES;
  
  return flg;
  UNUSED(grid);
  UNUSED(res_input);
}
