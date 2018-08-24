/*
// Alireza Rashti
// August 2018
*/

#include "solving_managings.h"
/* DeLimits */
#define DL_OC '{'
#define DL_CC '}'
#define COMMA ','
#define FLASH '>'

/* initializing Equation_T by null*/
void *init_eq(void)
{
  return 0;
}

/* add an equation to Equation_T data base.
// data_base is data base
// eq is the equation to be added to data base
// name is the name of equation that this equation is referred to.
// note: the end of data base is determined by null pointer
*/
void add_eq(sEquation_T ***const data_base, fEquation_T *const eq,const char *const name)
{
  sEquation_T **db = *data_base;
  unsigned ndb = 0;
  
  /* count number of data base */
  if (db != 0)
    for (ndb = 0; db[ndb] != 0; ++ndb);
  else
    ndb = 0;
  
  db = realloc(db,(ndb+2)*sizeof(*db));
  pointerEr(db);
  db[ndb+1] = 0;
  db[ndb] = calloc(1,sizeof(*db[ndb]));
  pointerEr(db[ndb]);
  
  sprintf(db[ndb]->name,name);
  db[ndb]->eq = eq;
  
  *data_base = db;
}

/* allocating and filling Solution_Man_T struct and its elements 
// in each patch according to input file.
*/
void populate_solution_man(Grid_T *const grid,sEquation_T **const field_eq,sEquation_T **const bc_eq)
{
  const char *par_f = GetParameterS_E("Fields");
  const char *par_o = GetParameterS_E("Solving_Order");
  char *par;
  char **field_name = 0;
  char **group = 0;
  unsigned nf = 0;
  unsigned ng = 0;
  char *tok,*save = 0;
  unsigned p,i;
  
  /* finding fields's name */
  par = dup_s(par_f);/* par = f1,f2,... */
  tok = tok_s(par,COMMA,&save);/* tok = f1 */
  while(tok)
  {
    field_name = realloc(field_name,(nf+1)*sizeof(*field_name));
    pointerEr(field_name);
    field_name[nf] = dup_s(tok);
    tok = tok_s(0,COMMA,&save);/* tok = f2 */
    nf++;
  }
  free(par);
  
  /* filling solution_man in each patch */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    patch->solution_man = calloc(1,sizeof(*patch->solution_man));
    pointerEr(patch->solution_man);
    patch->solution_man->field_name = 
        calloc(nf,sizeof(*patch->solution_man->field_name));
      pointerEr(patch->solution_man->field_name);
      
    for (i = 0; i < nf; ++i)
      patch->solution_man->field_name[i] = dup_s(field_name[i]);
      
    patch->solution_man->nf = nf;
  }
  free_2d_mem(field_name,nf);
  
  /* finding solving order */
  save = 0;
  par = dup_s(par_o);/* par = {f1,f2,...}->{f3,f4,...}->... */
  tok = tok_s(par,FLASH,&save);/* =>tok = {f1,f2,...}- */
  while(tok)
  {
    char *sub,*save_sub;
    char *tok2,*save2 = 0;
    if (!check_format_s(tok,"?{?}?"))
      abortEr(FORMAT_ER_PAR);
      
    ng = 0;
    group = 0;
    sub  = sub_s(tok,DL_OC,DL_CC,&save_sub);/* sub = f1,f2,... */
    tok2 = tok_s(sub,COMMA,&save2);/* tok2 = f1 */
    while(tok2)
    {
      group = realloc(group,(ng+1)*sizeof(*group));
      pointerEr(group);
      group[ng] = dup_s(tok2);
      tok2 = tok_s(0,COMMA,&save2);/* tok2 = f2 */
      ng++;
    }/* group = {f1,f2,...}*/
    
    /* filling solution structure */
    if (!ng)
      abortEr("No group fields to be solved are found!\n");
      
    fill_solve(grid,group,ng,field_eq,bc_eq);
    
    tok = tok_s(0,FLASH,&save);/* tok = {f3,f4,...} */
    free_2d_mem(group,ng);
  }
  free(par);
}

/* allocating memory and filling solve struct group by group */
static void fill_solve(Grid_T *const grid,char **const group,const unsigned ng,sEquation_T **const field_eq,sEquation_T **const bc_eq)
{
  unsigned p;
  
  assert(ng);
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Solve_T *solve;
    double *b = alloc_double(ng*patch->nn);
    unsigned ns,i;
    
    ns = ++patch->solution_man->ns;
    solve = alloc_solve(patch,ns);
    solve->b = b;
    solve->nf = ng;
    solve->field = calloc(ng,sizeof(*solve->field));
    pointerEr(solve->field);
    solve->f_occupy = calloc(ng,sizeof(*solve->f_occupy));
    pointerEr(solve->f_occupy);
    solve->field_eq = calloc(ng,sizeof(*solve->field_eq));
    pointerEr(solve->field_eq);
    solve->bc_eq = calloc(ng,sizeof(*solve->bc_eq));
    pointerEr(solve->bc_eq);
    solve->f_name = calloc(ng,sizeof(*solve->f_name));
    pointerEr(solve->f_name);
    
    for (i = 0; i < ng; ++i)
    {
      solve->f_name[i] = dup_s(group[i]);
      solve->f_occupy[i] = i*patch->nn;
      solve->field_eq[i] = get_field_eq(group[i],field_eq);
      solve->bc_eq[i]    = get_field_eq(group[i],bc_eq);
      solve->solver      = get_solver_method(GetParameterS_E("Linear_Solver"));
      solve->field[i]    = prepare_field(group[i],"(3dim)",patch);
    }
    
  }/* end of FOR_ALL_PATCHES(p,grid) */
}

/* given name of equation and data base of equations 
// it returns the corresponding equation to that given name.
// ->return value: found equation, error if not found
*/
fEquation_T *get_field_eq(const char *const name, sEquation_T **const db)
{
  fEquation_T *eq = 0;
  unsigned i;
  
  if (!db)
    return 0;
  
  for (i = 0; db[i] != 0; ++i)
  {
    if (strstr_i(db[i]->name,name))
    {
      eq = db[i]->eq;
      break;
    }
  }
  
  if (!eq)
    abortEr_s("No such equation for \"%s\" exists.\n",name);
  
  return eq;
}

/* given name of solver it returns its function for that solver
// ->return value: solver, error if not found
*/
fEquation_Solver_T *get_solver_method(const char *const solver)
{
  fEquation_Solver_T *solver_eq = 0;
  
  if (strcmp_i(solver,"UMFPACK"))
    solver_eq = direct_solver_umfpack_di;
  else if (strcmp_i(solver,"UMFPACK_long"))
    solver_eq = direct_solver_umfpack_dl;
  else
    abortEr_s("There is no such solver defined.\n",solver);
    
  return solver_eq;
}

/* given name, attribute and patch
// it returns a field with demanded properties.
// ->return value: new made field
*/
static Field_T *prepare_field(const char *const name,const char *const attr,Patch_T *const patch)
{
  Field_T *field = calloc(1,sizeof(*field));
  pointerEr(field);
  
  field->name = dup_s(name);
  field->attr = dup_s(attr);
  field->v = alloc_double(patch->nn);
  field->patch = patch;
  
  return field;
}
