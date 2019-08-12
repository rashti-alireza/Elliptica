/*
// Alireza Rashti
// August 2018
*/

#include "solving_managings.h"
/* DeLimits */
#define DL_OC '{'
#define DL_CC '}'
#define COMMA ','
#define ARROW '>'

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

/* initializing Solving_Man_T struct and its elements */
void initialize_solving_man(Grid_T *const grid,sEquation_T **const field_eq,sEquation_T **const bc_eq,sEquation_T **const jacobian_field_eq,sEquation_T **const jacobian_bc_eq)
{
  const char *par_f = GetParameterS_E("Solving_Order");
  char *par;
  char **field_name = 0;
  unsigned nf = 0;
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
  
  /* filling solving_man in each patch */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    /* allocations */
    patch->solving_man = calloc(1,sizeof(*patch->solving_man));
    pointerEr(patch->solving_man);
    patch->solving_man->field_name = 
        calloc(nf,sizeof(*patch->solving_man->field_name));
    pointerEr(patch->solving_man->field_name);
    patch->solving_man->field_eq = calloc(nf,sizeof(*patch->solving_man->field_eq));
    pointerEr(patch->solving_man->field_eq);
    patch->solving_man->bc_eq = calloc(nf,sizeof(*patch->solving_man->bc_eq));
    pointerEr(patch->solving_man->bc_eq);
    patch->solving_man->jacobian_field_eq = calloc(nf,sizeof(*patch->solving_man->jacobian_field_eq));
    pointerEr(patch->solving_man->jacobian_field_eq);
    patch->solving_man->jacobian_bc_eq = calloc(nf,sizeof(*patch->solving_man->jacobian_bc_eq));
    pointerEr(patch->solving_man->jacobian_bc_eq);
    
    for (i = 0; i < nf; ++i)
    {
      patch->solving_man->field_name[i]  = dup_s(field_name[i]);
      patch->solving_man->field_eq[i]    = get_field_eq(field_name[i],field_eq);
      patch->solving_man->bc_eq[i]       = get_field_eq(field_name[i],bc_eq);
      patch->solving_man->jacobian_field_eq[i] = get_field_eq(field_name[i],jacobian_field_eq);
      patch->solving_man->jacobian_bc_eq[i]    = get_field_eq(field_name[i],jacobian_bc_eq);
    }
      
    patch->solving_man->nf    = nf;
    patch->solving_man->patch = patch;
    patch->solving_man->Frms  = DBL_MAX;
  }
  
  free_2d_mem(field_name,nf);
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
    abortEr_s("No such equation for \"%s\" exists.\n"
    "Note: if it is a Jacobian E.Q. and supposed to be for couple of fields, \n"
    "the naming convenstion is jacobian_f1_f2_... \n"
    "i.e. every related field must be mentioned in the name of function.\n",name);
  
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
