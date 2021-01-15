/*
// Alireza Rashti
// August 2018
*/

#include "solve_manager.h"
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
  Uint ndb = 0;
  
  /* count number of data base */
  if (db != 0)
    for (ndb = 0; db[ndb] != 0; ++ndb);
  else
    ndb = 0;
  
  /* check if this is not already exists */
  for (Uint i = 0; i < ndb; ++i)
  {
    if (strcmp_i(db[i]->name,name))
      Errors("Equation '%s' has been already added!",name);
  }
  
  db = realloc(db,(ndb+2)*sizeof(*db));
  IsNull(db);
  db[ndb+1] = 0;
  db[ndb] = calloc(1,sizeof(*db[ndb]));
  IsNull(db[ndb]);
  
  sprintf(db[ndb]->name,name);
  db[ndb]->eq = eq;
  
  *data_base = db;
}

/* initializing Solving_Man_T struct and its elements */
void initialize_solving_man(Grid_T *const grid,
                            sEquation_T **const field_eq,
                            sEquation_T **const bc_eq,
                            sEquation_T **const jacobian_field_eq,
                            sEquation_T **const jacobian_bc_eq,
                            const char *const par_prefix/* prefix for eq param, ex: Eq_phi1 */)
{
  const char *par_f = Pgets("solve_order");
  char eq_fname[STR_LEN1]  = {'\0'};
  char par_fname[STR_LEN0] = {'\0'};
  char val_fname[STR_LEN0] = {'\0'};
  char **field_name = 0;
  char *tok  = 0;
  char *save = 0;
  char *par  = 0;
  char *aux  = 0;
  Uint nf = 0;
  Uint p,i;
  
  /* finding fields's name */
  par = dup_s(par_f);/* par = f1,f2,... */
  tok = tok_s(par,COMMA,&save);/* tok = f1 */
  while(tok)
  {
    field_name = realloc(field_name,(nf+1)*sizeof(*field_name));
    IsNull(field_name);
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
    if (!patch->solving_man)
    {
      patch->solving_man = calloc(1,sizeof(*patch->solving_man));
      IsNull(patch->solving_man);
    }
    
    if (patch->solving_man->field_name)
    {
      free_2d_mem(patch->solving_man->field_name,patch->solving_man->nf);
    }
    patch->solving_man->field_name = calloc(nf,sizeof(*patch->solving_man->field_name));
    IsNull(patch->solving_man->field_name);
    
    Free(patch->solving_man->field_eq);
    patch->solving_man->field_eq = calloc(nf,sizeof(*patch->solving_man->field_eq));
    IsNull(patch->solving_man->field_eq);
    
    Free(patch->solving_man->bc_eq);
    patch->solving_man->bc_eq = calloc(nf,sizeof(*patch->solving_man->bc_eq));
    IsNull(patch->solving_man->bc_eq);
    
    Free(patch->solving_man->jacobian_field_eq);
    patch->solving_man->jacobian_field_eq = calloc(nf,sizeof(*patch->solving_man->jacobian_field_eq));
    IsNull(patch->solving_man->jacobian_field_eq);
    
    Free(patch->solving_man->jacobian_bc_eq);
    patch->solving_man->jacobian_bc_eq = calloc(nf,sizeof(*patch->solving_man->jacobian_bc_eq));
    IsNull(patch->solving_man->jacobian_bc_eq);
    
    for (i = 0; i < nf; ++i)
    {
      /* get the prefix of equation name if par_prefix != 0
      // ex:
      // par_fname = "Eq_phi"
      // val_fname = "XCTS_curve,NS"
      // => eq_name = "XCTS_curve_phi". */
      if (par_prefix)
      {
        /* get the param name */
        sprintf(par_fname,"%s%s",par_prefix,field_name[i]);
        /* get param value */
        sprintf(val_fname,"%s",Pgets(par_fname));
        /* get rid of comma, since there must be a comma in param value */
        aux = strchr(val_fname,',');
        if (!aux)
        {
          Errors("Could not find comma in parameter '%s'!",par_fname);
        }
        aux[0] = '\0';
        aux = val_fname;
        sprintf(eq_fname,"%s_%s",aux,field_name[i]);
      }
      else
      {
       sprintf(eq_fname,"%s",field_name[i]); 
      }
      
      patch->solving_man->field_name[i] = dup_s(field_name[i]);
      patch->solving_man->field_eq[i]   = get_field_eq(eq_fname,field_eq,Prefix_EQ);
      patch->solving_man->bc_eq[i]      = get_field_eq(eq_fname,bc_eq,Prefix_BC);
      patch->solving_man->jacobian_field_eq[i] = get_field_eq(eq_fname,jacobian_field_eq,Prefix_Jac_EQ);
      patch->solving_man->jacobian_bc_eq[i]    = get_field_eq(eq_fname,jacobian_bc_eq,Prefix_Jac_BC);
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
static fEquation_T *get_field_eq(const char *const name, 
                                 sEquation_T **const db,
                                 const char *const prefix)
{
  fEquation_T *eq = 0;
  char eq_name[STR_LEN2] = {'\0'};
  Uint i;
  
  sprintf(eq_name,"%s%s",prefix,name);
  
  if (!db)
    Error0("Empty data base!");
  
  for (i = 0; db[i] != 0; ++i)
  {
    if (strcmp_i(db[i]->name,eq_name))
    {
      eq = db[i]->eq;
      break;
    }
  }
  
  if (!eq)
    Errors("No such equation \"%s\" found!",eq_name);
  
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
    Errors("There is no such solver defined.\n",solver);
    
  return solver_eq;
}

/* freeing Equation_T data base */
void free_db_eqs(sEquation_T **db)
{
  Uint ndb;
  
  if (!db) return;
  
  for (ndb = 0; db[ndb] != 0; ++ndb)
  {
    free(db[ndb]);
  }
  
  free(db);
}



