/*
// Alireza Rashti
// August 2018
*/

#include "solving_managing.h"

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

