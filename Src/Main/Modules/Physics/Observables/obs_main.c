/*
// Alireza Rashti
// September 2019
*/

/* synopsis:
// =========
//
// * initialize observable *
// Observable_T *obs = init_observable(grid,plan_items_func,free_items_func);
//
// * specifiy which obeservable *
// obs->quantity = "ADM_momentums" # means ADM momentums
//
// * plan observable *
// plan_observable(obs);# it finds out the related patches, physical metric etc.
//
// * calculate the observable *
// double Px_ADM = obs->Px_ADM(obs);# x component
// double Py_ADM = obs->Py_ADM(obs);# y component
// double Pz_ADM = obs->Pz_ADM(obs);# z component
// double Jx_ADM = obs->Jx_ADM(obs);# x component of angular momentum
// double Jy_ADM = obs->Jy_ADM(obs);# y component of angular momentum
// double Jz_ADM = obs->Jz_ADM(obs);# z component of angular momentum
//
// *free*
// free_observable(obs);
*/

#include "obs_main.h"

/* initialzing stuct Observable_T */
Observable_T *init_observable(void *grid,void (*plan_items)(struct OBSERVABLE_T *obs),void (*free_items)(struct OBSERVABLE_T *obs))
{
  Observable_T *const obs = calloc(1,sizeof(*obs));
  IsNull(obs);

  obs->grid = grid;
  obs->plan_items = plan_items;
  obs->free_items = free_items;
  
  return obs;
}

/* planning observable according to each grid and project */
void plan_observable(Observable_T *const obs)
{
  if (obs->plan_items)
    obs->plan_items(obs);
}

/* free stuct Observable_T and items */
void free_observable(Observable_T *obs)
{
  if (!obs)
    return;
  if (obs->free_items)
    obs->free_items(obs);
  free(obs);
}

