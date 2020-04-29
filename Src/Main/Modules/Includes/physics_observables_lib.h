#ifndef physics_observables_LIB_H
#define physics_observables_LIB_H

/* struct for physics observables */
typedef struct OBSERVABLE_T
{
  const char *quantity;/* which quantity is computed */
  Grid_T *grid;/* computational grid */
  void *items;/* this is general struct that composes 
               // the needed material and items to calculate 
               // the quantities of interest. this struct is populated
               // by plan and freed by free fucntions */
  unsigned Nitems;/* number of items */
  
  /* functions to calculate momentums in each direction */
  double (*Px)(struct OBSERVABLE_T *const obs);
  double (*Py)(struct OBSERVABLE_T *const obs);
  double (*Pz)(struct OBSERVABLE_T *const obs);
  double (*Jx)(struct OBSERVABLE_T *const obs);
  double (*Jy)(struct OBSERVABLE_T *const obs);
  double (*Jz)(struct OBSERVABLE_T *const obs);
  
  /* function to calculate mass */
  double (*M)(struct OBSERVABLE_T *const obs);
  
  /* preparing functions and freeing functions of items stuct
  // these functions are populated during initilization of Observable_T */
  void (*plan_items)(struct OBSERVABLE_T *obs);
  void (*free_items)(struct OBSERVABLE_T *obs);
}Observable_T;


Observable_T *init_observable(void *grid,void (*plan_items)(struct OBSERVABLE_T *obs),void (*free_items)(struct OBSERVABLE_T *obs));
void plan_observable(Observable_T *const obs);
void free_observable(Observable_T *obs);


#endif


