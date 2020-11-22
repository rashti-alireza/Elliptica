#ifndef physics_observables_LIB_H
#define physics_observables_LIB_H

#define OBSERVABLE_STR_LEN (99)

/* forward declaration */
struct OBJ_MAN_T;

/* struct for physics observables */
typedef struct OBSERVABLE_T
{
  char quantity[OBSERVABLE_STR_LEN];/* which quantity is computed */
  Grid_T *grid;/* computational grid */
  struct OBJ_MAN_T *obj;/* compact object */
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
}Observable_T;


Observable_T *init_observable(struct OBJ_MAN_T *const obj,const char *const sq);
void free_observable(Observable_T *obs);

#undef OBSERVABLE_STR_LEN

#endif


