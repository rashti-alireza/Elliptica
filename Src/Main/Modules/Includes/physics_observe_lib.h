#ifndef physics_observe_LIB_H
#define physics_observe_LIB_H

#define OBSERVE_STR_LEN (99)

/* forward declaration */
struct PHYSICS_T;

/* struct for physics observe */
typedef struct OBSERVE_T
{
  char quantity[OBSERVE_STR_LEN];/* which quantity is computed */
  Grid_T *grid;/* computational grid */
  struct PHYSICS_T *phys;/* physics pertinent to compact object of interest */
  void *items;/* this is general struct that composes 
               // the needed material and items to calculate 
               // the quantities of interest. this struct is populated
               // by plan and freed by free fucntions */
  unsigned Nitems;/* number of items */
  
  /* functions to calculate momentums in each direction */
  double (*Px)(struct OBSERVE_T *const obs);
  double (*Py)(struct OBSERVE_T *const obs);
  double (*Pz)(struct OBSERVE_T *const obs);
  double (*Jx)(struct OBSERVE_T *const obs);
  double (*Jy)(struct OBSERVE_T *const obs);
  double (*Jz)(struct OBSERVE_T *const obs);
  
  /* function to calculate mass */
  double (*M)(struct OBSERVE_T *const obs);
}Observe_T;


int observe(struct PHYSICS_T *const phys,const char *const sq,double *const save);

#undef OBSERVE_STR_LEN

#endif


