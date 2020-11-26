#ifndef physics_observe_LIB_H
#define physics_observe_LIB_H

#define OBSERVE_STR_LEN (99)

/* forward declaration */
struct PHYSICS_T;
struct GRID_T;


/* struct for physics observe */
typedef struct OBSERVE_T
{
  char quantity[OBSERVE_STR_LEN];/* which quantity is computed */
  double *ret;/* return value */
  struct GRID_T *grid;/* computational grid */
  struct PHYSICS_T *phys;/* physics pertinent to compact object of interest */
  void *items;/* this is general struct that composes 
               // the needed material and items to calculate 
               // the quantities of interest. this struct is populated
               // by plan and freed by free fucntions */
  unsigned Nitems;/* number of items */
  
}Observe_T;


int observe(struct PHYSICS_T *const phys,const char *const sq,double *const ret);

#undef OBSERVE_STR_LEN

#endif


