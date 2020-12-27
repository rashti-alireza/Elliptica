#ifndef physics_observe_LIB_H
#define physics_observe_LIB_H
#include "elliptica_system_lib.h"

#define OBSERVE_STR_LEN (99)

/* forward declaration */
struct PHYSICS_T;
struct GRID_T;


/* struct for physics observe */
typedef struct OBSERVE_T
{
  char quantity[OBSERVE_STR_LEN];/* which quantity is computed */
  const char *method;/* method used in calculation */
  double *ret;/* return value */
  struct GRID_T *grid;/* computational grid */
  struct PHYSICS_T *phys;/* physics pertinent to compact object of interest */
  void *items;/* this is general struct that composes 
               // the needed material and items to calculate 
               // the quantities of interest. this struct is populated
               // by plan and freed by free fucntions */
  Uint Nitems;/* number of items */
  
}Observe_T;

int observe_main(struct PHYSICS_T *const phys);
int observe(struct PHYSICS_T *const phys,const char *const sq,
            const char *const method,double *const ret);

#undef OBSERVE_STR_LEN

#endif


