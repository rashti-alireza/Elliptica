#include "core_lib.h"
#include "maths_equation_solvings_lib.h"
#include "maths_general_lib.h"

Root_Finder_T *init_root_finder(const Uint n);
double *execute_root_finder(Root_Finder_T *const root);
void plan_root_finder(Root_Finder_T *const root);
void free_root_finder(Root_Finder_T *root);
static double *root_finder_steepest_descent(Root_Finder_T *const root);
static double g_SD(double (**f)(void *params,const double *const x),void *params,const double *const x,Root_Finder_T *const root);
static double dg_dx_FD3M_SD(void *params,double *const x,const Uint dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const Uint dir),Root_Finder_T *const root);
static double dg_dx_FD3L_SD(void *params,double *const x,const Uint dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const Uint dir),Root_Finder_T *const root);
static double dg_dx_FD3R_SD(void *params,double *const x,const Uint dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const Uint dir),Root_Finder_T *const root);
static double dg_dx_of_df_dx_SD(void *params,double *const x,const Uint dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const Uint dir),Root_Finder_T *const root);
void print_root_finder_exit_status(const Root_Finder_T *const root);
static double *root_finder_bisect_single(Root_Finder_T *const root);


