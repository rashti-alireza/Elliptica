#ifndef physics_observables_lib
#define physics_observables_lib


/* struct for physics observables */
typedef struct OBSERVABLE_T
{
  const char *quantity;/* which quantity is computed */
  void *grid;/* computational grid */
  /* ADM momentums */
  struct
  {
    void *patch;/* the patch in which the following variables are defined */
    /* physical metric components */
    double *g00;
    double *g01;
    double *g02;
    double *g11;
    double *g12;
    double *g22;
    /* normal vector at the surface S, outward */
    double *n_U0;
    double *n_U1;
    double *n_U2;
    /* integration flags */
    unsigned surface_integration_flg: 1;/* if 1 means it measn 
                                        // we need surface integration 
                                        // on this patch as well, 
                                        // 0 means, no need */
    /* which hypersurface the surface integral is carried out */
    unsigned X_surface: 1;
    unsigned Y_surface: 1;
    unsigned Z_surface: 1;
    /* index of hypersurface for each X,Y and Z respectively */
    unsigned I;
    unsigned J;
    unsigned K;
    
  }**ADM;
  unsigned N_ADM;/* number of ADM struct */
  /* functions to calculate ADM momentums in each direction */
  
  double (*Px_ADM)(struct OBSERVABLE_T *const obs);
  double (*Py_ADM)(struct OBSERVABLE_T *const obs);
  double (*Pz_ADM)(struct OBSERVABLE_T *const obs);
  double (*Jx_ADM)(struct OBSERVABLE_T *const obs);
  double (*Jy_ADM)(struct OBSERVABLE_T *const obs);
  double (*Jz_ADM)(struct OBSERVABLE_T *const obs);
  
  /* function to calculate masses */
  double (*ADM_mass)(struct OBSERVABLE_T *obs);
  double (*Komar_mass)(struct OBSERVABLE_T *obs);
}Observable_T;

Observable_T *init_observable(void *grid);
void plan_observable(Observable_T *const obs);
void free_observable(Observable_T *obs);


#endif


