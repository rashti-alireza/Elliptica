#ifndef ELLIPTICA_ID_READER_LIB_H
#define ELLIPTICA_ID_READER_LIB_H

#include "elliptica_system_lib.h"

#define EIDR_MAX_STR (400)

// saves options and info to interact and control between the evo code and ID code
typedef struct ELLIPTICA_ID_READER_T
{
  char checkpoint_path[EIDR_MAX_STR];// path/to/elliptica/checkpoint/file
  const char *ifields;// input fields, e.g., "alpha,betax,adm_gxx"
  char *system;// the ID system, e.g., BHNS, NSNS, etc.
  char *option;// options for exportation, e.g., asymptotically_inertial
  Uint nfield;// total number of Field_Dictionary
  double **field;// value = field[indx(Field_Dictionary[i])][ijk]
  double npoints;// number of (x,y,z) coords
  double *x_coords;// Cartesian x coord values, x_coords[ijk]
  double *y_coords;// Cartesian y coord values, y_coords[ijk]
  double *z_coords;// Cartesian z coord values, z_coords[ijk]
  Uint (*indx)(const char *const fname);// find field index  
  void (*param)(const char *const lv,const char *const rv,
                struct ELLIPTICA_ID_READER_T *const idr);// set params from evo code
  char **params_lv;// left value,  params[0]="force_balance"
  char **params_rv;// right value, params[0]="on"
  Uint nparams;// number of params
}Elliptica_ID_Reader_T;

Elliptica_ID_Reader_T *elliptica_id_reader_init (
  const char *const checkpnt/* path/to/elliptica/checkpoint/file */,
  const char *const option/* the option for export */
  );

int elliptica_id_reader_interpolate(Elliptica_ID_Reader_T *const idr);
int elliptica_id_reader_free(Elliptica_ID_Reader_T *idr);


#undef EIDR_MAX_STR

#endif