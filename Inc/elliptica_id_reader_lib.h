#ifndef ELLIPTICA_ID_READER_LIB_H
#define ELLIPTICA_ID_READER_LIB_H

#define EIDR_MAX_STR (400)

// saves options and info to interact and control between the evo code and ID code
typedef struct ELLIPTICA_ID_READER_T
{
  char checkpoint_path[EIDR_MAX_STR];// path/to/elliptica/checkpoint/file
  const char *ifields;// input fields, e.g., "alpha,betax,adm_gxx"
  char *system;// the ID system, e.g., BHNS, NSNS, etc.
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

#undef EIDR_MAX_STR

#endif
