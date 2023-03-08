#ifndef ELLIPTICA_ID_READER_LIB_H
#define ELLIPTICA_ID_READER_LIB_H

// Expose it for C++ interfacing
#ifdef __cplusplus
extern "C" {
#endif
/* ---------------------------------------------------------------------- */


/* NOTE: this is the header to be used for an evolution code */
/* NOTE: unsigned is not used here in case there is any conflicts for evo codes. */

// saves options and info to interact and control between the evo code and ID code
typedef struct ELLIPTICA_ID_READER_T
{
  char *checkpoint_path;// path/to/elliptica/checkpoint/file
  const char *ifields;// input fields, e.g., "alpha,betax,adm_gxx"
  char *system;// the ID system, e.g., BHNS, NSNS, etc.
  char *option;// options for exportation, e.g., asymptotically_inertial
  unsigned nfield;// total number of Field_Dictionary
  double **field;// value = field[indx(Field_Dictionary[i])][ijk]
  int npoints;// length of 1D coord array, use int(not Uint) as most evo codes use int
  double *x_coords;// Cartesian x coord values, x_coords[ijk] with length npoints
  double *y_coords;// Cartesian y coord values, y_coords[ijk] with length npoints
  double *z_coords;// Cartesian z coord values, z_coords[ijk] with length npoints
  unsigned (*indx)(const char *const fname);// find field index  
  void (*set_param)(const char *const lv,const char *const rv,
                    struct ELLIPTICA_ID_READER_T *const idr);// set params from evo code side
  double (*get_param_dbl)(const char *const lv,
                          struct ELLIPTICA_ID_READER_T *const idr);// get a double param from a checkpoint file
  char **params_lv;// left value,  params[0]="force_balance"
  char **params_rv;// right value, params[0]="on"
  unsigned nparams;// number of params
}Elliptica_ID_Reader_T;

Elliptica_ID_Reader_T *elliptica_id_reader_init (
  const char *const checkpnt/* path/to/elliptica/checkpoint/file */,
  const char *const option/* the option for export */
  );

int elliptica_id_reader_interpolate(Elliptica_ID_Reader_T *const idr);
int elliptica_id_reader_free(Elliptica_ID_Reader_T *idr);


/* ---------------------------------------------------------------------- */
#ifdef __cplusplus
}
#endif

#endif
