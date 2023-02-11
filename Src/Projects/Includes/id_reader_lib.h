#ifndef ID_READER_LIB_H
#define ID_READER_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct GRID_T;
struct ELLIPTICA_ID_READER_T;

/* strcut for initial data exporting  */
typedef struct ID_EXPORT_T
{
  struct GRID_T *grid;
  double *x,*y,*z;/* (x,y,z) coords */
  double *X,*Y,*Z;/* (X,Y,Z) coords */
  Uint *patchn;/* patch number for each coord */
  Uint npoints;/* number of coords */
  int **f_index;/* field index for each patch and for each field
                // ex: f_index[p][f] = for patch p and field f. */
}ID_Export_T;

void 
  idexp_load_Cartesian_coordinates_from_file
    (const char *const coords_file_path,ID_Export_T *const pnt);

void *
  idexp_new_binary_file_to_write
    (const char *const file_path,const char *const fields_name);


void idexp_close_file(FILE *file);

void
  idexp_interpolate_fields_and_write_to_file
    (FILE *const file,ID_Export_T *const pnt,
     const char *const fields_name_str/* comma separated */,
     const char *const evo_fields_name_str/* comma separated */);


void idexp_find_XYZ_from_xyz(struct ELLIPTICA_ID_READER_T *const idr, 
                             ID_Export_T *const pnt,const double *const CM);


void 
  idexp_interpolate_fields_and_save_in_array
    (struct ELLIPTICA_ID_READER_T *const idr, ID_Export_T *const pnt,
     const char *const fields_name_str/* comma separated */,
     const char *const evo_fields_name_str/* comma separated */);

void idexp_free(ID_Export_T *pnt);
ID_Export_T *idexp_init(void);

#endif
