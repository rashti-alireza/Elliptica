#ifndef checkpoint_LIB_H
#define checkpoint_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct PHYSICS_T;
struct GRID_T;

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


void write_checkpoint(struct PHYSICS_T *const phys,const char *const out_dir);
int is_checkpoint_sound(const char *const file_path);
int can_we_use_checkpoint(const char *const cur_out_dir);

void *open_checkpoint_file_then_read_grid_and_params(struct PHYSICS_T *const phys);
void read_fields_from_checkpoint_file(struct PHYSICS_T *const phys,FILE *const file);

Parameter_T *parameter_query_from_checkpoint(const char *const par_name,FILE *const file);

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


void idexp_free(ID_Export_T *pnt);
ID_Export_T *idexp_init(void);

#endif




