#include "idr_header.h"
#include "fields_lib.h"
#include "maths_spectral_methods_lib.h"

#define STR_LEN_MAX (999)
/* NOTE: HEADER and FOOTER and END_MSG must be consistent
// while a file is read. */
#define HEADER "#{data#"
#define FOOTER "#}data#"
#define END_MSG "\n#file_completed#\n"

void 
  idr_load_Cartesian_coordinates_from_file
    (const char *const coords_file_path,ID_Reader_T *const pnt);
    
void *
  idr_new_binary_file_to_write
    (const char *const file_path,const char *const fields_name);


void idr_close_file(FILE *file);

void
  idr_interpolate_fields_and_write_to_file
    (FILE *const file,ID_Reader_T *const pnt,
     const char *const fields_name_str/* comma separated */,
     const char *const evo_fields_name_str/* comma separated */);


void idr_free(ID_Reader_T *pnt);
ID_Reader_T *idr_init(void);

void idr_find_XYZ_from_xyz(Elliptica_ID_Reader_T *const idr, 
                             ID_Reader_T *const pnt,const double *const CM);


void 
  idr_interpolate_fields_and_save_in_array
    (Elliptica_ID_Reader_T *const idr, ID_Reader_T *const pnt,
     const char *const fields_name_str/* comma separated */,
     const char *const evo_fields_name_str/* comma separated */);
 
 
