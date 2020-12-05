#include "core_lib.h"
#include "error_handling_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "fields_lib.h"
#include "maths_spectral_methods_lib.h"


#define STR_LEN_MAX (999)
/* NOTE: HEADER and FOOTER and END_MSG must be consistent
// while a file is read. */
#define HEADER "#{data#"
#define FOOTER "#}data#"
#define END_MSG "\n#file_completed#\n"

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

