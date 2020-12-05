#include "core_lib.h"
#include "error_handling_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "fields_lib.h"
#include "maths_spectral_methods_lib.h"


#define STR_LEN_MAX (1000)
#define HEADER "#{data#"
#define FOOTER "#}data#"
#define END_MSG "\n#file_completed#\n"

void load_Cartesian_coordinates_from_file
                        (const char *const coords_file_path,
                         ID_Export_T *const pnt);

