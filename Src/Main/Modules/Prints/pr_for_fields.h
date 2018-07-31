#include <silo.h>
#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "text_tools_lib.h"
#include "macros_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"

#define MAX_STR_LEN 500
#define FORMAT_ER_PAR "Incorrect parameter format.\n"\
    "Print_Fields_4D must be written like:\n"\
      "Print_Fields_4D = yes: format:..,{(field1,field2,...)vs(x,y,z)|coord}...\n"

struct Info_S
{
  char (*field)[MAX_STR_LEN];
  char axis[3][MAX_STR_LEN];
  char coord[MAX_STR_LEN];
  unsigned nf;/* number of fields */
};

void pr_fields(Pr_Field_T *const pr);
Pr_Field_T *init_PrField(const Grid_T *const grid);
void free_PrField(Pr_Field_T *pr);
static void read_parameter_4d(const char *const par,Pr_Field_T *const pr);
static void free_info_s(Pr_Field_T *const pr);
static void pr_fields_on_grid_HDF5_4d(Pr_Field_T *const pr);
static void pr_structured_mesh_3d_silo(const Pr_Field_T *const pr);
static void prepare_node_structured_mesh_3d_silo(const char *const type,const Patch_T *const patch,float **const x,float **const y,float **const z);
static void free_nodes_silo(float *x,float *y,float *z);
