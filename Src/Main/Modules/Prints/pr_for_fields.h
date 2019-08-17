#include <silo.h>
#include "core_lib.h"/* it includes prints_lib.h too */
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "macros_lib.h"
#include "maths_general_lib.h"
#include "memory_managing_lib.h"
#include "coordinates_lib.h"

#define MAX_STR_LEN 500
#define FORMAT_ER_PAR "Incorrect parameter format.\n"\
    "Print_Fields_4D must be written like:\n"\
      "Print_Fields_4D = yes: format:..,{field1,field2,...,(fx,fy,fz),...}\n"

struct Info_S
{
  char *field;
  char *comp[3];/* components name */
  unsigned vec_flg:1;/* 1 if it is vector, 0 otherwise */
};

void pr_fields(Pr_Field_T *const pr);
Pr_Field_T *init_PrField(const Grid_T *const grid);
void free_PrField(Pr_Field_T *pr);
static void read_parameter_4d(const char *const par,Pr_Field_T *const pr);
static void free_info_s(Pr_Field_T *const pr);
static void pr_fields_on_grid_HDF5_4d(Pr_Field_T *const pr);
static void pr_structured_mesh_3d_silo(const Pr_Field_T *const pr);
static void pr_scalar_on_structured_mesh_3d_silo(const Pr_Field_T *const pr);
static void pr_vector_on_structured_mesh_3d_silo(const Pr_Field_T *const pr);
static void *make_structured_mesh_3d_Cartesian(Pr_Field_T *const pr,const Patch_T *const patch);
static void *make_structured_mesh_3d_Curvilinear(Pr_Field_T *const pr,const Patch_T *const patch);
static void prepare_node_structured_mesh_3d_silo(const char *const type,const Patch_T *const patch,double *const x,double *const y,double *const z);
static void pr_multi_mesh(Pr_Field_T *const pr);
