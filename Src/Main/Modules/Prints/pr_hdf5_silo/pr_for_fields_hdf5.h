#include <silo.h>
#include "core_lib.h"/* it includes prints_lib.h too */
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "macros_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "fields_lib.h"
#include "../pr_for_fields_structs.h"

#define MAX_STR_LEN 500

Pr_Field_T *init_PrField(const Grid_T *const grid);
void pr_fields(Pr_Field_T *const pr);
void pr_hdf5_silo(Pr_Field_T *const pr);
void free_PrField(Pr_Field_T *pr);
static void free_info_s(Pr_Field_T *const pr);
static void pr_structured_mesh_3d_silo(const Pr_Field_T *const pr);
static void pr_scalar_on_structured_mesh_3d_silo(const Pr_Field_T *const pr);
static void pr_vector_on_structured_mesh_3d_silo(const Pr_Field_T *const pr);
static void *make_structured_mesh_3d_xyz(Pr_Field_T *const pr,const Patch_T *const patch);
static void *make_structured_mesh_3d_abc(Pr_Field_T *const pr,const Patch_T *const patch);
static void prepare_node_structured_mesh_3d_silo(const char *const type,const Patch_T *const patch,double *const x,double *const y,double *const z);
static void write_multi_mesh(const Pr_Field_T *const pr);
static void write_multi_vars(const Pr_Field_T *const pr);
static void make_multi_var(const Pr_Field_T *const pr,const char *const var);

