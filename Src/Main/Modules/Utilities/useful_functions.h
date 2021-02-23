#include "core_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "fields_lib.h"
#include "maths_spectral_methods_lib.h"

#include <sys/resource.h>

#define EPS (1E-11)
#define UF_OpenMP(x) _Pragma(#x)
#define MAX_STR_LEN (1000)

/* calculate difference between two fields with given index
// NOTE: it is used in function diff_3x3_symmetric_fields. */
#define CalcDiff(index) \
  {\
    double diffval = fabs(diff1_##index[ijk]-diff2_##index[ijk]);\
    max = (max < diffval ? diffval : max);\
    if (pr_points && GRT(diffval,0.))\
    {\
      printf(Pretty0"%s:(%0.2f,%0.2f,%0.2f) = %0.2e\n",\
             patch->name,\
             patch->node[ijk]->x[0],\
             patch->node[ijk]->x[1],\
             patch->node[ijk]->x[2],\
             diffval\
             );\
    }\
  }

void test_start(const char *const file,const int line);
Uint countf(void *const p);
Uint ijk_to_i_row_major_order(const Uint l, const Uint *const n);
Uint ijk_to_j_row_major_order(const Uint l, const Uint *const n);
Uint ijk_to_k_row_major_order(const Uint l, const Uint *const n);
Collocation_T get_collocation(const char *const coll);
Basis_T get_basis(const char *const basis);
int IsOnEdge(const Uint *const n,const Uint p);
int IsOnFace(const double *const x, const Patch_T *const patch,Uint *const f,const double precision_factor);
Uint node_onFace(const double *const x, const Uint f,const Patch_T *const patch);
static Uint check_interface(const double *const X, const Patch_T *const patch, const int u);
SubFace_T *get_paired_subface(const SubFace_T *const sub);
Uint total_nodes_grid(const Grid_T *const grid);
Uint total_nodes_patch(const Patch_T *const patch);
Coord_T find_coord(const char *const coordsys);
double random_double(const double initial,const double final,const Uint s);
void copy_subface(SubFace_T *const s2,const SubFace_T *const s1);
Uint subface_map_invers_id(const SubFace_T *const subface,const Uint n);
Uint *dup_UINT(const Uint *const s,const Uint N);
double max_Jacobian_dX_dx(Patch_T *const patch);
double spectral_derivative_max_error(const Field_T *const f,const Uint o);
void dbprint(const double *v,const Uint n,const char *const desc);
Patch_T *GetPatch(const char *const stem,const Grid_T *const grid);
double spectral_expansion_truncation_error(Field_T *const f);
void print_spectral_expansion_truncation_error(Grid_T *const grid);
void shell_command(const char *const cmd);
void free_2d_mem(void *mem0, const Uint long c);
void free_2d(void *mem0);
double **alloc_2D_double(const long Uint R,const long Uint C);
double *alloc_double(const Uint N);
Uint IsItFarthestOutermostPatch(const Patch_T *const patch);

double 
how_much_memory
  (
    const char *const unit/* gb,mb,kb */
  );


void header_and_clock(const char *const msg);
void footer_and_clock(const char *const msg);
double f_of_X(const char *const field,
              const double *const X/* patch coords */,
              Patch_T *const patch);

double diff_3x3_symmetric_fields(Grid_T *const grid,
                               const char *const stem1/* field1 */,
                               const char *const stem2/* field2 */,
                               const char *const rank/* [up/down] */,
                               const int pr_points/* print all points */);

void superimpose_simple(Grid_T *const grid,
                        const char *const f,
                        const char *const f1,
                        const char *const f2,
                        const double extra);



void interpolate_fields_from_old_grid_to_new_grid
     (Grid_T *const ogrid/* old */,Grid_T *const ngrid/* new */,
     const char *const field_names/* comma separated field names */,
     const int copy/* if 1 only copy, if 0 only 3d interpolation */);



