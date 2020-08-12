#include "core_lib.h"
#include "error_handling_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "fields_lib.h"
#include <sys/resource.h>

#define EPS 1E-11
#define UF_OpenMP(x) _Pragma(#x)
#define MAX_STR_LEN 1000

void test_start(const char *const file,const int line);
unsigned countf(void *const p);
unsigned L(const unsigned *const n, const unsigned i, const unsigned j, const unsigned k);
unsigned I(const unsigned l, const unsigned *const n);
unsigned J(const unsigned l, const unsigned *const n);
unsigned K(const unsigned l, const unsigned *const n);
Collocation_T get_collocation(const char *const coll);
Basis_T get_basis(const char *const basis);
int IsOnEdge(const unsigned *const n,const unsigned p);
int IsOnFace(const double *const x, const Patch_T *const patch,unsigned *const f);
unsigned node_onFace(const double *const x, const unsigned f,const Patch_T *const patch);
void IJK(const unsigned l, const unsigned *const n, unsigned *const i, unsigned *const j, unsigned *const k);
static unsigned check_interface(const double *const X, const Patch_T *const patch, const int u);
SubFace_T *get_paired_subface(const SubFace_T *const sub);
unsigned total_nodes_grid(const Grid_T *const grid);
unsigned total_nodes_patch(const Patch_T *const patch);
Coord_T find_coord(const char *const coordsys);
double random_double(const double initial,const double final,const unsigned s);
void copy_subface(SubFace_T *const s2,const SubFace_T *const s1);
unsigned subface_map_invers_id(const SubFace_T *const subface,const unsigned n);
unsigned *dup_UINT(const unsigned *const s,const unsigned N);
double max_Jacobian_dX_dx(Patch_T *const patch);
double spectral_derivative_max_error(const Field_T *const f,const unsigned o);
unsigned IsItHorizonPatch(const Patch_T *const patch);
unsigned IsItNSPatch(const Patch_T *const patch);
unsigned IsItNSSurroundingPatch(const Patch_T *const patch);
void dbprint(const double *v,const unsigned n,const char *const desc);
Patch_T *GetPatch(const char *const stem,const Grid_T *const grid);
unsigned IsItNSSurface(const Patch_T *const patch);
unsigned IsItOutermostPatch(const Patch_T *const patch);
unsigned IsItFillingBoxPatch(const Patch_T *const patch);
unsigned IsItInsideBHPatch(const Patch_T *const patch);
double spectral_expansion_truncation_error(Field_T *const f);
void print_spectral_expansion_truncation_error(Grid_T *const grid);
void shell_command(const char *const cmd);
void free_2d_mem(void *mem0, const unsigned long c);
void free_2d(void *mem0);
double **alloc_2D_double(const long unsigned R,const long unsigned C);
double *alloc_double(const unsigned N);
void _free(void *p);
unsigned IsItFarthestOutermostPatch(const Patch_T *const patch);

double 
how_much_memory
  (
    const char *const unit/* gb,mb,kb */
  );





