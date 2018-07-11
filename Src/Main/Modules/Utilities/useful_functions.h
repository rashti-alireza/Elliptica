#include "core_lib.h"
#include "error_handling_lib.h"
#include "coordinates_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"

#define EPS 1E-11

void test_start(const char *const file,const int line);
unsigned countf(void *const p);
unsigned L(const unsigned *const n, const unsigned i, const unsigned j, const unsigned k);
unsigned I(const unsigned l, const unsigned *const n);
unsigned J(const unsigned l, const unsigned *const n);
unsigned K(const unsigned l, const unsigned *const n);
Collocation_T get_collocation(const char *const coll);
int IsOnEdge(const unsigned *const n,const unsigned p);
int IsOnFace(const double *const x, const Patch_T *const patch,unsigned *const f);
unsigned node_onFace(const double *const x, const unsigned f,const Patch_T *const patch);
void IJK(const unsigned l, const unsigned *const n, unsigned *const i, unsigned *const j, unsigned *const k);
static unsigned check_interface(const double *const X, const Patch_T *const patch, const int u);
SubFace_T *get_paired_subface(const SubFace_T *const sub);
unsigned total_nodes_grid(const Grid_T *const grid);
unsigned total_nodes_patch(const Patch_T *const patch);
