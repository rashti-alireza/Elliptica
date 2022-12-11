#ifndef utilities_LIB_H
#define utilities_LIB_H
#include "elliptica_system_lib.h"


#include "manifold_lib.h"

/* tuple (i,j) format to linear format (row-major order). */
#define i_j_to_ij_row_major_order(nj,i,j) ((j)+(nj)*(i))

/* triple (i,j,k) format to linear format (row-major order). */
#define i_j_k_to_ijk_row_major_order(n,i,j,k) \
 ((k)+((n)[2])*((j)+((n)[1])*(i)))

/* converting linear format ijk to tuple (i,j,k) format */
#define ijk_to_i_j_k(ijk,n,i,j,k)  \
 (ijk_to_i_j_k_row_major_order((ijk),(n),(i),(j),(k)))

/* converting tuple (i,j,k) format to linear format ijk */
#define i_j_k_to_ijk(n,i,j,k)  \
 (i_j_k_to_ijk_row_major_order((n),(i),(j),(k)))

/* converting tuple (i,j) format to linear format ij */
#define i_j_to_ij(nj,i,j)  \
 (i_j_to_ij_row_major_order((nj),(i),(j)))

/* converting a symmetric tuple (i,j) format to a liner format with dimension n. */
#define i_j_to_ij_symmetric(n,i,j) i_j_to_ij_symmetric_row_major_order(n, i, j)

/* converting a symmetric tuple (i,j) format to a liner format with dimension n
// row major order. NOTE: we need second index, here j, to satisfy: j >= i >= 0. 
// in order to make the usage easy for the user, we impose this condition in the macro.
// example:
// the indices xx, xy, xz, yy, yz, and zz map to 0, 1, 2, 3, 4, and 5 (here n = 3). */
#define i_j_to_ij_symmetric_row_major_order(n,i,j) \
( j > i ? ((n)*(i) - (i)*((i)-1)/2 + ((j)-(i))) : ((n)*(j) - (j)*((j)-1)/2 + ((i)-(j))) )

#define TEST_START test_start(__FILE__,__LINE__);


/* forward declaration structures */
struct FIELD_T;

/* general function patch to void */
typedef void fFunc_PtoV_T (Patch_T *const patch);

/* general function grid to pointer to double */
typedef double *fFunc_Patch2Pdouble_T(Patch_T *const patch);

/* patch to void struct */
typedef struct sFUNC_PtoV_T
{
  char *task;
  Coord_T coord;
  fFunc_PtoV_T *f;
}sFunc_PtoV_T;

/* grid to double struct */
typedef struct sFUNC_PATCH2PDOUBLE_T
{
  char *name;
  fFunc_Patch2Pdouble_T *func;
  Uint flg: 1;/* used for different purposes */
}sFunc_Patch2Pdouble_T;

void test_start(const char *const file,const int line);
Uint countf(void *const p);
void init_func_PtoV(sFunc_PtoV_T ***const func);
void add_func_PtoV(sFunc_PtoV_T ***const func,void (*f)(Patch_T *const patch),const char *const task,const Coord_T coord);
void run_func_PtoV(sFunc_PtoV_T **const func,const char *const task,Patch_T *const patch);
Collocation_T get_collocation(const char *const coll);
Coord_T find_coord(const char *const coordsys);
Basis_T get_basis(const char *const basis);
INLINE void ijk_to_i_j_k_row_major_order(const Uint l, const Uint *const n, Uint *const i, Uint *const j, Uint *const k) INLINE_WARN_UNUSED_FUNC;
Uint ijk_to_i_row_major_order(const Uint l, const Uint *const n);
Uint ijk_to_j_row_major_order(const Uint l, const Uint *const n);
Uint ijk_to_k_row_major_order(const Uint l, const Uint *const n);
int IsOnEdge(const Uint *const n,const Uint p);
int IsOnFace(const double *const x, const Patch_T *const patch,Uint *const f,const double precision_factor);
SubFace_T *get_paired_subface(const SubFace_T *const sub);
Uint total_nodes_grid(const Grid_T *const grid);
Uint total_nodes_patch(const Patch_T *const patch);
void init_func_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const func);
void add_func_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const func,double *(*f)(Patch_T *const patch),const char *const name);
sFunc_Patch2Pdouble_T *get_func_Patch2Pdouble(const char *const name,sFunc_Patch2Pdouble_T **const func);
double random_double(const double initial,const double final,const Uint s);
void copy_subface(SubFace_T *const s2,const SubFace_T *const s1);
Uint subface_map_invers_id(const SubFace_T *const subface,const Uint n);
Uint *dup_UINT(const Uint *const s,const Uint N);
double max_Jacobian_dX_dx(Patch_T *const patch);
double spectral_derivative_max_error(const struct FIELD_T *const f,const Uint o);
void dbprint(const double *v,const Uint n,const char *const desc);
Patch_T *GetPatch(const char *const stem,const Grid_T *const grid);
double spectral_expansion_truncation_error(struct FIELD_T *const f);
void print_spectral_expansion_truncation_error(Grid_T *const grid);
void shell_command(const char *const cmd);
void free_2d_mem(void *mem0, const Uint long c);
void free_2d(void *mem0);
double **alloc_2D_double(const long Uint R,const long Uint C);
double *alloc_double(const Uint N);
void *alloc_sFunc_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const mem);
void free_func_PtoV(sFunc_PtoV_T **func);
void *alloc_sFunc_PtoV(sFunc_PtoV_T ***const mem);

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

void add_3x3_symmetric_field(Grid_T *const grid,const char *const stem,const char *rank);

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

/* inline function definitions */

/* linear format to triple (i,j,k) format (row-major order). */
INLINE void ijk_to_i_j_k_row_major_order(const Uint l, const Uint *const n, Uint *const i, Uint *const j, Uint *const k)
{
  Uint tmp;
  
  tmp = l % (n[2]*n[1]);
  *i  = l / (n[2]*n[1]);
  *j  = tmp / n[2];
  *k  = tmp % n[2];
}



#endif




