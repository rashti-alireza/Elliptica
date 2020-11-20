#ifndef utilities_LIB_H
#define utilities_LIB_H


#include "manifold_lib.h"

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
  unsigned flg: 1;/* used for different purposes */
}sFunc_Patch2Pdouble_T;

void test_start(const char *const file,const int line);
unsigned countf(void *const p);
void init_func_PtoV(sFunc_PtoV_T ***const func);
void add_func_PtoV(sFunc_PtoV_T ***const func,void (*f)(Patch_T *const patch),const char *const task,const Coord_T coord);
void run_func_PtoV(sFunc_PtoV_T **const func,const char *const task,Patch_T *const patch);
Collocation_T get_collocation(const char *const coll);
Coord_T find_coord(const char *const coordsys);
Basis_T get_basis(const char *const basis);
void IJK(const unsigned l, const unsigned *const n, unsigned *const i, unsigned *const j, unsigned *const k);
unsigned L(const unsigned *const n, const unsigned i, const unsigned j, const unsigned k);
unsigned I(const unsigned l, const unsigned *const n);
unsigned J(const unsigned l, const unsigned *const n);
unsigned K(const unsigned l, const unsigned *const n);
int IsOnEdge(const unsigned *const n,const unsigned p);
int IsOnFace(const double *const x, const Patch_T *const patch,unsigned *const f);
SubFace_T *get_paired_subface(const SubFace_T *const sub);
unsigned total_nodes_grid(const Grid_T *const grid);
unsigned total_nodes_patch(const Patch_T *const patch);
void init_func_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const func);
void add_func_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const func,double *(*f)(Patch_T *const patch),const char *const name);
sFunc_Patch2Pdouble_T *get_func_Patch2Pdouble(const char *const name,sFunc_Patch2Pdouble_T **const func);
double random_double(const double initial,const double final,const unsigned s);
void copy_subface(SubFace_T *const s2,const SubFace_T *const s1);
unsigned subface_map_invers_id(const SubFace_T *const subface,const unsigned n);
unsigned *dup_UINT(const unsigned *const s,const unsigned N);
double max_Jacobian_dX_dx(Patch_T *const patch);
double spectral_derivative_max_error(const struct FIELD_T *const f,const unsigned o);
void dbprint(const double *v,const unsigned n,const char *const desc);
Patch_T *GetPatch(const char *const stem,const Grid_T *const grid);
double spectral_expansion_truncation_error(struct FIELD_T *const f);
void print_spectral_expansion_truncation_error(Grid_T *const grid);
void shell_command(const char *const cmd);
void free_2d_mem(void *mem0, const unsigned long c);
void free_2d(void *mem0);
double **alloc_2D_double(const long unsigned R,const long unsigned C);
double *alloc_double(const unsigned N);
void *alloc_sFunc_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const mem);
void free_func_PtoV(sFunc_PtoV_T **func);
void _free(void *p);
void *alloc_sFunc_PtoV(sFunc_PtoV_T ***const mem);

double 
how_much_memory
  (
    const char *const unit/* gb,mb,kb */
  );

void header_and_clock(const char *const msg);
void footer_and_clock(const char *const msg);

void 
update_parameters_and_directories
  (
   const unsigned main_loop_iter,
   const char *const dir_name_format/* eg: "BBN_%s_%ux%ux%u" */
  );

void free_grid_and_its_parameters(Grid_T *grid,const int keep_grid);

#endif




