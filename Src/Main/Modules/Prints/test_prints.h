#include "core_lib.h"
#include "error_handling_lib.h"
#include "print_lib.h"
#include "utilities_lib.h"
#include "text_tools_lib.h"
#include "macros_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"

extern Parameter_T **parameters_global;

struct Archive_S
{
  SubFace_T *s1;
  SubFace_T *s2;
  char *n1;
  char *n2;
};

int test_print(const Print_T f);
void pr_parameters(void);
void pr_coords(const Grid_T *const grid);
void pr_interfaces(const Grid_T *const grid);
static void add_to_archive(struct Archive_S **const arch,SubFace_T *const s1,SubFace_T *const s2,unsigned *const n,const char *const desc);
static void free_archive(struct Archive_S *arch,const unsigned Narch);
void pr_derivatives_DiffByNode(const double *const numc, const double *const anac,const Patch_T *const patch,const char *const prefix);
