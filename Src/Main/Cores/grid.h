#include "core_lib.h"
#include "error_handling_lib.h"
#include "memory_managing_lib.h"

/* returning value */
struct Ret_S
{
  char s0[20],s1[20],s2[20];
};


static void fill_patches(Grid_T *grid);
static void fill_patches_cartesian(Grid_T *grid);
static void fill_nodes(Grid_T *grid);
static void fill_nodes_cartesian(Grid_T *grid);
static void make_keyword_parameter(struct Ret_S *ret,char *box,char *needle);
int make_coordinates(Grid_T *grid);