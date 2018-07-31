#include <silo.h>
#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "text_tools_lib.h"
#include "macros_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"

#define MAX_STR_LEN 500

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
