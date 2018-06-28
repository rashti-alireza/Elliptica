#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"

typedef enum MODE_T
{
  GUESS,
  FORCE_IN
}Mode_T;


static void IsConsistent(Needle_T *needle);
static void look(Needle_T *needle,Mode_T mode);
void needle_in(Needle_T *needle,Patch_T *patch);
