#include "core_lib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

extern Parameter_T **parameters_global;

void bbn_write_checkpoint(const Grid_T *const grid);