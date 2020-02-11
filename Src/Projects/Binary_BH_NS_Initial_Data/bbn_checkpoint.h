#include "core_lib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAX_ARR 500

extern Parameter_T **parameters_global;

void bbn_write_checkpoint(const Grid_T *const grid);
static void move_checkpoint_file(const char *const file_name);
static void write_parameters(const Grid_T *const grid,const char *const file_name);
static void write_fields(const Grid_T *const grid,const char *const file_name);
