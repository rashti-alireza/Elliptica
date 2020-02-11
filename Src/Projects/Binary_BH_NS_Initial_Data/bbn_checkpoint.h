#include "core_lib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAX_ARR 500
/* write macro: if the quantity exists put its size first and then the content otherwise put its size zero */
#define Write(x,y) if (x) {size_t size = y; fwrite(&size,sizeof(size),1,file);fwrite(x,sizeof(*x),size,file);}\
                   else   {size_t size = 0; fwrite(&size,sizeof(size),1,file);}

extern Parameter_T **parameters_global;

void bbn_write_checkpoint(const Grid_T *const grid);
static void move_checkpoint_file(const char *const file_name);
static void write_parameters(const Grid_T *const grid,const char *const file_name);
static void write_fields(const Grid_T *const grid,const char *const file_name);
static void write_header(const Grid_T *const grid,const char *const file_name);
