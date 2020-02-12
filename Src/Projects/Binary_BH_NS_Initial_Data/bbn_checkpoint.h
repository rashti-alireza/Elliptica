#include "core_lib.h"
#include "memory_managing_lib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAX_ARR 500
#define checkpoint_file_name "checkpoint.dat"
#define HEADER_DONE "#header_done#"

/* write macro: if the quantity exists put its size first and then the content otherwise put its size zero */
#define Write(x,y) if (x) {size_t size = y; fwrite(&size,sizeof(size),1,file);fwrite(x,sizeof(*x),size,file);}\
                   else   {size_t size = 0; fwrite(&size,sizeof(size),1,file);}

/* read variable if any */ 
#define Read(x,y)   {unsigned Exsits = 0;\
                     fread(&Exsits, sizeof(exsits),1,file);\
                     if (Exsits){fread(&x,Exsits,1,file);}\
                     else       {x = 0;}}

extern Grid_T **grids_global;
extern Parameter_T **parameters_global;

void bbn_write_checkpoint(const Grid_T *const grid);
Grid_T *bbn_read_checkpoint(void);
static void move_checkpoint_file(void);
static void write_parameters(const Grid_T *const grid);
static void write_fields(const Grid_T *const grid);
static void write_header(const Grid_T *const grid);
