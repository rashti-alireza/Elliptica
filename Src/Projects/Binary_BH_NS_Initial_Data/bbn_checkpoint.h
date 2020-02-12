#include "core_lib.h"
#include "memory_managing_lib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAX_ARR 500
#define checkpoint_file_name "checkpoint.dat"
#define ALLOC_HEADER "#{ALLOCATION#"
#define ALLOC_FOOTER "#}ALLOCATION#"
#define PARAM_HEADER "#{PARAM#"
#define PARAM_FOOTER "#}PARAM#"
#define FIELD_HEADER "#{FIELD#"
#define FIELD_FOOTER "#}FIELD#"

/* this is how we write binary data: first write size and then value. 
// thus, when we wanna read the data the first one gives of the memory allocation and the next gives us value */
#define Write(x,y) if (x) {unsigned SIZE_ = y; fwrite(&SIZE_,sizeof(SIZE_),1,file);fwrite(x,sizeof(*x),SIZE_,file);}\
                   else   {unsigned SIZE_ = 0; fwrite(&SIZE_,sizeof(SIZE_),1,file);}

/* read pointer */
#define ReadP(x,y)  {unsigned SIZE_ = 0; fread(&SIZE_, sizeof(SIZE_),1,file); y = SIZE_;\
                     x = calloc(SIZE_,sizeof(*x)); fread(x,sizeof(*x),SIZE_,file);{if (!y) x = 0;}}

/* read variable */
#define ReadV(x,y)  {unsigned SIZE_ = 0; fread(&SIZE_, sizeof(SIZE_),1,file); y = SIZE_;\
                     fread(x,sizeof(*x),SIZE_,file);}
                     
extern Grid_T **grids_global;
extern Parameter_T **parameters_global;

void bbn_write_checkpoint(const Grid_T *const grid);
Grid_T *bbn_read_checkpoint(void);
static void move_checkpoint_file(void);
static void write_parameters(const Grid_T *const grid);
static void write_fields(const Grid_T *const grid);
static void write_header(const Grid_T *const grid);
