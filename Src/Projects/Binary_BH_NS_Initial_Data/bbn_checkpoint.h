#include "bbn_headers.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAX_ARR 500
#define ALLOC_HEADER "#{ALLOCATION#"
#define ALLOC_FOOTER "#}ALLOCATION#"
#define PARAM_HEADER "#{PARAM#"
#define PARAM_FOOTER "#}PARAM#"
#define FIELD_HEADER "#{FIELD#"
#define FIELD_FOOTER "#}FIELD#"
#define END_MSG      "\n#checkpoint_file_completed#\n"

/* this is how we write binary data: first write size and then value. 
// thus, when we wanna read the data the first one gives of the memory allocation and the next gives us value */
#define Write(x,y) \
if (x){\
  unsigned SIZE_ = (unsigned)(y);\
  assert(fwrite(&SIZE_,sizeof(SIZE_),1,file));\
  assert(fwrite(x,sizeof(*(x)),SIZE_,file));\
}else{\
  unsigned SIZE_ = 0;\
  assert(fwrite(&SIZE_,sizeof(SIZE_),1,file));\
}

/* read pointer */
#define ReadP(x) {\
  unsigned SIZE_ = 0;\
  assert(fread(&SIZE_, sizeof(SIZE_),1,file));\
  if (SIZE_) {\
    x = calloc(SIZE_,sizeof(*(x))),pointerEr(x);\
    assert(fread(x,sizeof(*(x)),SIZE_,file));}\
  else { x = 0;}}

/* read variable */
#define ReadV(x) {\
  unsigned SIZE_ = 0;\
  assert(fread(&SIZE_, sizeof(SIZE_),1,file));\
  assert(fread(x,sizeof(*(x)),SIZE_,file));}
                     
extern Grid_T **grids_global;
extern Parameter_T **parameters_global;

struct checkpoint_header
{
 Grid_T *grid;
 //unsigned npatch;
 unsigned npar;
 unsigned grid_number;
 char *grid_kind;
};

void bbn_write_checkpoint(const Grid_T *const grid);
Grid_T *bbn_initi_from_checkpoint(FILE *const file);
static void move_checkpoint_file(void);
static void write_parameters(const Grid_T *const grid);
static void write_fields(const Grid_T *const grid);
static void write_header(const Grid_T *const grid);
static void read_parameters(struct checkpoint_header *const alloc_info,FILE *const file);
static void read_fields(struct checkpoint_header *const alloc_info,FILE *const file);
static void read_header(struct checkpoint_header *const alloc_info,FILE *const file);
static void alloc_db(struct checkpoint_header *const alloc_info);
static int DoSaveField(const Field_T *const field);
static void init_mediate_field(Grid_T *const grid);
Parameter_T *bbn_parameter_query_from_checkpoint_file(const char *const par_name,FILE *const file);
int bbn_IsCheckpointFileCompleted(const char *const file_path);
