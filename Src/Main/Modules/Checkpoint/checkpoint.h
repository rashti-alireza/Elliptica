#include "headers.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "text_and_file_tools_lib.h"

#define MAX_ARR (500)
#define ALLOC_HEADER "#{ALLOCATION#"
#define ALLOC_FOOTER "#}ALLOCATION#"
#define PARAM_HEADER "#{PARAM#"
#define PARAM_FOOTER "#}PARAM#"
#define FIELD_HEADER "#{FIELD#"
#define FIELD_FOOTER "#}FIELD#"
#define END_MSG      "\n#checkpoint_file_completed#\n"

                     
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

void write_checkpoint(Grid_T *const grid);
void read_fields_from_checkpoint(Grid_T *const grid,FILE *const file);
Grid_T *initi_from_checkpoint(FILE *const file);
Parameter_T *parameter_query_from_checkpoint(const char *const par_name,FILE *const file);
int is_checkpoint_sound(const char *const file_path);
static void move_checkpoint_file(void);
static void write_parameters(const Grid_T *const grid);
static void write_fields(const Grid_T *const grid);
static void write_header(const Grid_T *const grid);
static void read_parameters(struct checkpoint_header *const alloc_info,FILE *const file);
static void read_header(struct checkpoint_header *const alloc_info,FILE *const file);
static void alloc_db(struct checkpoint_header *const alloc_info);
static void find_and_save_modified_checkpoint_pars(void);
static void incorporate_modified_checkpoint_par(void);
static void free_modified_checkpoint_par(void);


