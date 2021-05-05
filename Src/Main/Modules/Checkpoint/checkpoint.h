#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include "core_lib.h"
#include "fields_lib.h"
#include "physics_lib.h"
#include "checkpoint_lib.h"


#define MAX_ARR   (400)
#define MAX_ARRx2 (2*MAX_ARR)
#define MAX_ARRx3 (3*MAX_ARR)
#define MAX_ARRx4 (4*MAX_ARR)
#define MAX_ARRx5 (5*MAX_ARR)
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
 //Uint npatch;
 Uint npar;
 Uint grid_number;
 char *grid_kind;
};


void write_checkpoint(Physics_T *const phys,const char *const out_dir);
void read_fields_from_checkpoint_file(Physics_T *const phys,FILE *const file);
void *open_checkpoint_file_then_read_grid_and_params(Physics_T *const phys);
int can_we_use_checkpoint(const char *const cur_out_dir);
static Grid_T *read_grid_and_params(FILE *const file);
Parameter_T *parameter_query_from_checkpoint(const char *const par_name,FILE *const file);
int is_checkpoint_sound(const char *const file_path);
static void move_checkpoint_file(const char *const folder);
static void write_parameters(const Grid_T *const grid,const char *const folder);
static void write_fields(const Grid_T *const grid,const char *const folder);
static void write_header(const Grid_T *const grid,const char *const folder);
static void read_parameters(struct checkpoint_header *const alloc_info,FILE *const file);
static void read_header(struct checkpoint_header *const alloc_info,FILE *const file);
static void alloc_db(struct checkpoint_header *const alloc_info);
static void find_and_save_modified_checkpoint_pars(void);
static void incorporate_modified_checkpoint_par(void);
static void free_modified_checkpoint_par(void);


