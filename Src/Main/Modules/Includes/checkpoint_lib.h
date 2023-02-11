#ifndef checkpoint_LIB_H
#define checkpoint_LIB_H
#include "elliptica_system_lib.h"

/* the name of checkpoint file */
#define CHECKPOINT_FILE_NAME "checkpoint.dat"

/* to modify or set a parameter in checkpoint file data */
#define CHECKPOINT_SET_PARAM_ "modify:"

/* forward declaration */
struct PHYSICS_T;

void write_checkpoint(struct PHYSICS_T *const phys,const char *const out_dir);
int is_checkpoint_sound(const char *const file_path);
int can_we_use_checkpoint(const char *const cur_out_dir);
void *open_checkpoint_file_then_read_grid_and_params(struct PHYSICS_T *const phys);
void read_fields_from_checkpoint_file(struct PHYSICS_T *const phys,FILE *const file);
Parameter_T *parameter_query_from_checkpoint(const char *const par_name,FILE *const file);

#endif




