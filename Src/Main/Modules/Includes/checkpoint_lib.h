#ifndef checkpoint_LIB_H
#define checkpoint_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct PHYSICS_T;

void write_checkpoint(struct PHYSICS_T *const phys,const char *const out_dir);
void read_fields_from_checkpoint(Grid_T *const grid,FILE *const file);
Grid_T *init_from_checkpoint(FILE *const file);
Parameter_T *parameter_query_from_checkpoint(const char *const par_name,FILE *const file);
int is_checkpoint_sound(const char *const file_path);
Grid_T *load_checkpoint_file(void);
int can_we_use_checkpoint(const char *const cur_out_dir);

#endif




