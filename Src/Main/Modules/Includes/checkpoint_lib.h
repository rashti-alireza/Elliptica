#ifndef checkpoint_LIB_H
#define checkpoint_LIB_H

void write_checkpoint(Grid_T *const grid);
void read_fields_from_checkpoint(Grid_T *const grid,FILE *const file);
Grid_T *initi_from_checkpoint(FILE *const file);
Parameter_T *parameter_query_from_checkpoint(const char *const par_name,FILE *const file);
int is_checkpoint_sound(const char *const file_path);
Grid_T *load_checkpoint_file(void);
int can_we_use_checkpoint(void);

#endif




