void *alloc_parameter(Parameter_T ***mem);
void *alloc_project(Project_T ***mem);
void *alloc_grid(void);
void alloc_patches(Grid_T *grid);
void alloc_nodes(Grid_T *grid);
void free_2d(void **mem);
void free_matrix(void **mem, long int c);
