/* allocating */
void *alloc_parameter(Parameter_T ***const mem);
void *alloc_project(Project_T ***const mem);
void *alloc_grid(void);
void alloc_patches(Grid_T *const grid);
void alloc_nodes(Grid_T *const grid);
void *alloc_point(const unsigned s);
void alloc_interface(Patch_T *const patch);
void *alloc_sFunc_PtoV(sFunc_PtoV_T ***const mem);
void *alloc_needle(void);
void *alloc_sFunc_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const mem);
double *alloc_double(const unsigned N);
double **alloc_2D_double(const long unsigned R,const long unsigned C);
Matrix_T *alloc_matrix(const Matrix_SF_T type_e,const long row,const long col);
Solve_T *alloc_solve(Patch_T *const patch, const unsigned n);

/* free */
void free_needle(Needle_T *needle);
void free_2d(void *mem0);
void free_2d_mem(void *mem, const unsigned long c);
void free_points(Grid_T *const grid);
void free_func_PtoV(sFunc_PtoV_T **func);
void free_field(Field_T *f);
void free_func_Patch2Pdouble(sFunc_Patch2Pdouble_T **func);
void free_v2(Field_T *f);
void free_v(Field_T *f);
void free_info(Field_T *f);
void free_attr(Field_T *f);
void free_coeffs(Field_T *fld);
void free_db_eqs(sEquation_T **db);
void free_interpolation(Interpolation_T *interp_s);
void free_matrix(Matrix_T *m);
