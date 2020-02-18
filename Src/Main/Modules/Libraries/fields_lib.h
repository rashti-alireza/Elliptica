Field_T *add_field(const char *const name,const char *attribute,Patch_T *const patch,const Flag_T alloc_flg);
void remove_field(Field_T *f);
void add_attribute(Field_T *const fld,const char *const attribute);
int LookUpField(const char *const name,const Patch_T *const patch);
int LookUpField_E(const char *const name,const Patch_T *const patch);
double *make_coeffs_1d(Field_T *const f,const unsigned dir);
double *make_coeffs_2d(Field_T *const f,const unsigned dir1,const unsigned dir2);
double *make_coeffs_3d(Field_T *const f);
void enable_fields(Grid_T *const grid);
void free_v2(Field_T *f);
void free_attr(Field_T *f);
void free_v(Field_T *f);
void free_info(Field_T *f);
void free_field(Field_T *fld);
void empty_field(Field_T *fld);
void free_coeffs(Field_T *fld);





