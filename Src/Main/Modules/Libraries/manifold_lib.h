struct Collocation_s
{
  double min;
  double max;
  unsigned n;
  double stp;
  double a;
  double b;
  Collocation_T c;
};

double point_value(const unsigned i, const struct Collocation_s *const coll_s);
void initialize_collocation_struct(const Patch_T *const patch,struct Collocation_s *const colloc,const unsigned dir);
int make_patches(Grid_T *const grid);
int realize_geometry(Grid_T *const grid);
int X_of_x(double *const X,const double *const x,const Patch_T *const patch);
int x_of_X(double *const x,const double *const X,const Patch_T *const patch);
double x_coord(const unsigned i,const Patch_T *const patch);
double y_coord(const unsigned i,const Patch_T *const patch);
double z_coord(const unsigned i,const Patch_T *const patch);
double X_coord(const unsigned i,const Patch_T *const patch);
double Y_coord(const unsigned i,const Patch_T *const patch);
double Z_coord(const unsigned i,const Patch_T *const patch);
double dq2_dq1(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double JT_Cartesian_patch(Patch_T *const patch,const Dd_T q2_e, const Dd_T q1_e,const unsigned p);
double *normal_vec(Point_T *const point);
double General2ChebyshevExtrema(const double X,const unsigned dir,const Patch_T *const patch);
void grid_characteristics_example(Grid_T *const grid);
double R_interpolation_CS(Field_T *const R,const double *const X);
void SignAndIndex_permutation_CubedSphere(const Flag_T side,unsigned *const a,unsigned *const b,unsigned *const c,double *const s);
void test_CubedSpherical_Coordinates(Grid_T *const grid);
void test_dq_dN(Grid_T *const grid);
void needle_ex(Needle_T *const needle,const Patch_T *const patch);
void needle_in(Needle_T *const needle,const Patch_T *const patch);
void point_finder(Needle_T *const needle);
void populate_right_BH(Grid_T *const grid,const unsigned pn);
void populate_right_BH_central_box(Grid_T *const grid,const unsigned pn);
int make_nodes(Grid_T *const grid);
void alloc_nodes(Grid_T *const grid);
void alloc_interface(Patch_T *const patch);
void *alloc_point(const unsigned s);
void free_points(Grid_T *const grid);
void free_patch_interface(Patch_T *const patch);
void free_needle(Needle_T *needle);
void *alloc_needle(void);
Patch_T make_temp_patch(const Patch_T *const patch);
void free_temp_patch(Patch_T *const patch);
void *alloc_grid(void);
void alloc_patches(Grid_T *const grid);
void free_grid(Grid_T *grid);
void free_patch(Patch_T *patch);
void free_grid_db(void);





