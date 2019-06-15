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
double interpolation_2d_PH(Field_T *const R, const Patch_T *const patch,const double *const X);
void grid_characteristics_example(Grid_T *const grid);
double interpolation_2d_PH(Field_T *const R, const Patch_T *const patch,const double *const X);
double interpolation_2d_CS(Field_T *const R, const Patch_T *const patch,const double *const X);
void SignAndIndex_permutation_CubedSphere(const Flag_T side,unsigned *const a,unsigned *const b,unsigned *const c,double *const s);
void test_dNi_dxj(Grid_T *const grid);
