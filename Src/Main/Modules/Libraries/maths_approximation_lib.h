void rft_1d_ChebyshevExtrema_coeffs(double *const values ,double *const coeffs,const unsigned n);
void rft_1d_ChebyshevNodes_coeffs(double *const values ,double *const coeffs,const unsigned n);
void *r2cft_1d_EquiSpaced_coeffs(const double *const value,const unsigned n);
double *c2rft_1d_EquiSpaced_values(void *const coeffs,const unsigned N);
int derivative_tests(Grid_T *const grid);
int interpolation_tests(Grid_T *const grid);
int fourier_transformation_tests(Grid_T *const grid);
int Ylm_transformation_tests(Grid_T *const grid);
Interpolation_T *init_interpolation(void);
double execute_interpolation(Interpolation_T *const interp_s);
void plan_interpolation(Interpolation_T *const interp_s);
void get_Ylm_coeffs(double *const realClm,double *const imagClm,const double *const f,const unsigned Ntheta,const unsigned Nphi,const unsigned Lmax);
double interpolation_Ylm(const double *const realClm,const double *const imagClm,const unsigned Lmax, const double theta, const double phi);
double *df_dphi_Ylm(const double *const realClm,const double *const imagClm,const unsigned Ntheta, const unsigned Nphi,const unsigned Lmax);
double *df_dtheta_Ylm(const double *const realClm,const double *const imagClm,const unsigned Ntheta, const unsigned Nphi,const unsigned Lmax);
unsigned lm2n(const unsigned l, const unsigned m);
double *alloc_ClmYlm(unsigned Lmax);


