void fftw_1d_ChebyshevExtrema_coeffs(double *const values,double *const coeffs,const unsigned n);
void fftw_1d_ChebyshevExtrema_values(double *const values,double *const coeffs,const unsigned n);
void fftw_3d_ChebyshevExtrema_coeffs(double *const values,double *const coeffs,const unsigned *const n);
void fftw_3d_ChebyshevExtrema_values(double *const values,double *const coeffs,const unsigned *const n);
void rft_1d_ChebyshevExtrema_coeffs(double *const values ,double *const coeffs,const unsigned n);
void rft_1d_ChebyshevNodes_coeffs(double *const values ,double *const coeffs,const unsigned n);
int derivative_tests(Grid_T *const grid);
int interpolation_tests(Grid_T *const grid);
Interpolation_T *init_interpolation(void);
double execute_interpolation(Interpolation_T *const interp_s);
void plan_interpolation(Interpolation_T *const interp_s);

