#ifndef TOV_LIB_H
#define TOV_LIB_H

/* forward declaration */
struct PHYSICS_T;

/* struct for TOV stars */
typedef struct TOV_PROJECT_T
{
 struct PHYSICS_T *phys;/* physics */
 const char *description;
 double ADM_m;/* ADM mass of NS */
 double bar_m;/* baryonic mass of NS */
 double h_cent;/* enthalpy at the center of NS */
 Uint N;/* number of interpolation points, choose them to be odd */
 double *m;/* total mass at each point */
 double *r;/* radius at each point for metric inside the star:
           // ds^{2}=-e^{2\phi \left( r\right) }dt^{2}+\left( 1-\frac {2m\left( r\right) }{r}\right) ^{-1}dr^{2}+r^{2}\left( d\theta^{2}+sin^{2} \left(\theta \right)d\phi^{2}) */
 double *p;/* pressure at each point */
 double *h;/* enthalpy at each point */
 double *phi;/* at the spacetime metric g_00 = - exp[2phi] */
 double *rbar; /* rbar -> radius at each point for metric inside the star:
             // d^{2}_{s}=-e^{2\phi }dt^{2}+\psi^{4}\left( \overline {r}\right) \left( d\overline {r}^{2}+\overline {r}^{2}d\Omega \right) */
 double *psi;/* conformal factor */
 Uint status:1;/* status of exist: 0 = no error happens, 1 = some errors happened. */
 Uint exit_if_error:1;/* 1 = if error happens exit Elliptica (default), 
                          // 0 = if error happens don't exit Elliptica  */
}TOV_T;

TOV_T *TOV_solution(TOV_T *const TOV);
void TOV_free(TOV_T *TOV);
TOV_T *TOV_init(void);


#endif



