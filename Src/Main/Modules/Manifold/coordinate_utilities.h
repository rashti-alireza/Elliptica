#include "core_lib.h"
#include "error_handling_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "fields_lib.h"
#include "manifold_header.h"


/* dealing with coordinates and interpolation error between the boundaries
// the followings got from experiment. */
#define EPS_collocation   (1E-8)
#define EPS_coord_general (1E-8)
#define EPS_coord_OB_SCS1 (1E-1)
#define EPS_coord_OB_SCS2 (1E-3)
#define EPS_coord_OT_SCS1 (1E-1)
#define EPS_coord_OT_SCS2 (1E-3)
#define EPS_coord_LOW_n1  (1E0)
#define EPS_coord_LOW_n2  (1E-1)
#define LOW_n             (9)

/* h is an small distance, for instance grid space */
#define LSS_coord(x,y,h)    ((x) < (y)-(h))
#define GRT_coord(x,y,h)    ((x) > (y)+(h))
#define EQL_coord(x,y,h)    (((x) < (y)+(h)) && ((x) > (y)-(h)))
#define LSSEQL_coord(x,y,h) (LSS_coord(x,y,h) || EQL_coord(x,y,h))
#define GRTEQL_coord(x,y,h) (GRT_coord(x,y,h) || EQL_coord(x,y,h))

/* find if x takes place inside the interval [min,max] */
#define IsInside(x,min,max,h) \
  ( LSSEQL_coord(x[0],max[0],h) && GRTEQL_coord(x[0],min[0],h) && \
    LSSEQL_coord(x[1],max[1],h) && GRTEQL_coord(x[1],min[1],h) && \
    LSSEQL_coord(x[2],max[2],h) && GRTEQL_coord(x[2],min[2],h))


/* find a reasonable small distance for this patch and put it in h. */
#define set_h_coord(h,patch,scale) \
{\
  h[0] = scale*(patch->max[0]-patch->min[0])/patch->n[0];\
  h[1] = scale*(patch->max[1]-patch->min[1])/patch->n[1];\
  h[2] = scale*(patch->max[2]-patch->min[2])/patch->n[2];\
}


typedef enum MODE_T
{
  GUESS,
  FORCE_IN
}Mode_T;

/* DON'T MODIFY */
enum Limit
{
  MIN0 = 0,
  MAX0,
  MIN1,
  MAX1,
  MIN2,
  MAX2,
  TOT_Limit
};

void point_finder(Needle_T *const needle);
void needle_ex(Needle_T *const needle,const Patch_T *const patch);
void needle_in(Needle_T *const needle,const Patch_T *const patch);
void needle_guess(Needle_T *const needle,const Patch_T *const patch);
void needle_ans(Needle_T *const needle,const Patch_T *const patch);
static void find(Needle_T *const needle,Mode_T mode);
int X_of_x(double *const X,const double *const x,const Patch_T *const patch);
unsigned find_node(const double *const x, const Patch_T *const patch,Flag_T *const flg);
double x_coord(const unsigned i,const Patch_T *const patch);
double y_coord(const unsigned i,const Patch_T *const patch);
double z_coord(const unsigned i,const Patch_T *const patch);
double X_coord(const unsigned i,const Patch_T *const patch);
double Y_coord(const unsigned i,const Patch_T *const patch);
double Z_coord(const unsigned i,const Patch_T *const patch);
int x_of_X(double *const x,const double *const X,const Patch_T *const patch);
static int x_of_X_CS_coord(double *const x,const double *const X,const Patch_T *const patch,const int check_flg);
static int X_of_x_CS_coord(double *const X,const double *const cart,const Patch_T *const patch,const int check_flg);
static int x_of_X_Cartesian_coord(double *const x,const double *const X,const Patch_T *const patch);
static int X_of_x_Cartesian_coord(double *const X,const double *const x,const Patch_T *const patch);
void free_needle(Needle_T *needle);
void *alloc_needle(void);
void theta_phi_of_XY_CS(double *const theta,double *const phi,const double *const X,const Flag_T side);

int 
IsItCovering
  (
  const Patch_T *const patch,/* the patch */
  const char *const region,/* BH/NS etc. see the list above */
  const Flag_T Fside/* LEFT or RIGHT or CENTER or NONE (side of region, if any) */
  );
  
  

Patch_T **
collect_patches
  (
  Grid_T *const grid,/* the grid */
  const char *const region,/* see the list in IsItCovering function */
  const Flag_T side,/* LEFT or RIGHT or CENTER or NONE */
  unsigned *const Np/* number of patches found */
  );
 

Grid_Char_T *init_grid_char(Grid_T *const last_grid);
void free_grid_char(Grid_Char_T *g);
Grid_Kind_T set_grid_kind(const char *const grid_kind);



