#include "bbn_headers.h"

#define STR_LEN_MAX 1000
#define HEADER "#{data#"
#define FOOTER "#}data#"
#define END_MSG "\n#file_completed#\n"


/* strcut for point where interpolate taken place */
struct interpolation_points
{
  double *x,*y,*z;/* (x,y,z) coords */
  double *X,*Y,*Z;/* (X,Y,Z) coords */
  unsigned *patchn;/* patch number for each coord */
  unsigned npoints;/* number of coords */
};

void bbn_bam_export_id(void);
static void load_coords_from_coords_file(struct interpolation_points *const pnt);
//static Grid_T *load_grid_from_checkpoint_file(void);
//static void interpolate_and_write(Grid_T *const grid);


