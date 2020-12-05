/*
// Alireza Rashti
// December 2020
*/

/* reading and writing tools for exporting of initial data */

#include "id_export_read_write.h"


/* given file path, it reads cartesian coordinate from the file
// then populates pnt struct and allocs memory for X coords too. */
void load_Cartesian_coordinates_from_file
                        (const char *const coords_file_path,
                         ID_Export_T *const pnt)
{
  FUNC_TIC
  
  FILE *file = 0;
  Uint npoints = 0;
  char *match_str  = 0;
  char str[STR_LEN_MAX] = {'\0'};
  Uint i = 0;
  
  /* open and read coords file */
  file = Fopen(coords_file_path,"r");
  
  /* winding file */
  fgets(str,STR_LEN_MAX,file);
  FReadP_bin(match_str)
  if (strcmp(match_str,HEADER))
    Error2("It could not find the header");
  Free(match_str);
  
  /* allocating */
  FReadV_bin(npoints)
  pnt->npoints = npoints;
  pnt->x       = alloc_double(npoints);
  pnt->y       = alloc_double(npoints);
  pnt->z       = alloc_double(npoints);
  pnt->X       = alloc_double(npoints);
  pnt->Y       = alloc_double(npoints);
  pnt->Z       = alloc_double(npoints);
  pnt->patchn  = calloc(npoints,sizeof(*pnt->patchn));
  IsNull(pnt->patchn);
  printf("~> number of points to interpolate = %u\n",npoints);
  
  /* reading (x,y,z) */
  for (i = 0; i < npoints; ++i)
  {
    FReadV_bin(pnt->x[i]);
    if(!isfinite(pnt->x[i]))
      Error2("bad coordinate.");
    
    FReadV_bin(pnt->y[i]);
    if(!isfinite(pnt->y[i]))
      Error2("bad coordinate.");
    
    FReadV_bin(pnt->z[i])
    if (!isfinite(pnt->z[i]))
     Error2("bad coordinate.");
  }
  FReadP_bin(match_str)
  if (strcmp(match_str,FOOTER))
    Error2("It could not find the footer.\n");
  Free(match_str);
  
  Fclose(file);
  
  FUNC_TOC
}
