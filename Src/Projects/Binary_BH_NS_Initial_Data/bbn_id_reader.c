/*
// Alireza Rashti
// Feb 2019
*/

#include "bbn_id_reader.h"

/* local variables to this file */
static char coords_file_path[STR_LEN_MAX];
static char fields_file_path[STR_LEN_MAX];
static char checkpoint_path[STR_LEN_MAX];
static char fields_name[STR_LEN_MAX];
  
/* exporting initial data for bam.
// it writes the required fields into a file to be read by bam. */
void bbn_bam_export_id(void)
{
  /* print some description */
  pr_clock();
  pr_line_custom('=');
  printf("{ Exporting initial data for BAM ...\n");
  
  Grid_T *grid = 0;
  struct interpolation_points points[1] = {0};
  
  /* write parameters b/c they'll be overridden 
  // when checkpoint file is loaded */
  sprintf(coords_file_path,"%s",Pgets("bam_bbn_coords_file_path"));
  sprintf(fields_file_path,"%s",Pgets("bam_bbn_fields_file_path"));
  sprintf(fields_name,     "%s",Pgets("bam_bbn_fields_name"));
  sprintf(checkpoint_path, "%s",Pgets("bam_bbn_checkpoint_path"));
  
  /* read (x,y,x) points from bam file to be interpolated on them */
  load_coords_from_coords_file(points);
  
  /* load grid from the checkpoint file */
  //grid = load_grid_from_checkpoint_file();
  
  /* interpolate at the points and write into fields_file */
  //interpolate_and_write(grid);
  
  /* free */
  free_grid(grid);
  _free(points->x);
  _free(points->y);
  _free(points->z);
  _free(points->X);
  _free(points->Y);
  _free(points->Z);
  _free(points->patchn);
  
  /* print some description */
  printf("} Exporting initial data for BAM ==> Done. :)\n");
  pr_clock();
  pr_line_custom('=');
}

/* read (x,y,x) points from bam file to be interpolated on them,
// and then populate (x,y,z) part of for the interpolation */
static void load_coords_from_coords_file(struct interpolation_points *const pnt)
{
  FILE *file = 0;
  unsigned npoints = 0;
  char *match_str  = 0;
  char str[STR_LEN_MAX] = {'\0'};
  unsigned i = 0;
  
  /* open and read coords file */
  file = fopen(coords_file_path,"r");
  pointerEr(file);
  
  /*winding file */
  fgets(str,STR_LEN_MAX,file);
  FReadP_bin(match_str)
  if (strcmp(match_str,HEADER))
    abortEr("It could not find the header.\n");
  _free(match_str);
  
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
  pointerEr(pnt->patchn);
  /* reading (x,y,z) */
  for (i = 0; i < npoints; ++i)
  {
    FReadV_bin(pnt->x[i]);
    FReadV_bin(pnt->y[i]);
    FReadV_bin(pnt->z[i]);
  }
  FReadP_bin(match_str)
  if (strcmp(match_str,FOOTER))
    abortEr("It could not find the footer.\n");
  _free(match_str);
  
  fclose(file);
}


