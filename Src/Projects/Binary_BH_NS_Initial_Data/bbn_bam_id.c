/*
// Alireza Rashti
// Feb 2019
*/

#include "bbn_bam_id.h"

/* local variables to this file */
static char coords_file_path[STR_LEN_MAX];
static char fields_file_path[STR_LEN_MAX];
static char checkpoint_path[STR_LEN_MAX];
static char bam_fields_name[STR_LEN_MAX];
  
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
  sprintf(coords_file_path,"%s",Pgets("bbn_bam_coords_file_path"));
  sprintf(fields_file_path,"%s",Pgets("bbn_bam_fields_file_path"));
  sprintf(bam_fields_name, "%s",Pgets("bbn_bam_fields_name"));
  sprintf(checkpoint_path, "%s",Pgets("bbn_bam_checkpoint_path"));
  
  /* read (x,y,x) points from bam file to be interpolated on them */
  load_coords_from_coords_file(points);
  
  /* load grid from the checkpoint file */
  grid = load_grid_from_checkpoint_file();
  
  /* interpolate at the points and write into fields_file */
  interpolate_and_write(grid,points);
  
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
  printf("} Exporting initial data for BAM ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* interpolate at the points and write into fields_file.
// this file will be reading by BAM as the initilization its fields */
static void interpolate_and_write(Grid_T *const grid,struct interpolation_points *const pnt)
{
  FILE *file = 0;
  const unsigned npoints = pnt->npoints;
  char **fields_name = 0,**bam_fields = 0;
  char title_line[STR_LEN_MAX];
  char *const p_title_line = title_line;/* to avoid GCC warning for FWriteP_bin */
  char msg[STR_LEN_MAX];
  char *const p_msg = msg;/* to avoid GCC warning for FWriteP_bin */
  double *interp_v = 0;
  unsigned p,f;
  
  /* populating pnt->(X,Y,Z) and pnt->patchn */
  printf("~> Preparing points for the interpolation ...\n");
  fflush(stdout);
  
  /* to avoid race condition between threads write all coeffs */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Field_T *R1_f  = 0;
    Field_T *R2_f  = 0;
    
    /* surface fields also are used for the interpolation in X_of_x function */
    if (patch->coordsys == CubedSpherical)
    {
      R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f;
      R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
      if (R1_f)
        make_coeffs_2d(R1_f,0,1);/* X and Y direction */
      if (R2_f)
        make_coeffs_2d(R2_f,0,1);/* X and Y direction */
    }
    else if (patch->coordsys == Cartesian)
    {
      R1_f  = R2_f = 0;
    }
    else
      abortEr(NO_OPTION);
  }
  
  OpenMP_1d_Pragma(omp parallel for)
  for (p = 0; p < npoints; ++p)
  {
    Needle_T *needle = alloc_needle();
    double x[3],X[3];
    unsigned pn = 0;
    
    needle->grid = grid;
    FOR_ALL_PATCHES(pn,grid)
    {
      Patch_T *patch = grid->patch[pn];
      needle_in(needle,patch);
    }
    x[0] = pnt->x[p];
    x[1] = pnt->y[p];
    x[2] = pnt->z[p];
    needle->x = x;
    point_finder(needle);
    if (!needle->Nans)
    {
      fprintf(stderr,"(%g,%g,%g) is troublesome!\n",x[0],x[1],x[2]);
      abortEr("It could not find a point!\n");
    }
    else
    {
      pn = needle->ans[0];
      pnt->patchn[p] = pn;
      if(!X_of_x(X,x,grid->patch[pn]))
        abortEr("It could not find X of x!\n");
      pnt->X[p] = X[0];
      pnt->Y[p] = X[1];
      pnt->Z[p] = X[2];
    }
    free_needle(needle);
  }
  
  /* translate fields from BAM notation to Elliptica notation */
  fields_name = translate_fields_name(&bam_fields);
  
  /* set bam fields based on initial data to be usable for bam */
  bbn_bam_set_bam_fields(grid);
  
  /* to avoid race condition between threads write all coeffs */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    unsigned fn = 0;
    
    while(fields_name[fn])
    {
      Field_T *field = patch->pool[Ind(fields_name[fn])];
      make_coeffs_3d(field);
      fn++;
    }
  }
  
  /* open fields_file and start interpolating and writing */
  file = fopen(fields_file_path,"wb");
  pointerEr(file);
  fprintf(file,"# this file contains values of %s\n",bam_fields_name);
  sprintf(title_line,"%s",HEADER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  
  interp_v = alloc_double(npoints);
  f = 0;
  while(fields_name[f])
  {
    printf("~> Interpolating and writing into disk: %s\n",fields_name[f]);
    fflush(stdout);
    /* write it into the fields_file */
    FWriteP_bin(bam_fields[f],strlen(bam_fields[f])+1);
    /* interpolating each fields at the all given points */
    OpenMP_1d_Pragma(omp parallel for)
    for (p = 0; p < npoints; ++p)
    {
      Patch_T *patch  = grid->patch[pnt->patchn[p]];
      Interpolation_T *interp_s = init_interpolation();
      interp_s->field = patch->pool[Ind(fields_name[f])];
      interp_s->XYZ_dir_flag = 1;
      interp_s->X = pnt->X[p];
      interp_s->Y = pnt->Y[p];
      interp_s->Z = pnt->Z[p];
      plan_interpolation(interp_s);
      interp_v[p] = execute_interpolation(interp_s);
      free_interpolation(interp_s);
    }
    
    for (p = 0; p < npoints; ++p)
    {
      /* write it into the fields_file */
      FWriteV_bin(interp_v[p],1);
    }
    f++;
  }
  sprintf(title_line,"%s",FOOTER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  sprintf(msg,"%s",END_MSG);
  FWriteP_bin(p_msg,strlen(msg)+1);
  fclose(file);
  
  _free(interp_v);
  free_2d(fields_name);
  free_2d(bam_fields);
}

/* load grid from the checkpoint file */
static Grid_T *load_grid_from_checkpoint_file(void)
{
  printf("~> Loading grid ...\n");
  fflush(stdout);
  
  Grid_T *grid = 0;
  FILE *file   = 0;
  
  /* open checkpoint file */
  file = fopen(checkpoint_path,"r");
  pointerEr(file);
  grid = bbn_init_from_checkpoint(file);
  fclose(file);
  
  /* extrapolate metric fields inside the BH */
  Pseti("STOP",0);
  bbn_extrapolate_metric_fields_insideBH(grid);
  
  return grid;
}
  
/* read (x,y,x) points from bam file to be interpolated on them,
// and then populate (x,y,z) part of for the interpolation */
static void load_coords_from_coords_file(struct interpolation_points *const pnt)
{
  printf("~> Loading coordinates for interpolation ...\n");
  fflush(stdout);
  
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
    assert(isfinite(pnt->x[i]));
    
    FReadV_bin(pnt->y[i]);
    assert(isfinite(pnt->y[i]));
    
    FReadV_bin(pnt->z[i]);
    assert(isfinite(pnt->z[i]));
  }
  FReadP_bin(match_str)
  if (strcmp(match_str,FOOTER))
    abortEr("It could not find the footer.\n");
  _free(match_str);
  
  fclose(file);
}

/* translating fields name from BAM to Elliptica */
static char **translate_fields_name(char ***const ptr_bam_fields)
{
  char **fields_name = 0;
  char **bam_fields  = 0;
  unsigned nf = 0;
  
  bam_fields = read_separated_items_in_string(bam_fields_name,',');
  nf = 0;
  while (bam_fields[nf])
  {
    ifcmpM("alpha")
    {
      add2fieldsname_0ind_M("bam_alpha");
    }
    elseifcmpM("grhd_rho")
    {
      add2fieldsname_0ind_M("bam_grhd_rho");
    }
    elseifcmpM("grhd_p")
    {
      add2fieldsname_0ind_M("bam_grhd_p");
    }
    elseifcmpM("grhd_epsl")
    {
      add2fieldsname_0ind_M("bam_grhd_epsl");
    }
    elseifcmpM("beta")
    {
      add2fieldsname_1ind_M("beta","bam_Beta","U");
    }
    elseifcmpM("grhd_v")
    {
      add2fieldsname_1ind_M("grhd_v","bam_grhd_v","U");
    }
    elseifcmpM("adm_g")
    {
      add2fieldsname_2ind_M("adm_g","bam_adm_g","D","D");
    }
    elseifcmpM("adm_K")
    {
      add2fieldsname_2ind_M("adm_K","bam_adm_K","D","D");
    }
    else
      abortEr_s("No option has not been defined for %s.\n",bam_fields[nf]);
    printf("~> Translating %s --> %s\n",bam_fields[nf],fields_name[nf]);
    nf++;
  }
  
  (*ptr_bam_fields) = bam_fields;
  return fields_name;
}
