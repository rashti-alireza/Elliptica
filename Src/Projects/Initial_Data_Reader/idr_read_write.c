/*
// Alireza Rashti
// December 2020
*/

/* reading and writing tools for exporting of initial data */

#include "idr_read_write.h"


/* given file path, it reads cartesian coordinate from the file
// then populates pnt struct and finds also X coords. */
void 
  idr_load_Cartesian_coordinates_from_file
    (const char *const coords_file_path,ID_Reader_T *const pnt)
{
  FUNC_TIC
  
  FILE *file = 0;
  Grid_T *const grid = pnt->grid;
  Uint npoints    = 0;
  char *match_str = 0;
  char str[STR_LEN_MAX] = {'\0'};
  char *l = 0;
  Uint i,p;
  
  /* some checks */
  if (!grid)
    Error1("Grid is empty!");
  
  /* open and read coords file */
  file = Fopen(coords_file_path,"r");
  
  /* winding file */
  l = fgets(str,STR_LEN_MAX,file);
  FReadP_bin(match_str)
  if (strcmp(match_str,HEADER))
    Error1("It could not find the header");
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
  printf(Pretty0"number of points to interpolate = %u\n",npoints);
  
  /* reading (x,y,z) */
  for (i = 0; i < npoints; ++i)
  {
    FReadV_bin(pnt->x[i]);
    if(!isfinite(pnt->x[i]))
      Error1("bad coordinate.");
    
    FReadV_bin(pnt->y[i]);
    if(!isfinite(pnt->y[i]))
      Error1("bad coordinate.");
    
    FReadV_bin(pnt->z[i])
    if (!isfinite(pnt->z[i]))
     Error1("bad coordinate.");
  }
  FReadP_bin(match_str)
  if (strcmp(match_str,FOOTER))
    Error1("It could not find the footer.\n");
  Free(match_str);
  
  Fclose(file);
  
  /* populating pnt->(X,Y,Z) and pnt->patchn */
  printf(Pretty0"Preparing points for the interpolation ...\n");
  fflush(stdout);
  OpenMP_1d_Pragma(omp parallel for)
  for (p = 0; p < npoints; ++p)
  {
    Patch_T *patch = 0;
    double x[3],X[3];
    
    x[0] = pnt->x[p];
    x[1] = pnt->y[p];
    x[2] = pnt->z[p];
    patch = x_in_which_patch(x,grid->patch,grid->np);
    if (patch && X_of_x(X,x,patch))
    {
      pnt->X[p]      = X[0];
      pnt->Y[p]      = X[1];
      pnt->Z[p]      = X[2];
      pnt->patchn[p] = patch->pn;
    }
    else
    {
      char errmsg[STR_LEN_MAX] = {'\0'};
      sprintf(errmsg,"It could not find X(%f,%f,%f)!\n",x[0],x[1],x[2]);
      Error1(errmsg);
    }
  }

  UNUSED(l)
  FUNC_TOC
}


/* interpolate the given fields_name at the points and 
// write into fields_file. this file will be reading by an evolution code 
// as the initilization its fields. 
// NOTE: the order of fields_name_str and evo_fields_name_str MUST
// be the same.*/
void 
  idr_interpolate_fields_and_write_to_file
    (FILE *const file,ID_Reader_T *const pnt,
     const char *const fields_name_str/* comma separated */,
     const char *const evo_fields_name_str/* comma separated */)
{
  Grid_T *const grid = pnt->grid;
  const Uint npoints = pnt->npoints;
  char **fields_name = 
    read_separated_items_in_string(fields_name_str,',');
  char **evo_fields   = 
    read_separated_items_in_string(evo_fields_name_str,',');
  double *interp_v = 0;
  Uint count_f;
  Uint p,f;
  
  /* some checks */
  if (!grid)
    Error1("Grid is empty!");
  
  if (!fields_name)
    Error1("No fields given!");
    
  if (!evo_fields)
    Error1("No fields given!");
  
  /* to avoid race condition between threads write all coeffs */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Uint fn = 0;
    
    while(fields_name[fn])
    {
      Field_T *field = patch->fields[Ind(fields_name[fn])];
      make_coeffs_3d(field);
      fn++;
    }
  }
  
  /* set f_index, note: it must be set right before interpolation
  // to make sure all fields are added already. */
  pnt->f_index = calloc(grid->np,sizeof(*pnt->f_index)); 
  IsNull(pnt->f_index);
  /* count f */
  count_f = 0;
  while(fields_name[count_f])
    ++count_f;
  
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch  = grid->patch[p];
    assert(patch->pn == p);
    
    pnt->f_index[p] = calloc(count_f,sizeof(*pnt->f_index[p]));
    IsNull(pnt->f_index[p]);
    
    f = 0;
    while(fields_name[f])
    {
      pnt->f_index[p][f] = Ind(fields_name[f]);
      ++f;
    }
  }
  
  interp_v = alloc_double(npoints);
  f = 0;
  while(fields_name[f])
  {
    printf(Pretty0"Interpolating and writing into disk: %s\n",fields_name[f]);
    fflush(stdout);
    /* write it into the fields_file */
    FWriteP_bin(evo_fields[f],strlen(evo_fields[f])+1);
    /* interpolating each fields at the all given points */
    OpenMP_1d_Pragma(omp parallel for)
    for (p = 0; p < npoints; ++p)
    {
      Patch_T *patch  = grid->patch[pnt->patchn[p]];
      Interpolation_T *interp_s = init_interpolation();
      interp_s->field = patch->fields[pnt->f_index[patch->pn][f]];
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
      /* doc test */
      if (!isfinite(interp_v[p]))
      {
        fprintf(stdout,"%s[%s](%g,%g,%g)|x(%g,%g,%g)|X = %g\n",
                fields_name[f],
                grid->patch[pnt->patchn[p]]->name,
                pnt->x[p],pnt->y[p],pnt->z[p],
                pnt->X[p],pnt->Y[p],pnt->Z[p],interp_v[p]);
        fflush(stdout);
        fprintf(stderr,"%s[%s](%g,%g,%g)|x(%g,%g,%g)|X = %g\n",
                fields_name[f],
                grid->patch[pnt->patchn[p]]->name,
                pnt->x[p],pnt->y[p],pnt->z[p],
                pnt->X[p],pnt->Y[p],pnt->Z[p],interp_v[p]);
        fflush(stderr);
        Error1("Doctest failed!\n");
      }
      /* write it into the fields_file */
      FWriteV_bin(interp_v[p],1);
    }
    f++;
  }
  
  Free(interp_v);
  free_2d(fields_name);
  free_2d(evo_fields);
  free_2d_mem(pnt->f_index,grid->np);
  pnt->f_index = 0;
}

/* interpolate the given fields_name at the points and 
// save the values in idr struct.
// NOTE: the order of fields_name_str and evo_fields_name_str MUST
// be the same.*/
void 
  idr_interpolate_fields_and_save_in_array
    (Elliptica_ID_Reader_T *const idr, ID_Reader_T *const pnt,
     const char *const fields_name_str/* comma separated */,
     const char *const evo_fields_name_str/* comma separated */)
{
  Grid_T *const grid = pnt->grid;
  const Uint npoints = pnt->npoints;
  char **fields_name = 
    read_separated_items_in_string(fields_name_str,',');
  char **evo_fields   = 
    read_separated_items_in_string(evo_fields_name_str,',');
  double *interp_v = 0;
  Uint count_f;
  Uint p,f;
  
  /* some checks */
  if (!grid)
    Error1("Grid is empty!");
  
  if (!fields_name)
    Error1("No fields given!");
  
  if (!evo_fields)
    Error1("No fields given!");

  /* to avoid race condition between threads write all coeffs */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Uint fn = 0;
    
    while(fields_name[fn])
    {
      Field_T *field = patch->fields[Ind(fields_name[fn])];
      make_coeffs_3d(field);
      fn++;
    }
  }
  
  /* set f_index, note: it must be set right before interpolation
  // to make sure all fields are added already. */
  pnt->f_index = calloc(grid->np,sizeof(*pnt->f_index)); 
  IsNull(pnt->f_index);
  /* count f */
  count_f = 0;
  while(fields_name[count_f])
    ++count_f;
  
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch  = grid->patch[p];
    assert(patch->pn == p);
    
    pnt->f_index[p] = calloc(count_f,sizeof(*pnt->f_index[p]));
    IsNull(pnt->f_index[p]);
    
    f = 0;
    while(fields_name[f])
    {
      pnt->f_index[p][f] = Ind(fields_name[f]);
      ++f;
    }
  }
  
  interp_v = alloc_double(npoints);
  f = 0;
  // NOTE: fields_name and evo_fields have the same order.
  while(fields_name[f])
  {
    // alloc mem field for id reader
    idr->field[idr->indx(evo_fields[f])] = alloc_double(npoints);
    
    printf(Pretty0"Interpolating and saving: %s\n",fields_name[f]);
    fflush(stdout);

    /* interpolating each fields at the given points */
    OpenMP_1d_Pragma(omp parallel for)
    for (p = 0; p < npoints; ++p)
    {
      Patch_T *patch  = grid->patch[pnt->patchn[p]];
      Interpolation_T *interp_s = pnt->interp_s[p];
      // each point can take different field
      interp_s->field = patch->fields[pnt->f_index[patch->pn][f]];
      /*
      Interpolation_T *interp_s = init_interpolation();
      interp_s->field = patch->fields[pnt->f_index[patch->pn][f]];
      interp_s->XYZ_dir_flag = 1;
      interp_s->X = pnt->X[p];
      interp_s->Y = pnt->Y[p];
      interp_s->Z = pnt->Z[p];
      plan_interpolation(interp_s);
      */
      interp_v[p] = execute_interpolation(interp_s);
      //free_interpolation(interp_s);
    }
    
    for (p = 0; p < npoints; ++p)
    {
      /* doc test */
      if (!isfinite(interp_v[p]))
      {
        fprintf(stdout,"%s[%s](%g,%g,%g)|x(%g,%g,%g)|X = %g\n",
                fields_name[f],
                grid->patch[pnt->patchn[p]]->name,
                pnt->x[p],pnt->y[p],pnt->z[p],
                pnt->X[p],pnt->Y[p],pnt->Z[p],interp_v[p]);
        fflush(stdout);
        fprintf(stderr,"%s[%s](%g,%g,%g)|x(%g,%g,%g)|X = %g\n",
                fields_name[f],
                grid->patch[pnt->patchn[p]]->name,
                pnt->x[p],pnt->y[p],pnt->z[p],
                pnt->X[p],pnt->Y[p],pnt->Z[p],interp_v[p]);
        fflush(stderr);
        Error1("Doctest failed!\n");
      }

      /* save the field */
      idr->field[idr->indx(evo_fields[f])][p] = interp_v[p];
    }
    
    f++;
  }
  
  Free(interp_v);
  free_2d(fields_name);
  free_2d(evo_fields);
  free_2d_mem(pnt->f_index,grid->np);
  pnt->f_index = 0;
}

/* populate ifield coeffs for a MT safe interpolation.
// NOTE: this function itself is not MT safe. */
void idr_set_ifield_coeffs(Elliptica_ID_Reader_T *const idr)
{
  Grid_T *const grid = idr->grid;
  char **fields_name = idr->id_field_names;// elliptica field names we need
  Uint p;
  
  /* some checks */
  if (!grid)
    Error1("Grid is empty!");
  
  if (!fields_name)
    Error1("No fields given!");

  /* to avoid race condition between threads write all coeffs */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Uint fn = 0;
    
    while(fields_name[fn])
    {
      Field_T *field = patch->fields[Ind(fields_name[fn])];
      make_coeffs_3d(field);
      fn++;
    }
  }
}

/* this is an MT safe interpolatation when use for ID reader routine 
// note: x,y,z are Cartesian coords of the evolution code.
// ex: double gxx = idr->fieldx(idr,"adm_gxx",x,y,z);
// ->return: interpolated value of the field at the given points. */
double idr_interpolate_field_thread_safe(
  Elliptica_ID_Reader_T *const idr, 
  const char *const field_name, const double x,const double y, const double z)
{
  Grid_T *const grid = idr->grid;
  Patch_T *patch = 0;
  double X[3] = {0.};
  double x_ell[3] = {0.}; // x coords in the ID coord. system

  // shift coords and find X
  x_ell[0] = x + idr->id_CM[0];
  x_ell[1] = y + idr->id_CM[1];
  x_ell[2] = z + idr->id_CM[2];
  patch = x_in_which_patch(x_ell,grid->patch,grid->np);
  if (!patch || !X_of_x(X,x_ell,patch))
  {
    char errmsg[STR_LEN_MAX] = {'\0'};
    sprintf(errmsg,"It could not find X(%f,%f,%f)!\n",x_ell[0],x_ell[1],x_ell[2]);
    Error1(errmsg);
  }
  
  // find field
  int f_indx = Ind( idr->id_field_names[idr->indx(field_name)] );
  
  // interpolate and return
  return interpolate_at_this_pnt(patch->fields[f_indx],X);
}

/* given field and coordinates, interpolate the field at the given point
// note: the given point X should be in code coordinates, i.e., XYZ coords.
// ->return: interpolated value of the field at the given points */
static double interpolate_at_this_pnt(Field_T *const field, const double X[3])
{
  Interpolation_T *interp_s = init_interpolation();
  double ret = DBL_MAX;
  
  interp_s->field = field;
  interp_s->XYZ_dir_flag = 1;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  plan_interpolation(interp_s);
  ret = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return ret;
}


/* -> binary file to write
// open new binary file and add appropriate header for export purposes.
// fields_name will be the header of this file, for instance,
// this argument can be evaluation's fields name */
void *
  idr_new_binary_file_to_write
    (const char *const file_path,const char *const fields_name)
{
  FILE *file;
  char title_line[STR_LEN_MAX];
  char *const p_title_line = title_line;/* to avoid GCC warning for FWriteP_bin */
  
  /* open fields_file and start interpolating and writing */
  file = Fopen(file_path,"wb");
  fprintf(file,"# this file contains values of %s\n",fields_name);
  sprintf(title_line,"%s",HEADER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);

  return file;
}

/* close the given ID fields file with appropriate footer */
void idr_close_file(FILE *file)
{
  char title_line[STR_LEN_MAX];
  char *const p_title_line = title_line;/* to avoid GCC warning for FWriteP_bin */
  char msg[STR_LEN_MAX];
  char *const p_msg = msg;/* to avoid GCC warning for FWriteP_bin */
  sprintf(title_line,"%s",FOOTER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  sprintf(msg,"%s",END_MSG);
  FWriteP_bin(p_msg,strlen(msg)+1);
  Fclose(file);
}

/* alloc struct ID_Reader_T */
ID_Reader_T *idr_init(void)
{
  ID_Reader_T *pnt = calloc(1,sizeof(*pnt));
  IsNull(pnt);
  return pnt;
}

/* free struct ID_Reader_T */
void idr_free(ID_Reader_T *pnt)
{
  if (!pnt)
    return;
    
  Free(pnt->x);
  Free(pnt->y);
  Free(pnt->z);
  Free(pnt->X);
  Free(pnt->Y);
  Free(pnt->Z);
  Free(pnt->patchn);
  free_2d_mem(pnt->f_index,pnt->grid->np);
  pnt->f_index = 0;
	
	if (pnt->interp_s)
	{
		for (Uint p = 0; p < pnt->npoints; ++p)
		{
			free_interpolation(pnt->interp_s[p]);
		}
		Free(pnt->interp_s);
	}
  Free(pnt);
}

/* get (x,y,z) points from id_reader struct and find the corresponding (X,Y,Z) coords. 
// NOTE: no allocation done for (x,y,z), i.e., we assuemd idr struct already contains 
// (x,y,z) coords.
// NOTE: CM denotes shifting of GIVEN (x,y,z) coords as follows:
// (x,y,z) = (x,y,z)' - CM => (x,y,z) + CM = (x,y,z)', prime denotes elliptica coords. */
void idr_find_XYZ_from_xyz(Elliptica_ID_Reader_T *const idr, 
                           ID_Reader_T *const pnt,const double *const CM)
{
  FUNC_TIC

  Grid_T *const grid = pnt->grid;
  Uint npoints = 0;
  
   /* some checks */
  if (!grid)
    Error1("Grid is empty!");
  
  pnt->x       = idr->x_coords;
  pnt->y       = idr->y_coords;
  pnt->z       = idr->z_coords;
  pnt->npoints = npoints = (Uint)idr->npoints;
  pnt->X       = alloc_double(npoints);
  pnt->Y       = alloc_double(npoints);
  pnt->Z       = alloc_double(npoints);
  pnt->patchn  = calloc(npoints,sizeof(*pnt->patchn));
  IsNull(pnt->patchn);
  pnt->interp_s = calloc(npoints,sizeof(*pnt->interp_s));
  IsNull(pnt->interp_s);
  printf(Pretty0"number of points to interpolate = %u\n",npoints);

  /* populating pnt->(X,Y,Z) and pnt->patchn */
  printf(Pretty0"Preparing points for the interpolation ...\n");
  fflush(stdout);
  OpenMP_1d_Pragma(omp parallel for)
  for (Uint p = 0; p < npoints; ++p)
  {
    Patch_T *patch = 0;
    double x[3],X[3];
    
    x[0] = pnt->x[p] + CM[0];
    x[1] = pnt->y[p] + CM[1];
    x[2] = pnt->z[p] + CM[2];
    patch = x_in_which_patch(x,grid->patch,grid->np);
    if (patch && X_of_x(X,x,patch))
    {
      pnt->X[p]      = X[0];
      pnt->Y[p]      = X[1];
      pnt->Z[p]      = X[2];
      pnt->patchn[p] = patch->pn;
      // setup interpolation scheme
      pnt->interp_s[p] = init_interpolation();
      // choose the very first field of this patch, NOT ideal!
      assert(patch->fields);
      assert(patch->fields[0]);
      pnt->interp_s[p]->field = patch->fields[0];
      pnt->interp_s[p]->XYZ_dir_flag = 1;
      pnt->interp_s[p]->X = pnt->X[p];
      pnt->interp_s[p]->Y = pnt->Y[p];
      pnt->interp_s[p]->Z = pnt->Z[p];
      plan_interpolation(pnt->interp_s[p]);
    }
    else
    {
      char errmsg[STR_LEN_MAX] = {'\0'};
      sprintf(errmsg,"It could not find X(%f,%f,%f)!\n",x[0],x[1],x[2]);
      Error1(errmsg);
    }
  }

  FUNC_TOC
}
