/*
// Alireza Rashti
// February 2023
*/

/* exporting initial data for evolution codes */

/* tutorial
// --------

const char *checkpnt_path  = "path/to/elliptica/checkpoint/file"
const int Npnts = 16*16*16; // for all x,y,z coords

// initialize
Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(checkpnt_path);

idr->fields   = "alpha,beta^x,beta^y,beta^z,g_xx,g_xy"; // NOTE: order matters
idr->Npoints  = Npnts;
idr->x_coords = a pointer to double type 1D array of Cartesian x coord values;
idr->y_coords = a pointer to double type 1D array of Cartesian y coord values;
idr->z_coords = a pointer to double type 1D array of Cartesian z coord values;

// interpolate
elliptica_id_reader_interpolate(idr);

// get interpolated values, e.g.,
int ijk = 0;
for(i,j,k)
{
    double g_xx   = idr->field[idr->indx("g_xx")][ijk];
    double beta_x = idr->field[idr->indx("beta^x")][ijk];
    double alpha  = idr->field[idr->indx("g_xx")][ijk];
    
    ijk++;
}

// free
elliptica_id_reader_free(idr);

*/

/* -> create a struct for initial data reader. 
// note: checkpoint file includes all info the system. */
Elliptica_ID_Reader_T * elliptica_id_reader_init (
  const char *const checkpnt/* path/to/elliptica/checkpoint/file */
  )
{
  // some sanity checks
  if (!checkpnt)
  {
    Error1("No checkpoint file path.")
  }
  
  Elliptica_ID_Reader_T *idr = calloc(1,sizeof(*idr));
  IsNull(idr);
  
  sprintf(idr->checkpoint_path,"%s",checkpnt);
  
  // set id system???
  
  return idr;
}



/* ->return success. interpolate fields for the given coordinates */
int elliptica_id_reader_interpolate(Elliptica_ID_Reader_T *const idr)
{
  // read checkpoint file according

  // find X,Y,Z coords for the given x,y,z coords (considering CM translation)
  
  // populate the field of interest
  
  return EXIT_SUCCESS;
}

/* ->return success. free everything */
int elliptica_id_reader_free(Elliptica_ID_Reader_T *idr)
{
  return EXIT_SUCCESS;
}

