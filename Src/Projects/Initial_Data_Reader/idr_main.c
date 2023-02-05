/*
// Alireza Rashti
// February 2023
*/

#include "idr_main.h"

/* exporting initial data for evolution codes */

/* field dictionary */
static const char *const Field_Dictionary[] =
{
"alpha",/* lapse: alpha */
"betax","betay","betaz",/* shift: beta^i */
"adm_gxx","adm_gxy","adm_gxz",/* metric: g_ij */
"adm_gyy","adm_gyz","adm_gzz",/* metric: g_ij */
"adm_Kxx","adm_Kxy","adm_Kxz",/* extrinsic curvature: K_ij */
"adm_Kyy","adm_Kyz","adm_Kzz",/* extrinsic curvature: K_ij */

/* matter part */
"grhd_rho",/* primitive rho */
"grhd_p",/* primitive p */
"grhd_epsl",/* primitive epsilon: total_energy_density = grhd_rho(1+grhd_epsl)*/
"grhd_vx","grhd_vy","grhd_vz",/* primitive v, measured by an Eulerian observer, 
                              // v^i = u^i/(alpha u^0) + beta^i / alpha
                              // where u^{mu}=(u^0,u^i) is the 4-velocity of the fluid */
                         
  0/* --> detemine the last pointer */
};



/* tutorial
// --------

const char *checkpnt_path  = "path/to/elliptica/checkpoint/file"
const int Npnts = 16*16*16; // for all x,y,z coords

// initialize
Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(checkpnt_path);

// the list of fields to be interpolated. should be comma separated
idr->ifields  = "alpha,betax,betay,betaz,adm_gxx,adm_gxy";
idr->npoints  = Npnts;
idr->x_coords = a pointer to double type 1D array of Cartesian x coord values;
idr->y_coords = a pointer to double type 1D array of Cartesian y coord values;
idr->z_coords = a pointer to double type 1D array of Cartesian z coord values;

// set parameter for elliptica
idr->param("bh_filler","on",idr);

// interpolate
elliptica_id_reader_interpolate(idr);

// get interpolated values, e.g.,
int ijk = 0;
for(i,j,k)
{
    double g_xx  = idr->field[idr->indx("adm_gxx")][ijk];
    double betax = idr->field[idr->indx("betax")][ijk];
    double alpha = idr->field[idr->indx("alpha")][ijk];
    
    ijk++;
}

// free
elliptica_id_reader_free(idr);

*/

/* -> create a struct for initial data reader. 
// note: checkpoint file includes all info the system. */
Elliptica_ID_Reader_T *elliptica_id_reader_init (
  const char *const checkpnt/* path/to/elliptica/checkpoint/file */
  )
{
  FUNC_TIC
  
  Uint nf;
  
  // some sanity checks
  if (!checkpnt)
  {
    Error1("No checkpoint file path.")
  }
  
  Elliptica_ID_Reader_T *idr = calloc(1,sizeof(*idr));
  IsNull(idr);
  FILE *file = 0;
  const Parameter_T *par = 0;
  
  // set path
  sprintf(idr->checkpoint_path,"%s",checkpnt);
  
  // set index finder
  idr->indx = find_field_index;
  
  // settig param function
  idr->param = set_param_from_evo;
  
  // read checkpoint file
  file = Fopen(idr->checkpoint_path,"r");
  IsNull(file);
  
  // which ID system
  par = parameter_query_from_checkpoint("Project",file);
  idr->system = dup_s(par->rv);
  
  // alloc field
  for (nf = 0; Field_Dictionary[nf]; nf++);
  idr->nfield = nf;
  idr->field  = calloc(nf,sizeof(*idr->field));
  IsNullf(idr->field);
  
  Fclose(file);

  FUNC_TOC  
  
  return idr;
}

/* set parameter for reader from evo code side.
// e.g., set_param_from_evo("bh_filler", "on",idr); --> bh_filler = on. */
static void set_param_from_evo(
          const char *const lv/* e.g., force_balance */, 
          const char *const rv/* e.g., on */,
          Elliptica_ID_Reader_T *const idr)
{
  // checks
  if(!lv && !rv)
  {
    Error1("A param with empty content.");
  }
  
  Uint nparams = idr->nparams;
  
  idr->params_lv = realloc(idr->params_lv,(nparams+1)*(sizeof(*idr->params_lv)));
  IsNull(idr->params_lv);
  idr->params_lv[nparams] = dup_s(rv);

  idr->params_rv = realloc(idr->params_rv,(nparams+1)*(sizeof(*idr->params_rv)));
  IsNull(idr->params_rv);
  idr->params_rv[nparams] = dup_s(rv);
  
  idr->nparams++;
}

/* given the field name, it returns back the index number of the field.
// if could found, gives error.
// ->return: field index. */
static Uint find_field_index(const char *const fname)
{
  Uint indx = Uint_MAX;
  Uint nf   = 0;
  
  while (Field_Dictionary[nf])
  {
    if (strcmp_i(fname,Field_Dictionary[nf]))
    {
      indx = nf;
      break;
    }
    nf++;
  }
  
  if (indx == Uint_MAX)
  {
    Errors("could not find the field '%s'",fname);
  }
  
  return indx;
}

/* ->return success. interpolate fields for the given coordinates */
int elliptica_id_reader_interpolate(Elliptica_ID_Reader_T *const idr)
{
  FUNC_TIC

  if (strcmp_i(idr->system,"BH_NS_binary_initial_data"))
  {
    bhns_evo_exporting_initial_data(idr);
  }
  else
  {
    Error1(NO_OPTION);
  }
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* ->return success. free everything */
int elliptica_id_reader_free(Elliptica_ID_Reader_T *idr)
{
  FUNC_TIC
  Uint n;
  
  Free(idr->system);
  
  // free fields
  for (n = 0; n < idr->nfield; n++)
  {
    Free(idr->field[n]);
  }
  Free(idr->field);
  
  // free params
  for (n = 0; n < idr->nparams; n++)
  {
    Free(idr->params_lv[n]);
    Free(idr->params_rv[n]);
    
  }
  Free(idr->params_lv);
  Free(idr->params_rv);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

