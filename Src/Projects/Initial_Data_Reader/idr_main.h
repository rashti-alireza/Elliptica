#include "core_lib.h"
#include "id_reader_lib.h"
#include "elliptica_id_reader_lib.h"
#include "checkpoint_lib.h"

// add ID exporters for new projects here:
void bhns_export_id_generic(Elliptica_ID_Reader_T *const idr);
//void nsns_export_id_generic(Elliptica_ID_Reader_T *const idr);


int Initial_Data_Reader(void *vp);
Elliptica_ID_Reader_T *elliptica_id_reader_init (
  const char *const checkpnt/* path/to/elliptica/checkpoint/file */,
  const char *const option/* the option for exporter */
  );

int elliptica_id_reader_interpolate(Elliptica_ID_Reader_T *const idr);
int elliptica_id_reader_free(Elliptica_ID_Reader_T *idr);
int init_global_variables(const char *const path);
void free_parameter_db(void);
void free_grid_db(void);

// locals
static void set_param_from_evo(
          const char *const lv/* e.g., force_balance */, 
          const char *const rv/* e.g., on */,
          Elliptica_ID_Reader_T *const idr);

static double get_param_double_from_checkpoint( 
          const char *const lv/* e.g., "BHNS_angular_velocity" */,
          Elliptica_ID_Reader_T *const idr);

static Uint find_field_index(const char *const fname);



