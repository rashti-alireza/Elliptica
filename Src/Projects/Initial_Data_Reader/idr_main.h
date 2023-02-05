#include "core_lib.h"
#include "elliptica_id_reader_lib.h"

Elliptica_ID_Reader_T *elliptica_id_reader_init (
  const char *const checkpnt/* path/to/elliptica/checkpoint/file */);
int elliptica_id_reader_interpolate(Elliptica_ID_Reader_T *const idr);
int elliptica_id_reader_free(Elliptica_ID_Reader_T *idr);
void bhns_evo_exporting_initial_data(Elliptica_ID_Reader_T *const idr)

