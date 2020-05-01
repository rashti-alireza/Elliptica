#include "bbn_headers.h"
#include "physics_observables_lib.h"
#include <unistd.h>

#define MAX_STR_LEN 400
#define LINE_STR    "-------------------------------------------------------------------------"    
#define PR_PARAMETR_IN_FILE(x) \
        { \
          const double par_temp = PgetdEZ(#x); \
          if (par_temp != DBL_MAX) \
          { \
            fprintf(file,"%-30s = %+0.15f\n",#x,par_temp);\
            if (pr_flg) printf("%-30s = %+0.15f\n",#x,par_temp); \
          } \
        }
        
#define PR_PARAMETR_IN_FILE_s(x) \
        {\
          const char *par_temp = PgetsEZ(#x); \
          if (par_temp) \
          { \
            fprintf(file,"%-30s = %s\n",#x,par_temp);\
            if (pr_flg) printf("%-30s = %s\n",#x,Pgets(#x)); \
          } \
        }
        


void bbn_study_initial_data(Grid_T *const grid);
void bbn_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder);
void bbn_print_residual_norms(Grid_T *const grid,const unsigned iteration, const char *const folder);
void bbn_print_properties(Grid_T *const grid,const unsigned iteration, const char *const folder,const char *const open_file_mode,const int pr_flg);
void bbn_measures(Grid_T *const grid);
static double NS_r_average(Grid_T *const grid);



