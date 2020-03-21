#include "bbn_headers.h"
#include "physics_observables_lib.h"
#include <unistd.h>

#define MAX_STR_LEN 400
#define LINE_STR    "-------------------------------------------------------------------------"    
#define PR_PARAMETR_IN_FILE(x) \
        fprintf(file,"%-11s = %+0.15f\n",#x,Pgetd(#x));\
        if (pr_flg) printf("%-11s = %+0.15f\n",#x,Pgetd(#x));
        


void bbn_study_initial_data(Grid_T *const grid);
void bbn_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder);
void bbn_print_residual_norms(Grid_T *const grid,const unsigned iteration, const char *const folder);
void bbn_print_properties(Grid_T *const grid,const unsigned iteration, const char *const folder,const char *const open_file_mode,const int pr_flg);


