#include "core_lib.h"
#include "error_handling_lib.h"

#define MAX_LENGTH 70
#define L_SYM '>'
#define R_SYM '<'
#define LINE '-'

extern time_t initial_time_global;

void pr_comment(const char *const comment);
void pr_line(void);
void pr_clock(void);
double get_time(void);
void pr_spent_time(const double start,const char *const event);
