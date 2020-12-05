#include "core_lib.h"
#include "error_handling_lib.h"

void checkup_pointer_error(const void *const p, const char *const file, const int line);
void bad_input_error(const char *const file, const int line);
void null_path_error(const void *const path,const char *const file, const int line);
void abort_error(const char *const massage,const char *const file, const int line,const int pr_stdout);
void abort_error_string(const char *const massage,const char *const detail,const char *const file, const int line);
void check_parameter(const Flag_T flg,const char *const file, const int line);
static void pr_line_stderr(const char c);



