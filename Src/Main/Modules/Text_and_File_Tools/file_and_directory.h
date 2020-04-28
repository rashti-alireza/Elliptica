#include "core_lib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "error_handling_lib.h"
#include "utilities_lib.h"

#define MAX_ARR  400
#define MAX_ARR2 800

char *make_directory(const char *const path,const char *const name);
char *make_folder(const char *const folder);
char *open_folder(const char *const folder);
void *fopen_and_check(const char *const file_path,const char *const mode,const char *const file_dbg, const int line_dbg);
