#include "core_lib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "error_handling_lib.h"

#define MAX_ARR 400

char *make_directory(const char *const path,const char *const name);
char *make_folder(const char *const folder);
char *open_folder(const char *const folder);
