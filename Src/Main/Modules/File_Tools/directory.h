#include "core_lib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAX_ARR 400

void make_directory(char **const path,const char *const name,const Flag_T flg);
