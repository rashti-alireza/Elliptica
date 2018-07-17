#include "core_lib.h"
#include "error_handling_lib.h"

int strcmp_i(const char *const s1, const char *const s2);
int strstr_i(const char *const s1, const char *const s2);
char *dup_s(const char *const str);
char *tok_s(char *const str,const char delimit,char **const savestr);
