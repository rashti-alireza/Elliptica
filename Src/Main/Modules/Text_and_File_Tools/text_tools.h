#include "core_lib.h"
#include "error_handling_lib.h"
#include <regex.h>

int strcmp_i(const char *const s1, const char *const s2);
int strstr_i(const char *const s1, const char *const s2);
char *dup_s(const char *const str);
char *tok_s(char *const str,const char delimit,char **const savestr);
char *sub_s(char *const str,const char d1,const char d2,char **const save);
int check_format_s(const char *str,const char *const format);
Uint find_index_string(char **const heystack,const Uint N,const char *const needle);
int regex_search(const char *const regex_pattern,const char *const str);
char *regex_find(const char *const regex_pattern,const char *const str);
char **read_separated_items_in_string(const char *const string,const char delimiter);
int regex_replace(const char *const orig/* original */,
                  const char *const regex_pattern/* regex pattern */,
                  const char *const repl/* replace by this piece */,
                  char *const save/* write the result in save  */);
                


