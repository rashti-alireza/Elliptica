#ifndef text_and_file_tools_LIB_H
#define text_and_file_tools_LIB_H
#include "elliptica_system_lib.h"

/* this is how we write binary data: first write size and then value. 
// thus, when we wanna read the data the first one gives of the memory allocation 
// and the next gives us value: */

/* write pointer: */
#define FWriteP_bin(x,y) \
if (x){\
  Uint SIZE_ = (Uint)(y);\
  assert(fwrite(&SIZE_,sizeof(SIZE_),1,file));\
  assert(fwrite(x,sizeof(*(x)),SIZE_,file));\
}else{\
  Uint SIZE_ = 0;\
  assert(fwrite(&SIZE_,sizeof(SIZE_),1,file));\
}

/* write variable: */
#define FWriteV_bin(x,y) \
{\
  Uint SIZE_ = (Uint)(y);\
  assert(fwrite(&SIZE_,sizeof(SIZE_),1,file));\
  assert(fwrite(&(x),sizeof(x),SIZE_,file));\
}

/* read pointer */
#define FReadP_bin(x) {\
  Uint SIZE_ = 0;\
  assert(fread(&SIZE_, sizeof(SIZE_),1,file));\
  if (SIZE_) {\
    x = calloc(SIZE_,sizeof(*(x))),IsNull(x);\
    assert(fread(x,sizeof(*(x)),SIZE_,file));}\
  else { x = 0;}}

/* read variable */
#define FReadV_bin(x) {\
  Uint SIZE_ = 0;\
  assert(fread(&SIZE_, sizeof(SIZE_),1,file));\
  assert(fread(&(x),sizeof(x),SIZE_,file));}


int strcmp_i(const char *const s1, const char *const s2);
int strstr_i(const char *const s1, const char *const s2);
char *dup_s(const char *const str);
char *tok_s(char *const str,const char delimit,char **const savestr);
char *make_directory(const char *const path,const char *const name);
char *sub_s(char *const str,const char d1,const char d2,char **const save);
int check_format_s(const char *str,const char *const format);
char *make_folder(const char *const folder);
char *open_folder(const char *const folder);
Uint find_index_string(char **const heystack,const Uint N,const char *const needle);
int regex_search(const char *const regex_pattern,const char *const str);
char *regex_find(const char *const regex_pattern,const char *const str);
char **read_separated_items_in_string(const char *const string,const char delimiter);
void *fopen_and_check(const char *const file_path,const char *const mode,const char *const file_dbg, const int line_dbg);
int regex_replace(const char *const orig/* original */,
                  const char *const regex_pattern/* regex pattern */,
                  const char *const repl/* replace by this piece */,
                  char *const save/* write the result in save  */);


#endif



