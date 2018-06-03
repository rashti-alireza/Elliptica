#include "core_lib.h"
#include "error_handling_lib.h"
#include "print_lib.h"

#define COMMENT '#'
#define ENTER   '\n'
#define SPACE   ' '
#define END      0
#define EQUAL   '='

extern char *global_path;
extern char *global_inputfile_name;

enum FLOW {LEFT,RIGHT};

void add_parameter(char *lv, char *rv);
char *make_directory(char *path,char *name,Flag_T flg);
static void *make_buffer(FILE *input);
static void populate_parameters(char *const buff);
