#include "core_lib.h"
#include "error_handling_lib.h"
#include "print_lib.h"
#include "utilites.h"

#define COMMENT '#'
#define ENTER   '\n'
#define SPACE   ' '
#define TAB     '\t'
#define END      0
#define EQUAL   '='
#define MAX_NUM_CHAR 3000

extern char *path_global;
extern char *inputfile_name_global;

enum FLOW {LEFT,RIGHT};

void add_parameter(char *lv, char *rv);
char *make_directory(char *path,char *name,Flag_T flg);
static void *make_buffer(FILE *input);
static void populate_parameters(char *const buff);
