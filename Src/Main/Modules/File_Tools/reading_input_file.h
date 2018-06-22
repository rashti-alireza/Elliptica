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

enum FLOW {LEFT,RIGHT};

static void *make_buffer(FILE *input);
static void populate_parameters(char *const buff);
