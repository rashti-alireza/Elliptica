#include "core_lib.h"
#include "error_handling_lib.h"

#define COMMENT '#'
#define ENTER   '\n'
#define SPACE   ' '
#define END      0
#define EQUAL   '='

enum FLOW {LEFT,RIGHT};

void add_parameter(char *lv, char *rv);
static void *make_buff(FILE *input);
static void populate_parameters(char *const buff);
