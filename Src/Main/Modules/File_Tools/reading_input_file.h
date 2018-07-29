#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"

#define COMMENT '#'
#define ENTER   '\n'
#define SPACE   ' '
#define TAB     '\t'
#define END      0
#define EQUAL   '='
#define MAX_NUM_CHAR 3000

enum FLOW {LEFT,RIGHT};

void read_input_file(const char *const path);
static void *make_buffer(FILE *const input);
static void populate_parameters(const char *const buff);
void add_parameter(const char *const lv, const char *const rv);
