#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"

#define COMMENT '#'
#define BACK_SLASH '\\'
#define ENTER   '\n'
#define SPACE   ' '
#define TAB     '\t'
#define END      0
#define EQUAL   '='
#define MAX_NUM_CHAR (9999)

enum FLOW {e_Left,e_Right};

void read_input_file(const char *const path);
static void *make_buffer(FILE *const input);
static void populate_parameters(const char *const buff);
void add_parameter(const char *const lv, const char *const rv);
