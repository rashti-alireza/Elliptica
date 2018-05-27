
void checkup_pointer(void *const p, char *file, int line);
void bad_input_error(char *file, int line);

#define checkup(x) checkup_pointer(x,__FILE__,__LINE__)
#define bad_input() bad_input_error(__FILE__,__LINE__)