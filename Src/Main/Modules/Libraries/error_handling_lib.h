void checkup_pointer_error(void *const p, char *file, int line);
void bad_input_error(char *file, int line);
void null_path_error(void *path,char *file, int line);
void abort_error(char *massage,char *file, int line);

#define pointerEr(x)   checkup_pointer_error(x,__FILE__,__LINE__)
#define bad_inputEr()  bad_input_error(__FILE__,__LINE__)
#define null_pathEr(x) null_path_error(x,__FILE__,__LINE__)
#define abortEr(x)     abort_error(x,__FILE__,__LINE__)
#define ERROR_MASSAGE       "\t!!!!!!ERROR and ABORT!!!!!!\n"
#define ERROR_MASSAGE_ABORT "\t!!!!!!ABORT DUO TO ERROR!!!!!!\n"


