#define pointerEr(x)   checkup_pointer_error(x,__FILE__,__LINE__)
#define bad_inputEr()  bad_input_error(__FILE__,__LINE__)
#define null_pathEr(x) null_path_error(x,__FILE__,__LINE__)
#define abortEr(x)     abort_error(x,__FILE__,__LINE__)
#define abortEr_s(x,y)     abort_error_string(x,y,__FILE__,__LINE__)
#define ERROR_MASSAGE       "\n\nERROR and ABORT:\n"
#define for_all_patches_macro(x,y) for(x = 0; y->patch[x] != 0; x++)
#define parameterEr(x) check_parameter(x,__FILE__,__LINE__)

