#define ERROR_LINE "x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x"
#define pointerEr(x)   checkup_pointer_error(x,__FILE__,__LINE__)
#define bad_inputEr()  bad_input_error(__FILE__,__LINE__)
#define null_pathEr(x) null_path_error(x,__FILE__,__LINE__)
#define abortEr(x)     abort_error(x,__FILE__,__LINE__)
#define abortEr_s(x,y)     abort_error_string(x,y,__FILE__,__LINE__)
#define ERROR_MASSAGE       "\n"ERROR_LINE"\nERROR and ABORT:\n"
#define FOR_ALL(x,y) for(x = 0; y[x] != 0; x++)
#define parameterEr(x) check_parameter(x,__FILE__,__LINE__)
#define FOR_SURFACE(x,y,z,n0,n1,n2) z = n2;for (x = 0; x < n0; ++x)\
                                            for (y = 0,z=n2; y < n1; ++y)
#define FOR_ijk(x,y,z,x_i,x_f,y_i,y_f,z_i,z_f) \
                      for (x = x_i; x < x_f; ++x)\
                       for (y = y_i; y < y_f; ++y)\
                        for (z = z_i; z < z_f; ++z)

                                            