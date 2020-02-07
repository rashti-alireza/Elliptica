/* some useful messages */
#define INCOMPLETE_FUNC "Other options have not been developed yet for this part!\n"
#define NO_JOB "No job has been defined for this case."
#define NO_OPTION "No such option has been defined."
#define ERROR_MASSAGE  ":(\n\nERROR and EXIT:\n"

/* some useful functions and commands */
#define UNUSED(x) (void)(x);
#define pointerEr(x)   checkup_pointer_error(x,__FILE__,__LINE__)
#define bad_inputEr()  bad_input_error(__FILE__,__LINE__)
#define null_pathEr(x) null_path_error(x,__FILE__,__LINE__)
#define abortEr(x)     abort_error(x,__FILE__,__LINE__)
#define abortEr_s(x,y) abort_error_string(x,y,__FILE__,__LINE__)
#define Ind(x)  LookUpField_E(x,patch)/* gives error if field is not found */
#define _Ind(x) LookUpField(x,patch)/* gives minus if field is not found */
#define FOR_ALL_PATCHES(n,grid) for ((n) = 0; (n) < grid->np; ++(n))/* loop over all patches of the given grid */
#define FOR_ALL_POINTS(n,patch) for ((n) = 0; (n) < patch->nn; ++(n))/* loop over all points of the given patch */
#define FOR_ALL(x,y) for((x) = 0; y[(x)] != 0; (x)++)
#define TIMER_ON(x) double x = get_time_sec();
#define TIMER_OFF(x) pr_spent_time(x,#x);
#define FOR_SURFACE(x,y,z,n0,n1,n2) (z) = (n2); for ((x) = 0; (x) < (n0); ++(x))\
                                                for ((y) = 0; (y) < (n1); ++(y))
#define FOR_ijk(x,y,z,x_i,x_f,y_i,y_f,z_i,z_f) \
                      for ((x) = (x_i); (x) < (x_f); ++(x))\
                       for ((y) = (y_i); (y) < (y_f); ++(y))\
                        for ((z) = (z_i); (z) < (z_f); ++(z))

#define x_(i)	x_coord((i),patch)
#define y_(i)	y_coord((i),patch)
#define z_(i)	z_coord((i),patch)
#define X_(i)	X_coord((i),patch)
#define Y_(i)	Y_coord((i),patch)
#define Z_(i)	Z_coord((i),patch)

/* variables and fields macors */
#define ADD_AND_ALLOC_FIELD(xNAME)     add_field(#xNAME,0,patch,YES);/* add field to patch->pool and alloc memory */
#define ADD_FIELD(xNAME)               add_field(#xNAME,0,patch,NO);/* add field to patch->pool BUT not alloc memory */
#define REMOVE_FIELD(xNAME)            remove_field(xNAME);/* remove the field utterly */
#define DECLARE_FIELD(xNAME)           Field_T *const xNAME = patch->pool[Ind(#xNAME)];/* access to the whole field */
#define DECLARE_AND_EMPTY_FIELD(xNAME) DECLARE_FIELD(xNAME)/* declare field */\
                                       empty_field(xNAME);/* free v,v2 and info of field */
                                       
/* access to the memory values to modify */
#define WRITE_v(xNAME)  const int _field_index_of_##xNAME = Ind(#xNAME);\
                        free_coeffs(patch->pool[_field_index_of_##xNAME]);\
                        double *const xNAME = patch->pool[_field_index_of_##xNAME]->v;
                        
/* access to the memory values READ ONLY */
#define READ_v(xNAME)   const double *const xNAME = patch->pool[Ind(#xNAME)]->v;

/* it frees v,v2,info of field and alloc memory for v */
#define REALLOC_v_WRITE_v(xNAME)   const int _field_index_of_##xNAME = Ind(#xNAME);\
                                   Field_T *const _F_##xNAME         = patch->pool[_field_index_of_##xNAME];\
                                   empty_field(_F_##xNAME);\
                                   _F_##xNAME->v                     = alloc_double(patch->nn);\
                                   double *const xNAME               = patch->pool[_field_index_of_##xNAME]->v;

/* access to the memory values READ ONLY and unuse to avoid gcc warning */
#define READ_v_UNUSED(xNAME)  READ_v(xNAME)\
                              UNUSED(xNAME);

                                        
/* it compactifies the prepration of Jacobian of derivatives */
#define JACOBIAN_DERIVATIVE(xNAME) const char *types_##xNAME[] = {#xNAME,0};\
                                   prepare_Js_jacobian_eq(patch,types_##xNAME);\
                                   Matrix_T *j_##xNAME = get_j_matrix(patch,#xNAME);\
                                   fJs_T *xNAME        = get_j_reader(j_##xNAME);

/* parameters */                                  
//#define Pcmps(x)
//#define Pseti(x)
//#define Psets(x)
//#define Paddg(x)



/* get value of a string parameter */
#define Pgets(x)   get_parameter_value_S(x,__FILE__,__LINE__,FATAL)/* if not exists give error */
#define PgetsEZ(x) get_parameter_value_S(x,__FILE__,__LINE__,NONE)/* if not exist go easy */

/* get value of an integer parameter */
#define Pgeti(x)   get_parameter_value_I(x,__FILE__,__LINE__,FATAL)/* if not exists give error */
#define PgetiEZ(x) get_parameter_value_I(x,__FILE__,__LINE__,NONE)/* if not exist go easy */

/* get value of a double parameter */
#define Pgetd(x)   get_parameter_value_D(x,__FILE__,__LINE__,FATAL)/* if not exists give error */
#define PgetdEZ(x) get_parameter_value_D(x,__FILE__,__LINE__,NONE)/* if not exist go easy */

/* get value of an array of parameter in double */
#define Pgetdd(x)   get_parameter_array_format(x,__FILE__,__LINE__,FATAL)/* if not exists give error */
#define PgetddEZ(x) get_parameter_array_format(x,__FILE__,__LINE__,NONE)/* if not exist go easy */

#define AddParameterDoubleF(x,y)    update_parameter_double_format(x,y)
#define UpdateParameterDoubleF(x,y) update_parameter_double_format(x,y)

#define PgetdoubleF_E(x) get_parameter_double_format(x,__FILE__,__LINE__,FATAL)/* if not exists give error */
#define PgetdoubleF(x) get_parameter_double_format(x,__FILE__,__LINE__,NONE)/* if not exist go easy */



/* OpenMP for 2 dimension */
#ifdef Pragma_OpenMP_2d
#define OpenMP_2d_Pragma(x) _Pragma ( #x )
#else 
#define OpenMP_2d_Pragma(x)
#endif

/* OpenMP for 1 dimension */
#ifdef Pragma_OpenMP_1d
#define OpenMP_1d_Pragma(x) _Pragma ( #x )
#else 
#define OpenMP_1d_Pragma(x)
#endif

/* OpenMP for go over all patches */
#ifdef Pragma_OpenMP_Patch
#define OpenMP_Patch_Pragma(x) _Pragma ( #x )
#else 
#define OpenMP_Patch_Pragma(x)
#endif

/* some constants */
#define TEST_SUCCESSFUL 0
#define TEST_UNSUCCESSFUL 1
 