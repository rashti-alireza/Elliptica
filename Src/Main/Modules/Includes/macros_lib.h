#ifndef macros_LIB_H
#define macros_LIB_H
#include "elliptica_system_lib.h"

/* str length */
#define MACRO__STR__LEN1 (99)

/* field stem regex format for using with stems */
#define FIELD__STEM__REGEX__FORMAT  "^.+_(U|D)[012].*$"

/* some useful messages and prints*/
#define INCOMPLETE_FUNC "Other options have not been developed yet for this part!\n"
#define NO_JOB    "No job has been defined for this case."
#define NO_OPTION "No such option has been defined. Implement me!"
#define ERROR_MASSAGE    ":(\n\nERROR and EXIT:\n"
#define WARNING_MASSAGE  "WARNING! WARNING! WARNING!\n\n"

#define Pretty0        "|--> "
#define Pretty1        "~> "
#define Pretty2        "--> "

/* some general macros */
#define FOLDER_AFFIX "_%02d"/* that's how folders are named */

/* some useful functions and commands */
#define UNUSED(x) (void)(x);
#define IsNull(x)   checkup_pointer_error(x,__FILE__,__LINE__)
#define Fopen(x,y) fopen_and_check(x,y,__FILE__,__LINE__)
#define Fclose(x)  (x ? fclose(x),(x) = NULL : NULL)
#define bad_inputEr()  bad_input_error(__FILE__,__LINE__)
#define null_pathEr(x) null_path_error(x,__FILE__,__LINE__)
#define Warning(x)   printf(WARNING_MASSAGE"%s\nFile = %s\nLine = %d\n",x,__FILE__,__LINE__)
#define Error0(x)    abort_error(x,__FILE__,__LINE__,0)
#define Error1(x)    abort_error(x,__FILE__,__LINE__,1)
#define Errors(x,y)  abort_error_string(x,y,__FILE__,__LINE__)
#define Ind(x)  LookUpField_E(x,patch)/* gives error if field is not found */
#define _Ind(x) LookUpField(x,patch)/* gives minus if field is not found */
#define FOR_ALL_PATCHES(n,grid) for (n = 0; (n) < (grid)->np; ++(n))/* loop over all patches of the given grid */
#define FOR_ALL_POINTS(n,patch) for ((n) = 0; (n) < (patch)->nn; ++(n))/* loop over all points of the given patch */
#define FOR_ALL(x,y) for((x) = 0; y[(x)] != 0; (x)++)
#define FOR_ALL_ijk   for (Uint ijk = 0; ijk < patch->nn; ++ijk)/* define ijk and loop */
#define FOR_ALL_p(N)  for (Uint p = 0; p < (N); ++p)/* define p and loop */
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

/* function is called  */
#define FUNC_TIC header_and_clock(__func__);

/* function ends  */
#define FUNC_TOC footer_and_clock(__func__);

/* relative (x,y,z) and r respect to center of patch */
#define DEF_RELATIVE_x  double x=patch->node[ijk]->x[0]-patch->c[0];
#define DEF_RELATIVE_y  double y=patch->node[ijk]->x[1]-patch->c[1];
#define DEF_RELATIVE_z  double z=patch->node[ijk]->x[2]-patch->c[2];
#define DEF_RELATIVE_r  double r=sqrt(Pow2(x)+Pow2(y)+Pow2(z));

/* variables and fields macors */
#define ADD_AND_ALLOC_FIELD(xNAME)     add_field(#xNAME,0,patch,YES);/* add field to patch->fields and alloc memory */
#define ADD_FIELD(xNAME)               add_field(#xNAME,0,patch,NO);/* add field to patch->fields BUT not alloc memory */
#define REMOVE_FIELD(xNAME)            remove_field(xNAME);/* remove the field utterly */
#define DECLARE_FIELD(xNAME)           Field_T *const xNAME = patch->fields[Ind(#xNAME)];/* access to the whole field */
#define DECLARE_AND_EMPTY_FIELD(xNAME) DECLARE_FIELD(xNAME)/* declare field */\
                                       empty_field(xNAME);/* free v,v2,v3 and info of field */
                                       
/* access to the memory values to modify */
#define WRITE_v(xNAME)  const int _field_index_of_##xNAME = Ind(#xNAME);\
                        free_coeffs(patch->fields[_field_index_of_##xNAME]);\
                        double *const xNAME = patch->fields[_field_index_of_##xNAME]->v;

/* access to field->v with specified stem => so the indices are adjusted. 
// note: stem is a pointer to char.
// in effect, it takes the indices (if any) of xNAME and trails to stem.
// ex: xNAME = g_D0 and stem = "K" => read field with name "K_D0" */
#define WRITE_v_STEM(xNAME,stem) \
 char field__name__##xNAME[MACRO__STR__LEN1] = {'\0'};\
 if (regex_search(FIELD__STEM__REGEX__FORMAT,#xNAME))\
     {\
     const char *field__index__##xNAME = strrchr(#xNAME,'_');\
     field__index__##xNAME = (field__index__##xNAME ? field__index__##xNAME : "???");\
     sprintf(field__name__##xNAME,"%s%s",stem,field__index__##xNAME);\
     }\
 else{sprintf(field__name__##xNAME,"%s"  ,stem);}\
 const int _field_index_of_##xNAME = Ind(field__name__##xNAME);\
 free_coeffs(patch->fields[_field_index_of_##xNAME]);\
 double *const xNAME = patch->fields[_field_index_of_##xNAME]->v;

                        
/* access to the memory values READ ONLY */
#define READ_v(xNAME)   const double *const xNAME = patch->fields[Ind(#xNAME)]->v;

/* read only of field->v with specified stem => so the indices are adjusted. 
// note: stem is a pointer to char.
// in effect, it takes the indices (if any) of xNAME and trails to stem.
// ex: xNAME = g_D0 and stem = "K" => read field with name "K_D0" */
#define READ_v_STEM(xNAME,stem) \
 char field__name__##xNAME[MACRO__STR__LEN1] = {'\0'};\
 if (regex_search(FIELD__STEM__REGEX__FORMAT,#xNAME))\
     {\
     const char *field__index__##xNAME = strrchr(#xNAME,'_');\
     field__index__##xNAME = (field__index__##xNAME ? field__index__##xNAME : "???");\
     sprintf(field__name__##xNAME,"%s%s",stem,field__index__##xNAME);\
     }\
 else {sprintf(field__name__##xNAME,"%s"  ,stem);}\
 const double *const xNAME = patch->fields[Ind(field__name__##xNAME)]->v;

/* empty_field and alloc memory for v with specified stem => so the indices are adjusted.
// note: stem is a pointer to char. */
#define REALLOC_v_WRITE_v_STEM(xNAME,stem) \
 char field__name__##xNAME[MACRO__STR__LEN1] = {'\0'};\
 if (regex_search(FIELD__STEM__REGEX__FORMAT,#xNAME))\
     {\
     const char *field__index__##xNAME = strrchr(#xNAME,'_');\
     field__index__##xNAME = (field__index__##xNAME ? field__index__##xNAME : "???");\
     sprintf(field__name__##xNAME,"%s%s",stem,field__index__##xNAME);\
     }\
 else {sprintf(field__name__##xNAME,"%s"  ,stem);}\
 const int _field_index_of_##xNAME = Ind(field__name__##xNAME);\
 empty_field(patch->fields[_field_index_of_##xNAME]);\
 patch->fields[_field_index_of_##xNAME]->v = alloc_double(patch->nn);\
 double *const xNAME = patch->fields[_field_index_of_##xNAME]->v;

/* empty_field and alloc memory for v */
#define REALLOC_v_WRITE_v(xNAME)   const int _field_index_of_##xNAME = Ind(#xNAME);\
                                   Field_T *const _F_##xNAME         = patch->fields[_field_index_of_##xNAME];\
                                   empty_field(_F_##xNAME);\
                                   _F_##xNAME->v                     = alloc_double(patch->nn);\
                                   double *const xNAME               = patch->fields[_field_index_of_##xNAME]->v;

/* access to the memory values READ ONLY and unuse to avoid gcc warning */
#define READ_v_UNUSED(xNAME)  READ_v(xNAME)\
                              UNUSED(xNAME);

/* take partial derivatives, 
// NOTE: no semicolon at the end to be more flexible */
#define dField_di(xNAME) partial_derivative(patch->fields[Ind(#xNAME)])

/* take partial derivatives with arbitrary stem
// note: stem is a point to char. */
#define dField_di_STEM(xNAME,stem) \
 char field__name__##xNAME[MACRO__STR__LEN1] = {'\0'};\
 const char *const field__index__##xNAME   = strrchr(#xNAME,'_');\
 sprintf(field__name__##xNAME,"%s%s",stem,field__index__##xNAME);\
 partial_derivative(patch->fields[Ind(field__name__##xNAME)]);

/* add a general parameter */
#define Paddg(x,y) add_parameter(x,y)

/* comparing value of the parameter x with y using strcmp_i */
#define Pcmps(x,y)   strcmp_i(Pgets(x),y)
#define Pcmpss(x,y)  strstr_i(Pgets(x),y)
#define PcmpsEZ(x,y) strcmp_i(PgetsEZ(x),y)

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

/* add or update a double type parameter */
#define Psetd(x,y)  update_parameter_double_format(x,(y),1)
#define Psetd_0print(x,y)  update_parameter_double_format(x,(y),0)

/* add or update an integer type parameter */
#define Pseti(x,y)  update_parameter_integer(x,(y))

/* add or update a string type parameter */
#define Psets(x,y)  update_parameter_string(x,y)

/* set default parameter */
#define Pset_default(x,y) set_default_parameter(x,y)

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

/* relax fashion update of a field by a function 
// func_updator = the function updates the field (with its arguments)
// patch        = computing patch
// fld          = field name (not in string format)
// w            = update weight (<=1.).
// NOTE: it frees field->v and assume func_updator allocate memory. */
#define RELAX_UPDATE_FUNC(func_updator,patch,fld,w) \
{\
 double w2_##fld = 1.-(w);/* weight */\
 Uint ijk__##fld;/* dummy index */\
 /* find field and take care of values */\
 Field_T *const new_field__##fld = patch->fields[Ind(#fld)];\
 free_coeffs(new_field__##fld);\
 /* if we can relax update */\
 if (new_field__##fld->v)\
 {\
  double *const old_value__##fld = new_field__##fld->v;\
  new_field__##fld->v = 0;\
  /* update */\
  func_updator;\
  for ((ijk__##fld) = 0; (ijk__##fld) < patch->nn; ++(ijk__##fld))\
  {\
    new_field__##fld->v[(ijk__##fld)] = \
      (w)*new_field__##fld->v[(ijk__##fld)]+\
      (w2_##fld)*old_value__##fld[(ijk__##fld)];\
  }\
  free(old_value__##fld);\
 }\
 else {func_updator;}\
}

/* free and set it to NULL. */
#define Free(ptr) {if(ptr){free(ptr);ptr = NULL;}}

#endif









