#include "bbn_headers.h"

#define Big_value (100)
#define STR_LEN_MAX (1000)
#define HEADER "#{data#"
#define FOOTER "#}data#"
#define END_MSG "\n#file_completed#\n"
/* the followings are some handy macro to deal with strings translation */
#define ifcmpM(x)     if(strstr(bam_fields[nf],x))

#define elseifcmpM(x) else if(strstr(bam_fields[nf],x))

#define str2index(s,i) {\
  if(s[0] == 'x')      i = 0;\
  else if(s[0] == 'y') i = 1;\
  else if(s[0] == 'z') i = 2;\
  else Error0(NO_OPTION);}

/* add x to fields_name */   
#define add2fieldsname_0ind_M(x) {\
  fields_name = realloc(fields_name,(nf+2)*sizeof(*fields_name));\
  IsNull(fields_name);\
  fields_name[nf+1] = 0;\
  fields_name[nf] = dup_s(x);}

/* add stem_Ui or stem_Di to fields_name */
#define add2fieldsname_1ind_M(needle,stem,UD) {\
  fields_name = realloc(fields_name,(nf+2)*sizeof(*fields_name));\
  IsNull(fields_name);\
  fields_name[nf+1] = 0;\
  char fname_[STR_LEN_MAX] = {'\0'};\
  unsigned index_ = 0;\
  unsigned len_ = (unsigned)strlen(bam_fields[nf]);assert(len_);\
  const char *aux_ = &bam_fields[nf][len_-1];/* Fx => aux->x */\
  str2index(aux_,index_);\
  sprintf(fname_,"%s_%s%u",stem,UD,index_);\
  fields_name[nf] = dup_s(fname_);}

/* add stem_UUi or stem_DDi etc to fields_name */
#define add2fieldsname_2ind_M(needle,stem,UD0,UD1) {\
  fields_name = realloc(fields_name,(nf+2)*sizeof(*fields_name));\
  IsNull(fields_name);\
  fields_name[nf+1] = 0;\
  char fname_[STR_LEN_MAX] = {'\0'};\
  unsigned index0_ = 0,index1_ = 0;\
  unsigned len_    = (unsigned)strlen(bam_fields[nf]);assert(len_>1);\
  const char *aux_ = &bam_fields[nf][len_-2];/* Fxy => aux->x */\
  str2index(aux_,index0_);\
  aux_++;\
  str2index(aux_,index1_);\
  sprintf(fname_,"%s_%s%u%s%u",stem,UD0,index0_,UD1,index1_);\
  fields_name[nf] = dup_s(fname_);}

/* strcut for point where interpolate taken place */
struct interpolation_points
{
  double *x,*y,*z;/* (x,y,z) coords */
  double *X,*Y,*Z;/* (X,Y,Z) coords */
  unsigned *patchn;/* patch number for each coord */
  unsigned npoints;/* number of coords */
};

void bbn_bam_export_id(void);
static void load_coords_from_coords_file(struct interpolation_points *const pnt);
static Grid_T *load_grid_from_checkpoint_file(void);
static void interpolate_and_write(Grid_T *const grid,struct interpolation_points *const pnt);
static char **translate_fields_name(char ***const bam_fields);
void bbn_bam_error(const char *const msg,const char *const file,const int line);

static void 
bam_output_doctest
  (
  Grid_T *const grid/* loaded grid */,
  struct interpolation_points *const pnt/* where interpolation takes place */
  );




