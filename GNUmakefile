##################
# Alireza Rashti #
#                #
# April 2020     #
##################

# Using GNU make to compile a software written in C language.
# It uses MyConfig file to get the source files and then 
# detects the dependencies automatically and finally after constructing
# the libraries (shared or static) it makes the executable output.
# A complete documentation can be found:
# https://www.gnu.org/software/make/manual/
#
#
# The hierarchy of the source files and libraries are like the followings;
# Thus, one can use this chart to adjust their project by the variables
# are defined at this Makefile in the section "paths and variables".
#
#
#                                  TOP
#                                   |
#                     +-------------+-------------+
#                     |             |             |
#	             Exe(exe)      Src(*.c,*.h)  Lib(*.o,*.so)
#	                            |
#	              +-------------+-------------+
#	              |                           |
#	            Main/         	       Projects/  
#	              |                           |
#	       +------+------+                +---+-----+
#	       |      |      |                |         |
#	       |    Core/  Includes/          |      Includes/
#              |      |                       |
#              |    main.c,makefile,...       |
#	       |                              |
#	  +----+-----+-----------+            +-------+------+
#	  |          |           |            |       |      |
#      Module1/   Module2/    ...             P1/     P2/    ...
#         |          |                        |       |
#      makefile   makefile                    |       |
#                                             |       |
#                                         makefile  makefile
#
#
#
#
# usage:
# $ make target
#
# example:
# --------
# $ make         # make the default target
# $ make install # it installs the software
# $ make clean   # it cleans the libraries and exe and junks
# $ make -f makefile_name # it uses makefile_name for make
# $ make -k # it runs and ignores the errors
# $ make -n # it only shows the sketch of make and doesn't make anything
# $ make -C dir # means cd to dir and then invoke make
# $ make -j4    # using 4 processors to build the target
#
# NOTE: you can find the master sub-make file for all modules and projects 
# at ./Doc/master_submake, this is the DEFAULT.
#
#
# To report bug or feedback email me at "rashti.alireza@gmail.com" !
########################################################################
#######################
## paths and variables:
#######################
# Top directory. see the above sketch.
TOP :=$(shell pwd)
# if TOP does not exist
ifeq ($(TOP),)
$(error $(PR_NL)"Could not find the top level directory!"$(PR_NL))
endif
# program name exe
EXEC := elliptica
# exe directory:
EXEC_DIR := $(TOP)/Exe
# library directory, see the above sketch
LIB_DIR := $(TOP)/Lib
# projects dir
PROJECT_DIR := $(TOP)/Src/Projects
# modules dir
MODULE_DIR := $(TOP)/Src/Main/Modules
# core directory, where main.c is, see the above sketch
CORE_DIR := Core
# ID reader library
IDR_TOP := $(TOP)/ID_Reader
IDR_LIB_DIR := $(IDR_TOP)/lib
# ID reader library
IDR_INC_DIR := $(TOP)/ID_Reader/include
# ID reader library name
IDR_LIB_NAME := libelliptica_id_reader.a
# master sub-make file path
MASTER_SUB_MAKE_FILE := $(TOP)/Doc/Make/master_submake
# sub-make file name stem in each directory
SUB_MAKE_NAME_STEM := .sub_makefile
# sub-make options:
SUB_MAKE_FLAGS := --no-print-directory
SUB_MAKE_FLAGS += --warn-undefined-variables
########################################################################
##############
## C compiler:
##############
CC     = gcc
# some default flags for the compiler
OFLAGS =# to be determined in MyConfig
WARN   =# to be determined in MyConfig
DFLAGS =# to be determined in MyConfig
# finding c files inter-dependencies using DEPFLAGS flag of the compiler
DEPFLAGS = -M
# ar command to archive the object files
AR = ar
# archive flags
ARFLAGS = rcs
# include path
INCS  = -I$(MODULE_DIR)/Includes
INCS += -I$(PROJECT_DIR)/Includes
INCS += -I$(TOP)/Src/Main/$(CORE_DIR)
# special includes
SPECIAL_INCS =
# libs from compiling of c files of modules and projects
C_LIBS=
# search path for C_LIBS
C_LIBS_PATH = -L$(LIB_DIR)
# special libs
SPECIAL_LIBS =
# system library
SYSTEM_LIBS = -lm
########################################################################
############
## MyConfig:
############
# inlcude MyConfig for more options and c source files
MODULE  =# to be determined in MyConfig
PROJECT =# to be determined in MyConfig
PROJECT_REPO=# to be set in MyConfig
include MyConfig
# MyConfig file path
MyConfig_FILE := $(TOP)/MyConfig
########################################################################
####################
## auto generations:
####################

# name of the auto generated c file
auto_gen_c_file_name := projects_data_base_MADE_BY_MAKE.c
# where we put the auto generated c file
auto_gen_c_file := $(TOP)/Src/Main/$(CORE_DIR)/$(auto_gen_c_file_name)
# project names added in MyConfig, NOTE: the name of the project function
# and the name of the project directory assumed to be the same.
PROJECT_NAMES := $(foreach d,$(PROJECT),$(notdir $d))
PROJECT_NAMES := $(strip $(PROJECT_NAMES))
# the following module(s) or project are mandatory for compilation;
# thus, if they don't exist auto generate a function with their 
# folder name but, the generated function does nothing, only return.
# NOTE: one must adjust $(auto_gen_c_file) target as well.
# NOTE: project returns EXIT_FAILURE.

## modules:
SPECIAL_PR_MODULE_NAMES =
# if the following required module(s) has not been added:
required_mod = $(MODULE_DIR)/Prints/pr_hdf5_silo
# get those module which are not listed in $(MODULE)
SPECIAL_PR_MODULE_NAMES = $(notdir $(filter-out $(MODULE),$(required_mod)))

## required Projects
SPECIAL_PR_PROJECT_NAMES =
# if the following required project(s) has not been added:
required_pro += $(PROJECT_DIR)/Approximate_Killing_Vector
required_pro += $(PROJECT_DIR)/BH_NS_Binary_Initial_Data
required_pro += $(PROJECT_DIR)/NS_NS_Binary_Initial_Data
required_pro += $(PROJECT_DIR)/BH_BH_Binary_Initial_Data
required_pro += $(PROJECT_DIR)/Single_NS_Initial_Data

SPECIAL_PR_PROJECT_NAMES = $(notdir $(filter-out $(PROJECT),$(required_pro)))

########################################################################
################################
## c sources and object sources:
################################
# all c file directories
C_DIRS  = $(TOP)/Src/Main/$(CORE_DIR)
C_DIRS += $(MODULE)
C_DIRS += $(PROJECT)
C_DIRS := $(strip $(C_DIRS))# strip extra spaces
# strip the last slash if any:
C_DIRS := $(foreach d,$(C_DIRS),$(patsubst %/,%,$d))
# all c file paths
C_FILES = $(foreach d,$(C_DIRS),$(wildcard $d/*.c))
# all header files in INCLUDES directories
H_FILES  = $(wildcard $(MODULE_DIR)/Includes/*.h)
H_FILES += $(wildcard $(PROJECT_DIR)/Includes/*.h)
# obj directories:
# the convention is we take the name of each dir in C_DIRS 
# (its last directory name) as the obj directory which contains 
# the object files compiled from all of the c files in that c_dir.
# this is the parent dir where all of O_DIRS located
O_TOP := $(LIB_DIR)/Obj
# extract the last directory name
O_DIRS := $(notdir $(C_DIRS))
# make the full path for the object directories
O_DIRS := $(foreach d,$(O_DIRS),$(O_TOP)/$(d))
# strip extra spaces
O_DIRS := $(strip $(O_DIRS))
# sub make files:
SUB_MAKE_FILES := $(foreach d,$(C_DIRS),$(join $d/,$(SUB_MAKE_NAME_STEM)))
########################################################################
########################################
## recap flags and export all variables:
########################################
# making all C_LIBS strings:
C_LIBS := $(foreach d,$(C_DIRS),$(addprefix -l, $(notdir $d)))
# Note: to resolve inter library dependenciesI added C_LIBS few times.
# you can add more if the linking at the last step fails.
C_LIBS += $(C_LIBS)
C_LIBS += $(C_LIBS)
# compiler flags:
CFLAGS = $(DFLAGS) $(OFLAGS) $(WARN) $(INCS) $(SPECIAL_INCS)
# all linking flags, Note: the order matters!
LDFLAGS = $(C_LIBS_PATH) $(C_LIBS) $(SPECIAL_LIBS) $(SYSTEM_LIBS)
# global crucial dependency files:
DEPENDENCY_FILES = $(H_FILES) $(MyConfig_FILE) 
# exporting to other submakes
export
########################################################################
#####################
## rules and targets:
#####################
## main.o corresponds to main.c
MAIN := $(O_TOP)/$(CORE_DIR)/main.o
##
## default rule to construct EXEC
install: $(EXEC)
	@echo $(PR_L0)
	@echo $(PR_F1) "successful compilation for '$(EXEC)'"
	@echo $(PR_F1) "find '$(EXEC)' at '$(EXEC_DIR)'"
	@echo $(PR_F1) "Thanks!"
	
	@true
.PHONY: install
##
## make the executable out of the object files
$(EXEC): $(DEPENDENCY_FILES) $(auto_gen_c_file) $(SUB_MAKE_FILES) | $(LIB_DIR) $(EXEC_DIR) MyConfig
# --> print
	@echo $(PR_F0) "compiling '$(EXEC)':"
	@echo $(PR_L0)
# --> invoke submakes with the default target
	@for d in $(C_DIRS); \
	  do \
	    $(MAKE) -f $(SUB_MAKE_NAME_STEM) $(SUB_MAKE_FLAGS) -C $$d; \
	  done
# --> link all of the libaries to build the EXEC:
	@$(call cmd_and_pr_func, $(CC) $(CFLAGS) -o $(EXEC_DIR)/$@ $(MAIN) $(LDFLAGS),$(EXEC))
# ------------------------------------------------------------------- #
##
## adding all of the determined projects at Myconfig 
## into a c file in Core to be compiled. 
## Note: this depends on how the automation is desinged for the code.
$(auto_gen_c_file): $(DEPENDENCY_FILES)
# --> if file exists delete it:
	@if [ -f $(auto_gen_c_file) ];\
	 then \
	   rm -rf $(auto_gen_c_file); \
	 fi
# --> headers and declarations:
	@echo "/* this is a generated code by make */" >> $(auto_gen_c_file)
	@echo "#include \"core_lib.h\"" >> $(auto_gen_c_file)
	@echo "void add_project(ProjFunc *const projfunc, const char *const name, const char *const des);" >> $(auto_gen_c_file)
	@echo "int create_db_projects(void);" >> $(auto_gen_c_file)
	@for p in $(PROJECT_NAMES); \
	  do \
	    echo "int " $$p "(void *vp);" >> $(auto_gen_c_file) ;\
	   done;
# --> add function to add the projects:
	@echo "int create_db_projects(void){" >> $(auto_gen_c_file)
	@for p in $(PROJECT_NAMES); \
	  do \
	    echo "  add_project(" "$$p,""\"$$p\",0);" >> $(auto_gen_c_file) ; \
	   done;
	@echo "  return EXIT_SUCCESS;" >> $(auto_gen_c_file)
	@echo "}" >> $(auto_gen_c_file)
# --> add function if spacial project does not exist:
	@for m in $(SPECIAL_PR_PROJECT_NAMES); \
	  do \
	   echo "int " $$m "(void *vp);" >> $(auto_gen_c_file) ; \
	   echo "int " $$m "(void *vp){" >> $(auto_gen_c_file) ; \
	   echo "UNUSED(vp);" >> $(auto_gen_c_file) ; \
	   echo "return EXIT_FAILURE;}" >> $(auto_gen_c_file) ; \
	 done
# --> add function if spacial print module does not exist:
	@for m in $(SPECIAL_PR_MODULE_NAMES); \
	  do \
	   echo "void " $$m "(Pr_Field_T *const pr);" >> $(auto_gen_c_file) ; \
	   echo "void " $$m "(Pr_Field_T *const pr){" >> $(auto_gen_c_file) ; \
	   echo "UNUSED(pr);}" >> $(auto_gen_c_file) ; \
	 done
	 	
##
## if EXEC_DIR does not exist make it.	
$(EXEC_DIR):
	@$(call PR_TASK_relPATH,"mkdir",$@)
	@mkdir -p $@
##
## if LIB_DIR does not exist make it.
$(LIB_DIR): | $(O_DIRS)
	@$(call PR_TASK_relPATH,"mkdir",$@)
	@mkdir -p $@
##
## if O_DIRS does not exist make it.
$(O_DIRS):
	@$(call PR_TASK_relPATH,"mkdir",$@)
	@mkdir -p $@
##
## update sub make files:
%/$(SUB_MAKE_NAME_STEM): $(MASTER_SUB_MAKE_FILE)
	@cp $(MASTER_SUB_MAKE_FILE) $@

##
## if there is no MyConfig file, use the prototype
MyConfig:
	@if [ ! -f MyConfig ];\
	then \
          cp Doc/Make/MyConfig.example MyConfig; \
        fi


##
## git_clone the specified projects
git_clone:
# --> invoke git clone
	@for p in $(PROJECT_REPO); \
	  do \
	    echo $(PR_F0) "git clone '$$p':";\
	    cd $(PROJECT_DIR);\
	    git clone $$p; \
	  done
.PHONY: git_clone

##
## clean Lib, auto generated files and dependency files in submake:
clean:
#	@echo $(PR_F0) "cleaning '$(EXEC)':"
#	@echo $(PR_L0)
	@$(call PR_TASK_relPATH,"rm -rf",$(LIB_DIR))
	@-rm -rf $(LIB_DIR)
	@$(call PR_TASK_relPATH,"rm -rf",$(IDR_TOP))
	@-rm -rf $(IDR_TOP)
#	@$(call PR_TASK_relPATH,"rm -rf",$(EXEC_DIR))
#	@-rm -rf $(EXEC_DIR)
	@$(call PR_TASK_relPATH,"rm -rf",$(auto_gen_c_file))
	@-rm -rf $(auto_gen_c_file)
# --> invoke submakes to clean dependency files:
	@for d in $(C_DIRS); \
	  do \
	    $(MAKE) -f $(SUB_MAKE_NAME_STEM) $(SUB_MAKE_FLAGS) -C $$d $@; \
	  done	 	
.PHONY: clean

id_reader:$(EXEC)
	@mkdir -p $(IDR_LIB_DIR)
	@mkdir -p $(IDR_INC_DIR)
	@cp $(PROJECT_DIR)/Includes/elliptica_id_reader_lib.h $(IDR_INC_DIR)
## --> making lib from all object files and update the lib for each new object file
	@for d in $(O_DIRS); \
	  do \
	  	o=$$(find $$d -type f); \
	  	$(AR) $(ARFLAGS) $(IDR_LIB_DIR)/$(IDR_LIB_NAME) $$o; \
	   done
## --> empty command for print
	@$(call cmd_and_pr_func, , $@)
	@echo $(PR_L0)
	@echo $(PR_F1) "find the library and include directories at '$(IDR_TOP)'"
	@echo $(PR_F1) "Thanks!"
.PHONY: id_reader

#######################################################################
####################################
## some tools for print and compile:
####################################
SED := sed

# new line variable
define PR_NL


endef

# print line with -
define PR_L0
"-------------------------------------------------------------------------"
endef

# print half a line with -
define PR_L1
"------------------------------------"
endef

# print ERROR message
define PR_ERROR
":(\nError:\n"
endef

# print -->
define PR_F0
"~~>"
endef

# print ==>
define PR_F1
"==>"
endef

# print <==
define PR_F2
"<=="
endef

# print many newlines and a line, like a divider,
# its important to be written in one line
define PR_DIVIDER
"\n\n\n\n\n\n"$(PR_L0)
endef

# get a string s1 and a path p1, it makes the p1 relative path and
# then prints --> "s1 p1"
define PR_TASK_relPATH
p=$$(echo $(2) | $(SED) 's,$(TOP),.,g' | $(SED) "s,\(.*\),'\1',g");\
if [ $$? -eq 0 ]; \
then \
echo $(PR_F0) $(1) $$p; \
else  \
echo $(PR_F0) $(1) $(2); \
fi;
endef

# command and print function, print error if happens otherwise succinctly.
# useage:
#
# $(call cmd_and_pr_func, command, print_this)
# $(1) refers to command and $(2) refers to print_this.
define cmd_and_pr_func
$(strip $(1)) 2> $@.COMPILE_ERROR ; \
ret=$$? ; \
if [ $$ret -eq 0 ]; \
then \
printf "%s %s %-45s %s\n" $(PR_F0) "make" $(strip $(2)) "=> [SUCCESSFUL] <="; \
else \
echo $(PR_L0) ; \
echo $(PR_ERROR) ; \
cat $@.COMPILE_ERROR ; \
printf "\n%s %s %-45s %s\n" $(PR_F0) "make" $(strip $(2)) "=> [FAILED] <="; \
echo $(PR_L0) ;\
fi ; \
rm -rf $@.COMPILE_ERROR; \
exit $$ret ;
endef

unexport cmd_and_pr_func
#######################################################################

