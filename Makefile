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
# The hierarchy of the source files and libraries are like the followings:
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
#
#
# NOTE: in the following the capital words are exported and 
# small words are local.
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
########################################################################
##############
## C compiler:
##############
CC     = gcc
# some default flags for the compiler
OFLAGS = -g
WARN   = -Wall
DFLAGS =
# finding c files inter-dependencies using DEPFLAGS flag of the compiler
DEPFLAGS = -M
# ar command to archive the object files
AR = ar
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
include MyConfig
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

# exporting to other submake
export
########################################################################
#####################
## rules and targets:
#####################
## main.o corresponds to main.c
MAIN := $(O_TOP)/$(CORE_DIR)/main.o
##
## default rule to construct EXEC
all: $(EXEC)
	@echo $(PR_L0)"\n\n"
	@echo $(PR_F1) "successful compilation for '$(EXEC)'"
	@echo $(PR_F1) "find '$(EXEC)' at '$(EXEC_DIR)'""\n\n"
	@true
.PHONY: all
##
## make the executable out of the object files
$(EXEC): MyConfig $(H_FILES) | $(LIB_DIR) $(EXEC_DIR)
	@echo $(PR_F0) "build '$(EXEC)':"
	@for d in $(C_DIRS); do $(MAKE) -C $$d; done
	$(CC) $(CFLAGS) -o $(EXEC_DIR)/$@ $(MAIN) $(LDFLAGS)
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
## if there is no MyConfig file, use the prototype
MyConfig:
	-if [ ! -f MyConfig ];\
	then \
        cp Doc/MyConfig.example MyConfig; \
        fi
#######################################################################
#############################
## some tools for nice print:
#############################
SED := sed
# new line variable
define PR_NL


endef
# print line with -
define PR_L0
"--------------------------------------------------------------------------"
endef
# print -->
define PR_F0
"-->"
endef
# print ==>
define PR_F1
"==>"
endef
# print <==
define PR_F2
"<=="
endef
# get a string s1 and a path p1, it makes the p1 relative path and
# then prints "s1 p1"
define PR_TASK_relPATH
p=$$(echo $(2) | $(SED) 's,$(TOP),.,g' | $(SED) "s,\(.*\),'\1',g");\
if [ $$? -eq 0 ]; \
then \
echo $(PR_F0) $(1) $$p; \
else  \
echo $(PR_F0) $(1) $(2); \
fi
endef
#######################################################################

