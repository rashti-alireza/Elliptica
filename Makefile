##################
# Alireza Rashti #
#                #
# April 2020     #
##################

# Using GNU make to compile a software written in C language.
# It uses MyConfig file to get the source files and then 
# detects the dependencies automatically and finally after constructing
# the libraries (shared or static) it makes the exeutable output.
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
#	       |   Cores/  Includes/          |      Includes/
#	       |                              |
#	  +----+-----+-----------+            +-------+------+
#	  |          |           |            |       |      |
#      Module1/   Module2/    ...             P1/     P2/    ...
#         |          |                        |       |
#      Makefile   Makefile                    |       |
#                                             |       |
#                                         Makefile  Makefile
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
EXEC := Elliptica
# exe directory:
EXEC_DIR := $(TOP)/Exe
# library directory, see the above sketch
LIB_DIR := $(TOP)/Lib
# projects dir
PROJECT_DIR := $(TOP)/Src/Projects
# modules dir
MODULE_DIR := $(TOP)/Src/Main/Modules
########################################################################
##############
## C compiler:
##############
CC = gcc
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
INCS += -I$(TOP)/Src/Main/Cores
# special includes
SPECIAL_INCS =
# library path
LIBS = -L$(TOP)/$(LIB_DIR)
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
###########################
## c files and directories:
###########################
# all c file directories
C_DIRS  = $(TOP)/Src/Main/Cores
C_DIRS += $(MODULE)
C_DIRS += $(PROJECT)
C_DIRS := $(strip $(C_DIRS))# strip extra spaces
# all c file paths
C_FILES = $(foreach d,$(C_DIRS),$(wildcard $(d)/*.c))
########################################################################
################################
## object files and directories:
################################
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
# making o_files corresponding to their c files mirror
#o_files := $(C_FILES:%.c=%.o);
#o_files:= \
#	$(foreach f,$(C_FILES),\
#	  $(join\
#	    $(addprefix \
#	      $(LIB_DIR)/$(O_DIR)/, $(notdir $(C_DIRS))\
#	     )/,\
#	    $(notdir $(f:.c=.o))\
#	   )\
#	 )
########################################################################
####################
## summery of flags:
####################
CFLAGS = $(DFLAGS) $(OFLAGS) $(WARN) $(INCS) $(SPECIAL_INCS)
# all linking flags:
LDFLAGS = $(LIBS) $(SPECIAL_LIBS) $(SYSTEM_LIBS)
########################################################################
#####################
## rules and targets:
#####################

export

## default rule to construct EXEC
all: $(EXEC)| $(EXEC_DIR) 
	@true
.PHONY: all

## make the executable out of the object files
$(EXEC): MyConfig | $(LIB_DIR)
	for x in $(C_DIRS); do $(MAKE) -C $$x;echo $$x; done
#	@echo $(PR_F1) $< $(PR_F2)

#%.o : %.c
#$(C_DIRS)/%.o: 
#	@echo $(PR_F1) $(C_DIRS) $(PR_F2)
#	@echo $(PR_F1) $(C_FILES) $(PR_F2)
#	@echo $(PR_F1) $(O_DIRS) $(PR_F2)
#	@echo $(PR_F1) $@ $(PR_F2)
#	#
#	#$(CC) $(CFLAGS) -o $@ -c $<

#%.o:%.c
	
## make object file
#.PHONY: $(compile.o)
#%.c:
#	@echo $(PR_F1) $*.o $(PR_F2)
#	
#	$(CC) $(CFLAGS) -o $@ -c $<
	
#@echo $(PR_F1) %.c $(PR_F2)
#@echo $@
## if EXEC_DIR does not exist make it.	
$(EXEC_DIR):
	@echo $(PR_F0)" mkdir $@"
	@mkdir -p $@

## if LIB_DIR does not exist make it.
$(LIB_DIR): | $(O_DIRS)
	@echo $(PR_F0)" mkdir $@"
	@mkdir -p $@

## if O_DIRS does not exist make it.
$(O_DIRS):
	@echo $(PR_F0)" mkdir $@"
	@mkdir -p $@
	
#.PHONY: compile
#%.o : %.c
#%.o: %.c
#	$(CC) $(CFLAGS) -c $< -o $@
	
#compile: $(obj)
#	$(CC) $(CFLAGS) -o $(EXEC_DIR)/$(EXEC) $(obj) $(LDFLAGS)

#.PHONY: find_interdependency
#find_interdependency: $(depend_paths)
#$(depend_paths):%.depend:%.c MyConfig# should I use myconfig here????
#	@echo $@
#	@set -e; rm -f $@;\
#	$(CC) $(DEPFLAGS) $(CFLAGS) $< > $@.$$$$; \
#	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
#	rm -f $@.$$$$

#.PHONY:print
#print:find_interdependency
#	@echo "starts here:\n"$(depend_paths) | tr " " "\n"

# now include all of the inter-dependency files
#include $(depend_paths)
#$(depend_files):
#include $(wildcard $(depend_paths))

# if there is no MyConfig file, use the prototype
MyConfig:
	-if [[ ! -f MyConfig ]];\
	then \
        cp Doc/MyConfig.example MyConfig; \
        fi


# figure out the dependencies


# make the library

# make the executable

#####################
# make it more general in case if I want to remove a module or project
# deleted


#all:
#	@echo "Lib=" $(Lib)

#######################################################################
################################	
## some variable for nice print:
################################
# new line variable
define PR_NL


endef
# print line with -
define PR_L0
"-------------------------------------------------------------------------"
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
#######################################################################

