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
#                                             |       |
#                                             +       +
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

########################################################################
################################
## path and name configurations:
################################
# Top directory. see the above sketch.
TOP :=$(shell pwd)
# if TOP does not exist
ifeq ($(TOP),)
$(error $(pr_nl)"Could not find the top level directory!"$(pr_nl))
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
module  =# to be determined in MyConfig
project =# to be determined in MyConfig
include MyConfig
########################################################################
###########################
## c files and directories:
###########################
# all c file directories
c_dirs  = $(TOP)/Src/Main/Cores
c_dirs += $(module)
c_dirs += $(project)
c_dirs := $(strip $(c_dirs))# strip extra spaces
# all c file paths
c_files = $(foreach d,$(c_dirs),$(wildcard $(d)/*.c))
########################################################################
################################
## object files and directories:
################################
# obj directories:
# the convention is we take the name of each dir in c_dirs 
# (its last directory name) as the obj directory which contains 
# the object files compiled from all of the c files in that c_dir.
# this is the parent dir where all of o_dirs located
o_parent := Obj
# extract the last directory name
o_dirs := $(notdir $(c_dirs))
# make the full path for the object directories
o_dirs := $(foreach d,$(o_dirs),$(LIB_DIR)/$(o_parent)/$(d))
# strip extra spaces
o_dirs := $(strip $(o_dirs))
# making o_files corresponding to their c files mirror
o_files:= \
	$(foreach f,$(c_files),\
	  $(join\
	    $(addprefix \
	      $(LIB_DIR)/$(o_parent)/, $(notdir $(c_dirs))
	     )/,\
	    $(notdir $(f:.c=.o))\
	   )\
	 )
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
## Default rule executed
all: $(EXEC)| $(EXEC_DIR) 
	@true
.PHONY: all

## make the executable out of the object files
$(EXEC): $(o_files)
	@echo $(pr_f1) $(EXEC) $(pr_f2)

$(o_dirs)/%.o: $(c_dirs)/%.c | $(LIB_DIR)
	@echo $(pr_f1) $(c_dirs)/$*.c $(pr_f2)
	$(CC) $(CFLAGS) -o $@ -c $<
	@touch $@

#%.o:%.c
	
## make object file
#.PHONY: $(compile.o)
#%.c:
#	@echo $(pr_f1) $*.o $(pr_f2)
#	
#	$(CC) $(CFLAGS) -o $@ -c $<
	
#@echo $(pr_f1) %.c $(pr_f2)
#@echo $@
## if EXEC_DIR does not exist make it.	
$(EXEC_DIR):
	@echo $(pr_f0)" mkdir $@"
	@mkdir -p $@

## if LIB_DIR does not exist make it.
$(LIB_DIR): | $(o_dirs)
	@echo $(pr_f0)" mkdir $@"
	@mkdir -p $@

## if o_dirs does not exist make it.
$(o_dirs):
	@echo $(pr_f0)" mkdir $@"
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
	
# some variable for nice print:
# new line variable
define pr_nl


endef
# print line with -
define pr_l0
"-------------------------------------------------------------------------"
endef
# print -->
define pr_f0
"-->"
endef
# print ==>
define pr_f1
"==>"
endef
# print <==
define pr_f2
"<=="
endef

