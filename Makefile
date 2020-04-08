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

# Top directory. see the above sketch.
TOP :=$(shell pwd)
# if TOP does not exist
ifeq ($(TOP),)
$(error $(n)"Could not find the top level directory!"$(n))
endif

# program name
EXEC := Elliptica

# exe directory:
EXEC_DIR := $(TOP)/Exe

# projects dir
PROJECT_DIR := $(TOP)/Src/Projects

# modules dir
MODULE_DIR := $(TOP)/Src/Main/Modules

# C compiler
CC = gcc

# some default flags for the compiler
OFLAGS = -g
WARN   = -Wall
DFLAGS =

# ar command to archive the object files
AR = ar

# include path
INCS  = -I$(MODULE_DIR)/Includes
INCS += -I$(PROJECT_DIR)/Includes
INCS += -I$(TOP)/Src/Main/Cores

# special includes
SPECIAL_INCS =

# library path
LIBS = -L$(TOP)/Lib

# special libs
SPECIAL_LIBS =

# system library
SYSTEM_LIBS = -lm

# inlcude MyConfig for more options and c source files
modules_path  =# to be determined in MyConfig
projects_path =# to be determined in MyConfig
include MyConfig

# all c file directories
c_dirs  = $(TOP)/Src/Main/Cores
c_dirs += $(modules_path)
c_dirs += $(projects_path)

# all c files patch
c_paths = $(foreach dir,$(c_dirs),$(wildcard $(dir)/*.c))

# all compiler flags
CFLAGS = $(DFLAGS) $(OFLAGS) $(WARN) $(INCS) $(SPECIAL_INCS)

# all linking flags:
LDFLAGS = $(LIBS) $(SPECIAL_LIBS) $(SYSTEM_LIBS)

# finding c files inter-dependencies using DEPFLAGS flag of the compiler
DEPFLAGS = -M
depend_paths = $(c_paths:.c=.depend)# this substitude .c with .depend
obj = $(c_paths:.c=.o)# this substitude .c with .depend



#.PHONY: compile
compile: find_interdependency $(obj)
	$(CC) $(CFLAGS) -o $(EXEC_DIR)/$(EXEC) $(obj) $(LDFLAGS)

.PHONY: find_interdependency
find_interdependency: $(depend_paths)
$(depend_paths):%.depend:%.c MyConfig# should I use myconfig here????
	@echo $@
	@set -e; rm -f $@;\
	$(CC) $(DEPFLAGS) $(CFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

.PHONY:print
print:find_interdependency
	@echo "starts here:\n"$(depend_paths) | tr " " "\n"

# now include all of the inter-dependency files
#include $(depend_paths)
include $(wildcard $(depend_paths))

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
	
# new line variable
define n


endef
