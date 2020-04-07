##################
# Alireza Rashti #
#                #
# April 2020     #
##################

# Using GNU make to compile a software written in C language.
# It uses MyConfig file to get the source files and then 
# detects the dependencies automatically and finally after constructing
# the libraries (shared or static) it makes the exeutable output.
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
ifeq ($(LIB),)
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



all:
	@echo "Lib=" $(Lib)
	
# new line variable
define n


endef
