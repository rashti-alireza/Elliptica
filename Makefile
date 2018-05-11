#Makefile
#Alireza Rashti, May 2018

#Makefile Directory
TOP := $(shell pwd)

#Compiler
CC = gcc

#gcc flags
GCCFLAGS = -Wall -g

#Include header files from different libraries
INCLUDE = -I $(TOP)/Src/Main/Modules/Libraries

#Compiler flags
CFLAGS  = $(GCCFLAGS)
CFLAGS += $(INLCUDE)

#Linking Libraries
LDFLAGS  = -lumfpack -lblas -lgfortran -llapack
LDFLAGS += -lfftw3
LDFLAGS += -lm

#Sources:
src = $(TOP)

#Compiling abc - default target
.PHONY: abc
abc: $(obj)
	$(CC) $(CFLAGS) -o $@ $? $(LDFLAGS)

#Cleaning the whole object and binary files 
.PHONY: clean
clean:
	-rm -rf $(TOP)/Bin/* $(TOP)/Lib/*