#Makefile
#Alireza Rashti, May 2018

#Makefile Directory
TOP := $(shell pwd)

#Lib directory
Lib := $(TOP)/Lib
#Including MyConfig which has projects and gcc flags among others
include MyConfig

#Including Modules:
include MyModules

#Searching path for module libraries
INCLUDE  = -I$(TOP)/Src/Main/Modules/Libraries
INCLUDE2 = -I$(TOP)/Src/Main/Cores

#Compiler flags
CFLAGS = $(GCCFLAGS) #in MyConfig

#Definition flags
DFLAGS = $(DEFFLAGS) #in MyConfig

#Linking Libraries
LDFLAGS  = -lumfpack -lblas -lgfortran -llapack
LDFLAGS += -lfftw3_omp -lfftw3
LDFLAGS += -lm
LDFLAGS += -lsiloh5

##Note: the following making object is very rudimentary
## I need to work on it more later
#Making Object files
c_src += $(foreach dir,$(modules_path),$(wildcard $(dir)/*.c))
c_src += $(foreach dir,$(projects_path),$(wildcard $(dir)/*.c))
obj = $(c_src:.c=.o)
	
###Targets###

#Compiling abc - default target
.PHONY: abc
$(EXE): $(c_src)
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDE) $(INCLUDE2) -o $(EXEDIR)/$(EXE) $? $(LDFLAGS)

#Cleaning the whole object and binary files 
.PHONY: clean
clean:
	-rm -rf $(TOP)/Bin/$(EXE) $(TOP)/$(Lib)/*
	$(foreach dir, $(modules_path), rm -rf $(dir)/*.o)
	$(foreach dir, $(projects_path), rm -rf $(dir)/*.o)
	$(foreach dir, $(modules_path), rm -rf $(dir)/*.*~)
	$(foreach dir, $(projects_path), rm -rf $(dir)/*.*~)

