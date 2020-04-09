##################
# Alireza Rashti #
#                #
# April 2020     #
##################

# a general make file for each module or project.
# this make file gathers all of c files and make object files
# from them and put them in the specified directory. finally,
# by using (AR) command, it archive all of the object files to later
# be linked.
# NOTE: local variables are in small words and exported variables
# by the main Makefile are with capital words.
#######################################################################
#############
## variables:
#############
# where I am:
top   := $(shell pwd)
# c directory
c_dir := $(top)
# find all of c files to be compiled in this directory
c_files := $(wildcard $(top)/*.c)
# make directory name correspond to this directory
o_dir := $(notdir $(top))
o_dir := $(O_TOP)/$(o_dir)
# library name for this directory
lib_name := lib$(notdir $(top)).a
# collect all o files related to all of the c files in this directory
o_files:= \
       $(foreach f,$(c_files),\
         $(join\
           $(addprefix \
             $(O_TOP)/, $(notdir $(top))\
            )/,\
           $(notdir $(f:.c=.o))\
          )\
        )
# dependency directory
d_dir := $(top)/.dep
# dependency files
d_files := \
       $(foreach f,$(c_files),\
         $(join\
	   $(d_dir)/,\
           $(notdir $(f:.c=.d))\
          )\
        )

#######################################################################
###########
## targets:
###########
# default target to make libaries out of object files
make_lib: $(o_files) $(d_files)
	$(AR) rcs $(LIB_DIR)/$(lib_name) $(o_files)
	
compile_o: $(o_files)
#	@echo $(PR_F1) $(c_dir) $(PR_F2)
#	@echo $(PR_F1) $(c_files) $(PR_F2)
	@echo $(PR_F1) $(o_files) $(PR_F2)
#	@echo $(PR_F1) $(o_dir) $(PR_F2)

# using string % to make the object file according to its c file.
# then put the resultant into $(o_dir)/$*.o.
$(o_dir)/%.o :$(c_dir)/%.c  $(d_files) 
	$(CC) $(CFLAGS) -o $(o_dir)/$*.o -c $<
#
#
#	@echo Making object file for :
#	@echo $(PR_F1); echo $(top) | grep -o '/Src/*' $(PR_F2)
#	@echo
#	
# figuring out the inter dependencies
#%.o : %.c
$(d_dir)/%.d : $(c_dir)/%.c | $(d_dir)
	set -e; rm -f $@;\
	$(CC) $(DEPFLAGS) $(CFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# making dep_dir if does not exist
$(d_dir):
	@mkdir -p $@

include $(wildcard $(d_dir)/*.d)

#.PHONY: compile_o make_lib
########################################################################