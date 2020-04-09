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
top=$(shell pwd)
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
dep_dir := $(top)/.dep
# dependency files
d_files := \
       $(foreach f,$(c_files),\
         $(join\
	   $(dep_dir)/,\
           $(notdir $(f:.c=.d))\
          )\
        )

#######################################################################
###########
## targets:
###########
# default target to make libaries out of object files
make_lib: compile_o
	@echo $(PR_F1)arching:
	$(AR) rcs $(LIB_DIR)/$(lib_name) $(o_files)
	
compile_o: $(o_files)
	@echo $(pr_f1) "object files made" $(pr_f2)

# using string % to make the object file accoding to its c file.
# then put the resultant into $(o_dir)/$*.o.
$(o_dir)/%.o: $(d_files) | $(dep_dir)
#	@echo Making object file for :
#	@echo $(pr_f1); echo $(top) | grep -o '/Src/*' $(pr_f2)
#	@echo
	$(CC) $(CFLAGS) -o $(o_dir)/$*.o -c $<

# figuring out the inter dependencies
#%.o : %.c
%.d:%.c
	set -e; rm -f $@;\
	$(CC) $(DEPFLAGS) $(CFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# making dep_dir if does not exist
$(dep_dir):
	@mkdir -p $@

include $(wildcard $(dep_dir)/*.d)

.PHONY: compile_o make_lib
########################################################################