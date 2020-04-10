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
# library name for this directory this is a standard name
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

## default target to make libaries out of object files
make_lib: $(o_files) $(d_files)
	@$(call build_and_print_func, $(AR) $(ARFLAGS) $(LIB_DIR)/$(lib_name) $(o_files), $(lib_name))

## using string % to make the object file according to its c file.
# then put the resultant into $(o_dir)/$*.o.
$(o_dir)/%.o :$(c_dir)/%.c  $(d_files)
	@$(call build_and_print_func, $(CC) $(CFLAGS) -o $(o_dir)/$*.o -c $<, $*.o)
	
## figuring out the inter dependencies
%.o : %.c
$(d_dir)/%.d : $(c_dir)/%.c | $(d_dir)
	@set -e; rm -f $@;\
	$(CC) $(DEPFLAGS) $(CFLAGS) $< > $@.$$$$; \
	$(SED) 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# making dep_dir if does not exist
$(d_dir):
	@mkdir -p $@

include $(wildcard $(d_dir)/*.d)

.PHONY: make_lib
########################################################################
#############
## functions:
#############

# build and print objects, if error happens it prints full info,
# otherwise it prints succinctly.
# useage:
#
# $(call build_and_print_func, command, print_this)
# $(1) refers to command and $(2) refers to print_this.
define build_and_print_func
$(strip $(1)) 2> $@.COMPILE_ERROR ; \
ret=$$? ; \
if [ $$ret -eq 0 ]; \
then \
printf "%s %s %-39s %s\n" $(PR_F0) "building" $(strip $(2)) "=> [SUCCESSFUL] <="; \
else \
echo $(PR_DIVIDER) ; \
printf "%s %s %-39s %s\n" $(PR_F0) "building" $(strip $(2)) "=> [FAILED] <="; \
cat $@.COMPILE_ERROR ; \
fi ; \
rm -rf $@.COMPILE_ERROR; \
exit $$ret ;
endef

# no export 
unexport
