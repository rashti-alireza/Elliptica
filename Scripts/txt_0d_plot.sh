# Alireza Rashti July 2022
#
# find out the columns of interest from the txt format 0d-outputs to 
# be plotted. the header (THE VERY FIRST LINE) of the output file should be like:
#
# # iteration field1|Linf field1|L1 field1|L2 field2|Linf ....
#
# NOTE: there should be a space between '#' and 'iteration' and 
# all other names.
#
# usage:
## to plot "ham1|L2" vs "iteration" for all resolutions at all "left_NS_around_front*" files:
# $ txt_0d_plot.sh  <dir_output_name> "ham1|L2" "left_NS_around_front.+"
#

#!/bin/bash

## load utils
source "./postprocess_utils.sh"

suffix0d="0d.txt"
outdir="Diagnostics_00"

argc=$#
argv=("$@")
## find column1 (col1) and coulmn2 (col2): (col2 vs col1)
col1="1"
col2="" 
## v vs x
x=""
v=""

# check if it needs help
if [[ $argc -le 1 || $1 =~ --h.? ]];
then
        printf \
"usage:\n"\
"------\n"\
"## to plot \"ham1|L2\" vs \"iteration\" for all resolutions at all \"left_NS_around_front*\" files \n"\
"$ txt_0d_plot.sh  <dir_output_name> \"ham1|L2\" \"left_NS_around_front.+\" \n\n"
        exit 1
fi

topdir=${argv[0]}
field=${argv[1]}

# check if topdir exists
if [[ ! -d ${topdir} ]]; 
then
        printf "Could not open topdir \"%s\"\n" ${topdir}
        exit 2
fi

## collect all subdirs inside the topdir inside outdir
subdirs=$(find "${topdir}" -type d -name "${outdir}")

## collect all matched files in each subdir

for subdir in ${subdirs[@]}
do
	matched_files=$(find "${subdir}" -type f -regex ".+${argv[ $(($argc -1)) ]}${suffix0d}$" )
	if [[ ${#matched_files} -eq 0 ]];
	then
		printf "Could not find any match for \"${argv[ $(($argc -1)) ]}\" file!\n"
		exit 2
	fi

	## find col2 in each file and save it
	colms=()
	for f in ${matched_files[@]}
	do
		echo $f
		c=$(find_field_position_header ${field} $f)
		if [[ ${#c} -eq 0 ]];
		then
			printf "cound not find ${field} inside $f\n"
			exit 2
		fi
		# Add new element at the end of the array
		echo $c
		#colms+=$c
	done
	
	## print pertinent found info
	for c in ${colms[@]}
	do
		echo $c
	done
	
	echo "----"
	
done


exit 0


## if 3 args given
counter=0
if [[ $argc -eq 3 ]];
then
	header=$(head -n1 ${file})
	coord=${argv[0]}
	field=${argv[1]}
	x=$coord
	v=$field
	
	echo $header
	for i in $header
	do
		if [[ "$coord" == "$i" ]];
		then
			col1=$counter
		elif [[ "$field" == "$i" ]]
		then
			col2=$counter
		fi
		## test
		echo $i $counter
		
		counter=$((counter+1))
	done
        if [[ "$col2" == "" ]];
        then
                echo "could not find the field name \"${field}\""
                exit 2
        fi

fi

counter=0
if [[ $argc -eq 2 ]];
then
	col1=1
	header=$(head -n1 ${file})
	field=${argv[0]}
	x="line_coordinate"
	v=$field
	
	echo $header
	for i in $header
	do
		if [[ "$field" == "$i" ]]
		then
			col2=$counter
		fi
		## test
		echo $i $counter
		
		counter=$((counter+1))
	done
	if [[ "$col2" == "" ]];
        then
                echo "could not find the field name \"${field}\""
                exit 2
        fi
fi

##
echo "---"
echo "plot: $v  vs. $x ==> column = $col2 vs. column = $col1"
echo "---"

tgraph.py -m -c ${col1}:${col2} ${file}




