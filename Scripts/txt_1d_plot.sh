# Alireza Rashti July 2022
#
# find out the columns of interest from the txt format 1d-outputs to 
# be plotted. the header (THE VERY FIRST LINE) of the output file should be like:
#
# line_coordinate x(X,Y,Z) y(X,Y,Z) z(X,Y,Z) field1 field2 ....
#
# NOTE: the should be a space between '#' and 'line_coordinate' and 
# all other names.
#
# usage:
# $ txt_1d_plot.sh --help
#

#!/bin/bash

## load utils
source "postprocess_utils.sh"

## some defs:
suffix0d="1d.txt"
outdir="Diagnostics_00"
coord_ref_file="reference_coords_help.txt"
coord_default="line_coordinate"
outdir="Diagnostics_00"

argc=$#
argv=("$@")

## find column1 (col1) and coulmn2 (col2):
col1=""
col2="" 


# check if it needs help
if [[ $argc -le 1 || $1 =~ --h.? ]];
then
        printf \
"usage:\n"\
"------\n"\
"## to plot psi vs x(X,Y,Z) for all resolutions at all \"left_NS_(around_)?front.+\" files:\n"\
"$ txt_1d_plot.sh <dir_output_name> x psi \"left_NS_(around_)?front.+\"\n\n"\
"## to plot psi vs referece coords. for all resolutions at all \"left_NS_front.+\" files:\n"\
"$ txt_1d_plot.sh <dir_output_name> psi \"left_NS_front.+\" \n\n"

	cat "$coord_ref_file"
        exit 1
fi

## parse and check inputs:
topdir=${argv[0]}
# check if topdir exists
if [[ ! -d ${topdir} ]]; 
then
        printf "Could not open topdir \"%s\"\n" ${topdir}
        exit 2
fi

if [[ $argc -eq 4 ]];
then
	coord=${argv[1]}
	## translate x,y,z
	case ${coord} in
		"x") coord="x(X,Y,Z)";;
		"y") coord="y(X,Y,Z)";;
		"z") coord="z(X,Y,Z)";;
		*) printf "Not expected ${coord}!\n"; exit 2;;
		esac
	field=${argv[2]}
elif [[ $argc -eq 3 ]];
then
	coord="${coord_default}"
	field=${argv[1]}
else
	printf "Too few arguments!\n"
        exit 2
fi

## collect all outdir inside the topdir
subdirs=$(find "${topdir}" -type d -name "${outdir}")

## collect all matched files in each outdir
field_colms=()
coord_colms=()
files=()
for subdir in ${subdirs[@]}
do
	matched_files=$(find "${subdir}" -type f -regex ".+${argv[ $(($argc -1)) ]}${suffix0d}$" )
	if [[ ${#matched_files} -eq 0 ]];
	then
		printf "!!\nCould not find any match for \"${argv[ $(($argc -1)) ]}\" in\n${subdir}\n"
		continue
	fi

	## find and save col2 in each file
	for f in ${matched_files[@]}
	do
		c=$(find_field_position_header ${field} $f)
		if [[ ${#c} -eq 0 ]];
		then
			printf "\n!!\ncound not find \"${field}\" inside:\n$f\n"
			exit 2
		fi
		# Add a new element at the end of the array
		field_colms+=($c)
		
		c=$(find_field_position_header ${coord} $f)
		if [[ ${#c} -eq 0 ]];
		then
			printf "\n!!\ncound not find \"${coord}\" inside:\n$f\n"
			exit 2
		fi
		# Add a new element at the end of the array
		coord_colms+=($c)
		
		# now add the file
		files+=($f)
	done
	
done
exit 0




file=${argv[ $(($argc-1)) ]}
# check if file exists
if [[ ! -f ${file} ]]; 
then
        printf "Could not open file %s\n" ${file}
        exit 2
fi

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




