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
# NOTE: it assumes posix-extended for find regextype.
#

#!/bin/bash

## load utils
source "plot_utils.sh"

## some defs:
suffix1d="1d.txt"
outdir="Diagnostics_00"
coord_default="line_coordinate"
outdir="Diagnostics_00"
res_regex="[[:digit:]]+x[[:digit:]]+x[[:digit:]]+"

argc=$#
argv=("$@")

## find column1 (col1) and column2 (col2):
coord=""
field=""

# check if it needs help
if [[ $argc -le 1 || $1 =~ --hel.? ]];
then
	pr_header "help"
	printf "Reading the plot files with '${suffix1d}' in the '${outdir}' directory.\n\n"
	
	pr_header "usage"
	printf "$ txt_1d_plot.sh <dir_output_name> <coord> <quantity> <region>\n\n"
	
	pr_header "examples"
	
	printf "## The following plots \"psi\" vs \"x(X,Y,Z)\" over all lines.\n"
	printf "## The region is all \"left_NS_(around_)?front.+\" files for all available resolutions.\n"
	printf "$ txt_1d_plot.sh bns_00 x psi \"left_NS_(around_)?front.+\"\n\n"
	
	printf "Below plots only for 18 and 20 resolutions on the line (0.5,0.5,Z):\n"
	printf "$ txt_1d_plot.sh bns_00 x psi \"(18|20).+left_NS_front.+_0.5_0.5_Z.+\"\n\n"
	
	printf "## Plotting psi vs referee coords:\n"
	printf "$ txt_1d_plot.sh bns_00 psi \"left_NS_front.+\" \n\n"

	pr_header "extra"
	
	printf "A rough translation of the reference coordinate (X,Y,Z) used in each\n"\
"cubed spherical patch to the Cartesian coordinates listed below.\n"\
"Note: Z always increases in the radial direction w.r.t the slice.\n"\
"\n"\
"up    : X = x, Y = y, Z = z\n"\
"down  : X = y, Y = x, Z = z\n"\
"\n"\
"left  : X = x, Y = z, Z = y\n"\
"right : X = z, Y = x, Z = y\n"\
"\n"\
"back  : X = z, Y = y, Z = x\n"\
"front : X = y, Y = z, Z = x\n\n"\

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

if [[ $argc -eq 4 ]];## ==> the coord is given
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
	matched_files=$(find "${subdir}" -type f -regextype posix-extended \
	               -regex ".+${argv[ $(($argc -1)) ]}${suffix1d}$" )
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

## for a quick legend
files_tmp=()
for ((i=0; i < ${#files[@]}; i++))
do
	fname=$(echo ${files[$i]} | grep -E -o "${res_regex}.+"    )
	fname=$(echo ${fname}     | sed -E "s/_[[:digit:]]+\/${outdir}\//_/g" )
	fname="${field}___${fname}__"
	tmp_file=$(mktemp ".${fname}XXXXX.txt")
	## soft link
	ln -fs ${files[$i]} ${tmp_file}
	files_tmp+=("${tmp_file}")
done

echo "---"
## print pertinent info and save command for tgraph
graph_cmds=""
for ((i=0; i < ${#files[@]}; i++))
do
	echo "${files[$i]}:"
	echo "plot: ${field} vs. ${coord} ==> column = ${field_colms[$i]} vs. column = ${coord_colms[$i]}"
	echo "---"
	graph_cmds+=" -c ${coord_colms[$i]}:${field_colms[$i]} ${files_tmp[$i]}"
	
done

## plot
tgraph.py -m ${graph_cmds}

## remove tmp files
for ((i=0; i < ${#files[@]}; i++))
do
	rm -rf ${files_tmp[$i]}
done

## end
