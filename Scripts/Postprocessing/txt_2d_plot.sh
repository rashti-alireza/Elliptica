# Alireza Rashti July 2022
#
# find out the columns of interest from the txt format 2d-outputs to 
# be plotted. the header (THE VERY FIRST LINE) of the output file should be like:
#
# plane_coordinate1 plane_coordinate2 x(X,Y,Z) y(X,Y,Z) z(X,Y,Z) field1 field2 ....
#
# NOTE: the should be a space between '#' and 'line_coordinate' and 
# all other names.
#
# usage:
# $ txt_2d_plot.sh --help
#
# NOTE: it assumes posix-extended for find regextype.
#

#!/usr/bin/env bash

## load utils
mypath=$(realpath $0)
mydir=$(dirname ${mypath})
source "${mydir}/plot_utils.sh"

## some defs:
suffix2d="2d.txt"
outdir="Diagnostics_00"
coord_default1="plane_coordinate1"
coord_default2="plane_coordinate2"
outdir="Diagnostics_00"
res_regex="[[:digit:]]+x[[:digit:]]+x[[:digit:]]+"

argc=$#
argv=("$@")

## find column1 (col1), column2 (col2) and column3 (col3)
coord1=""
coord2=""
field=""

# check if it needs help
if [[ $argc -le 1 || $1 =~ --hel.? ]];
then
	pr_header "help"
	printf "Reading the plot files with '${suffix2d}' in the '${outdir}' directory.\n\n"
	
	pr_header "usage"
	printf "$ txt_2d_plot.sh <dir_output_name> <coord1> <coord2> <quantity> <region>\n\n"
	
	pr_header "examples"
	
	printf "The following plots \"psi\" vs \"x(X,Y,Z)\" and \"y(X,Y,Z)\" over all surfaces.\n"
	printf "The region is all \"left_NS_(around_)?front.+\" files for all available resolutions.\n"
	printf "$ txt_2d_plot.sh bns_00 x y psi \"left_NS_(around_)?front.+\"\n\n"
	
	printf "Below, it plots only for 18 and 20 resolutions over the surface (X,Y,1):\n"
	printf "$ txt_2d_plot.sh bns_00 x y psi \"(18|20).+left_NS_(around_)?front.+_X_Y_1.+\"\n\n"
	
	printf "Plotting psi vs reference coords:\n"
	printf "$ txt_2d_plot.sh bns_00 psi \"(18|20).+left_NS_around_front.+_X_Y_1.+\" \n\n"

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

if [[ $argc -eq 5 ]];## ==> the coord is given
then
	coord1=${argv[1]}
	coord2=${argv[2]}
	
	## translate x,y,z
	case ${coord1} in
		"x") coord1="x(X,Y,Z)";;
		"y") coord1="y(X,Y,Z)";;
		"z") coord1="z(X,Y,Z)";;
		*) printf "Not expected ${coord1}!\n"; exit 2;;
		esac
	## translate x,y,z
	case ${coord2} in
		"x") coord2="x(X,Y,Z)";;
		"y") coord2="y(X,Y,Z)";;
		"z") coord2="z(X,Y,Z)";;
		*) printf "Not expected ${coord2}!\n"; exit 2;;
		esac
	field=${argv[3]}
elif [[ $argc -eq 3 ]];
then
	coord1="${coord_default1}"
	coord2="${coord_default2}"
	field=${argv[1]}
else
	printf "Too few arguments!\n"
        exit 2
fi

## collect all outdir inside the topdir
subdirs=$(find "${topdir}" -type d -name "${outdir}")

## collect all matched files in each outdir
field_colms=()
coord1_colms=()
coord2_colms=()
files=()
for subdir in ${subdirs[@]}
do
	matched_files=$(find "${subdir}" -type f -regextype posix-extended \
	               -regex ".+${argv[ $(($argc -1)) ]}${suffix2d}$" )
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
		
		c=$(find_field_position_header ${coord1} $f)
		if [[ ${#c} -eq 0 ]];
		then
			printf "\n!!\ncound not find \"${coord1}\" inside:\n$f\n"
			exit 2
		fi
		# Add a new element at the end of the array
		coord1_colms+=($c)
		
		c=$(find_field_position_header ${coord2} $f)
		if [[ ${#c} -eq 0 ]];
		then
			printf "\n!!\ncound not find \"${coord2}\" inside:\n$f\n"
			exit 2
		fi
		# Add a new element at the end of the array
		coord2_colms+=($c)
		
		# now add the file
		files+=($f)
	done
	
done

## remove tmp files
function rm_tmps(){
	for ((i=0; i < ${#files[@]}; i++))
	do
	       rm -rf ${files_tmp[$i]}
	done
}

## clean up after yourself
trap rm_tmps SIGHUP SIGINT SIGQUIT SIGABRT

## for a quick legend
files_tmp=()
for ((i=0; i < ${#files[@]}; i++))
do
	fname=$(echo ${files[$i]} | grep -E -o "${res_regex}.+"    )
	fname=$(echo ${fname}     | sed -E "s/_[[:digit:]]+\/${outdir}\//_/g" )
	fname="${field}___${fname}__"
	tmp_file=$(mktemp ".${fname}XXXXX")
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
	echo "plot: ${field} vs. (${coord1}, ${coord2}) ==> column = ${field_colms[$i]} vs. columns = (${coord1_colms[$i]}, ${coord2_colms[$i]})"
	echo "---"
	graph_cmds+=" -c ${coord1_colms[$i]}:${coord2_colms[$i]}:${field_colms[$i]} ${files_tmp[$i]}"
	
done

## plot
tgraph.py -m ${graph_cmds}

## remove tmp files
rm_tmps

## end
