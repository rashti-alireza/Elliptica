# Alireza Rashti July 2022
#
# find out the column of interest from the txt format 0d-outputs to 
# be plotted. the header (THE VERY FIRST LINE) of the output file should be like:
#
# # iteration field1|Linf field1|L1 field1|L2 field2|Linf ....
#
# NOTE: there should be a space between '#' and 'iteration' and 
# all other names.
# NOTE: it requires the bash file "postprocess_utils.sh" so add it into the PATH search
#
# usage:
# $ txt_0d_plot.sh --help 
# 
# NOTE: it assumes posix-extended for find regextype.
#

#!/usr/bin/env bash

## load utils
mypath=$(realpath $0)
mydir=$(dirname ${mypath})
source "${mydir}/plot_utils.sh"

## some defs
suffix0d="0d.txt"
outdir="Diagnostics_00"
col1_name="iteration"
res_regex="[[:digit:]]+x[[:digit:]]+x[[:digit:]]+"

argc=$#
argv=("$@")
## find column1 (col1) and coulmn2 (col2): (col2 vs col1)
col1="1" ## 1 is the iteration column
col2=""  ## we want to find this

# check if it needs help
if [[ $argc -le 1 || $1 =~ --h.? ]];
then
	pr_header "help"
	printf "Reading the plot files with '${suffix0d}' in the '${outdir}' directory.\n\n"
	
	pr_header "usage"
        printf "$ txt_0d_plot.sh <dir_output_name> <quantity> <region>\n\n"
        
	pr_header "examples"
	printf "## The following plots \"ham1|L2\" (L2 norm of Hamiltonian) vs \"iteration\".\n"
	printf "## The region is the \"left_NS_around_front.+\" files for all available resolutions.\n"
        printf "$ txt_0d_plot.sh bhns_00 \"ham1|L2\" \"left_NS_around_front.+\"\n\n"
        printf "## To only pick specific resolutions, e.g., 18 and 20, we have:\n"
        printf "$ txt_0d_plot.sh bhns_00 \"ham1|L2\" \"(18|20).+left_NS_around_front.+\"\n\n"
      
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

## collect all outdir inside the topdir
subdirs=$(find "${topdir}" -type d -name "${outdir}")

## collect all matched files in each outdir
colms=()
files=()
for subdir in ${subdirs[@]}
do
	matched_files=$(find "${subdir}" -type f -regextype posix-extended \
	                -regex ".+${argv[ $(($argc -1)) ]}${suffix0d}$" )
	if [[ ${#matched_files} -eq 0 ]];
	then
		printf "\n!!\nCould not find any match for \"${argv[ $(($argc -1)) ]}\" in:\n${subdir}\n"
		continue
	fi

	## find and save col2 in each file
	for f in ${matched_files[@]}
	do
		c=$(find_field_position_header ${field} $f)
		if [[ ${#c} -eq 0 ]];
		then
			printf "!!\ncound not find ${field} inside:\n$f\n"
			exit 2
		fi
		# Add new element at the end of the array
		colms+=($c)
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
	fname=$(echo ${fname}     | sed -E "s/\/${outdir}\/.+/_/g" )
	fname=$(echo ${fname}     | sed -E "s/_[[:digit:]]+_$//g"  )
	fname="${field}___${fname}__"
	tmp_file=$(mktemp ".${fname}XXXXX")
	## soft link
	ln -fs ${files[$i]} ${tmp_file}
	files_tmp+=("${tmp_file}")
done

echo "---"
## print pertinent info and save command for tgraph
graph_cmds=""
for ((i=0; i < ${#colms[@]}; i++))
do
	echo "${files[$i]}:"
	echo "plot: ${field} vs. ${col1_name} ==> column = ${colms[$i]} vs. column = ${col1}"
	echo "---"
	graph_cmds+=" -c ${col1}:${colms[$i]} ${files_tmp[$i]}"
	
done

## plot
tgraph.py -m ${graph_cmds}

## remove tmp files
rm_tmps

## end
