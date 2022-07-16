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
# $ txt_1d_plot.sh "x(X,Y,Z)" "psi" <file_name> ## plot "psi" vs "x(X,Y,Z)"
# $ txt_1d_plot.sh "psi" <file_name>            ## plot "psi" vs line_coordinate
#

#!/bin/bash

argc=$#
argv=("$@")
## find column1 (col1) and coulmn2 (col2):
col1=""
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
"$ txt_1d_plot.sh \"x(X,Y,Z)\" \"psi <file_name>\" #=> plot psi vs x(X,Y,Z)\n"\
"$ txt_1d_plot.sh \"psi\" <file_name> #=> plot psi vs line_coordinate\n"
        exit 1
fi

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




