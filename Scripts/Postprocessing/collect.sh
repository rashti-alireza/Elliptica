# Alireza Rashti May 2020
# 
# This small script collects the given quantity from the given file
# and then put it into a table to be plotted. It also plots it if
# the tgraph.py app is available on the system.
# The format of the file should be like the followings:
# q0 = v1
# q1 = v2
# q0 = v3
# q2 = v4
# q0 = v5
# q0 = v6
# .
# .
# .
#
# thus it collects let's say all of q0's, and the put them into a file
# such that the first column is row number and the second is value 
# pertinent to that row.
#
# Note1: It is case insensitive.
# Note2: This script tends to mainly use *_properties.txt files.
# Note3: you should be in the directory that contains the file.
#
# ToDo: fix Note3.
# 
# usage:
# $collect.sh quantity_name file_name    # collect and plot
# $collect.sh -s quantity_name file_name # collect, scale and plot
#

#!/bin/bash

argc=$#
argv=("$@")
prefix=${argv[0]}
quant="\b${argv[0]}\b"
file=${argv[1]}
out=${prefix}_${file}
scale=0

# check if it needs help
if [[ $argc -le 1 || $1 =~ --h.? ]]; then
        printf \
"usage:\n"\
"$ collect.sh quantity_name file_name    # --> (collect and plot)\n"\
"$ collect.sh -s quantity_name file_name # --> (collect, scale and plot)\n"
	exit 1
# if asked scaling
elif [[ $1 =~ -s.? ]]; then
	scale=1
	file=${argv[2]}
	prefix=${argv[1]}
	quant="\b${argv[1]}\b"
	out=scaled_${prefix}_${file}
fi

# check if file exists
if [[ ! -f ${file} ]]; 
then
	printf "Could not open file %s\n" ${file}
	exit 2
fi

# now make the table:
printf "# iteration %s\n" ${prefix} > ${out}

# if scale required
if [[ ${scale} -eq 1 ]];then
	data=$(grep -i -E ${quant} ${file} | \
		sed 's/+//g' | \
		awk '{print $3}'| tr '\n' ',')
	IFS=',' read -r -a array <<< "${data}"
	# find the max to scale
	max=0
	for (( i=0; i < ${#array[@]}; ++i ));
	do
		r=${array[i]}
		if (( $(echo "$r < 0" | bc -l) ));
	        then
        	 	r=$(echo "$r * (-1)" | bc -l)
        	fi
        	
        	if (( $(echo "$r > $max" | bc -l) ));
        	then
        		max=$r
        	fi
	done
	# scale by dividing by max and print
	for (( i=0; i < ${#array[@]}; ++i ));
	do
		array[i]=$(echo "${array[i]}/${max}" | bc -l)
		echo $((i+1)) ${array[i]} >> ${out}
	done
# just print to file
else
	grep -i -E ${quant} ${file} | \
	awk 'BEGIN{v=0;}{ v = v+1; print v " " $3}' >> ${out} ;
fi

# print:
printf "\n%s  = %s\n\n" "==> output file" ${out}

# plot with tgraph:
tgraph.py -m ${out}


# End

