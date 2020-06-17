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
# Note: It is case insensitive.
# Note: This script tends to mainly use *_geometry_and_physics.txt files.
#
# usage:
# $gather_quantity.sh quantity_name file_name
#

#!/bin/bash

argc=$#
argv=("$@")
prefix=${argv[0]}
quant="\b${argv[0]}\b"
file=${argv[1]}
out=${prefix}_${file}

# check if it needs help
if [[ $argc -le 1 || $1 =~ --h.? ]]; then
        printf \
"usage:\n"\
"$ gather_quantity.sh quantity_name file_name\n"
	exit 1
fi

# check if file exists
if [[ ! -f ${file} ]]; 
then
	printf "Could not open file %s\n" ${file}
	exit 2
fi

# now make the table:
printf "# iteration %s\n" ${quant} > ${out}
grep -i -E ${quant} ${file} | \
awk 'BEGIN{v=0;}{ v = v+1; print v " " $3}' >> ${out} ;

# print:
printf "\n%s  = %s\n\n" "==> output file" ${out}

# plot with tgraph:
tgraph.py -m ${out}




# End

