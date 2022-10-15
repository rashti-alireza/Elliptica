# Alireza Rashti July 2022
#
# variety of utilities for postprocessing to be used by other scripts
#


## find the field position $1 on the header of given file $2
## using echo, it returns the column number starting from 0
function find_field_position_header {
	local counter=0
	local header=$(head -n1 $2)
	local ret=""
	
	for i in $header
	do
		if [[ "$1" == "$i" ]]
		then
			ret=$counter
		fi
		counter=$((counter+1))
	done
	echo -n "${ret}"
	
	return 0
}

#export -f find_field_position_header

## printing fancy headers!
function pr_header(){
        local header_line="#########################################################"
        local len
        (( len=${#1} + 6))
        echo "${header_line:0:$len}"
        ## to capitalize use ^^
        echo "${header_line:0:2} ${1^^} ${header_line:0:2}"
        echo "${header_line:0:$len}"
}

