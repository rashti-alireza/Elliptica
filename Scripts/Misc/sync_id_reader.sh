#!/usr/bin/env bash

###########
## usage ##
###########
## $ ./me /path/to/Elliptica /path/to/Elliptica_ID_Reader/dir
##
## a quick ad-hoc method to collect a subset of source files needed for ID reader
##
## NOTE: check that all projects be in the working git branch (not testing ones)
##

E_TOP="$1"
IDR_TOP="$2"

##### TEMP
rm -rf ${IDR_TOP}/src/*

SRC=()

## Core
CORE_DIR="Src/Main/Core"
SRC+=("${CORE_DIR}")

## modules
MODULE_DIR="Src/Main/Modules"
SRC+=("${MODULE_DIR}/Manifold")
SRC+=("${MODULE_DIR}/Fields")
SRC+=("${MODULE_DIR}/Error_Handlings")
SRC+=("${MODULE_DIR}/Text_and_File_Tools")
SRC+=("${MODULE_DIR}/Prints")
SRC+=("${MODULE_DIR}/Utilities")
SRC+=("${MODULE_DIR}/Checkpoint")
#SRC+=("${MODULE_DIR}/Maths/Linear_Algebra")
SRC+=("${MODULE_DIR}/Maths/Field_Analysis")
SRC+=("${MODULE_DIR}/Maths/Equation_Solvings")
#SRC+=("${MODULE_DIR}/Maths/Matrix_Solvers")
SRC+=("${MODULE_DIR}/Maths/General")
SRC+=("${MODULE_DIR}/Maths/Spectral_Methods")
SRC+=("${MODULE_DIR}/Maths/Special_Functions")
SRC+=("${MODULE_DIR}/Maths/Calculus")
SRC+=("${MODULE_DIR}/Maths/Complex")
SRC+=("${MODULE_DIR}/Maths/Diff_Geom")
SRC+=("${MODULE_DIR}/Physics")
SRC+=("${MODULE_DIR}/Physics/EoS")
SRC+=("${MODULE_DIR}/Physics/Observe")
SRC+=("${MODULE_DIR}/Physics/Transformation")
SRC+=("${MODULE_DIR}/Physics/StressEnergy_Tensor")
SRC+=("${MODULE_DIR}/Physics/Star")
SRC+=("${MODULE_DIR}/Physics/BlackHole")
#SRC+=("${MODULE_DIR}/Physics/System")
SRC+=("${MODULE_DIR}/Physics/Free_Data")
SRC+=("${MODULE_DIR}/Physics/ADM")
#SRC+=("${MODULE_DIR}/Physics/Equation")
SRC+=("${MODULE_DIR}/Includes")

## projects
PROJ_DIR="Src/Projects"
SRC+=("${PROJ_DIR}/Includes")
SRC+=("${PROJ_DIR}/Initial_Data_Reader")
SRC+=("${PROJ_DIR}/BH_NS_Binary_Initial_Data")
#SRC+=("${PROJ_DIR}/NS_NS_Binary_Initial_Data")
SRC+=("${PROJ_DIR}/TOV_star")



## src dir
mkdir -vp ${IDR_TOP}/src

## find all *.c and *.h files
fs=()
fh=()
for d in ${SRC[@]}
do
	fs+=($(find "${E_TOP}/$d" -maxdepth 1 -type f -name "*.c" ))
	fh+=($(find "${E_TOP}/$d" -maxdepth 1 -type f -name "*.h" ))
done

for f in ${fs[@]}
do
	cp -uv $f ${IDR_TOP}/src/
done

for f in ${fh[@]}
do
	cp -uv $f ${IDR_TOP}/src/
done

## remove redundant files
rm -v ${IDR_TOP}/src/main.?
rm -v ${IDR_TOP}/src/bhns_analyze.?
rm -v ${IDR_TOP}/src/bhns_main.?
rm -v ${IDR_TOP}/src/bhns_solve_eqs.?
rm -v ${IDR_TOP}/src/projects_data_base_MADE_BY_MAKE.?
rm -v ${IDR_TOP}/src/pr_for_fields.?
rm -v ${IDR_TOP}/src/solve_eqs_ddm_schur_complement.?
rm -v ${IDR_TOP}/src/solve_eqs.?
rm -v ${IDR_TOP}/src/solvings_tests.?
rm -v ${IDR_TOP}/src/solve_manager.?
rm -v ${IDR_TOP}/src/dfs_df.?

## remove uncalling func
cd ${IDR_TOP}/src

sed -i 's/test_print(PRINT_COORDS)/0/g' *.c
sed -i '/move_dfdu_jacobian_patch/d' *.c
sed -i '/free_patch_SolMan_jacobian/d' *.c
sed -i '/free_patch_SolMan_method_Schur/d' *.c
sed -i '/alloc_matrix/d' *.c
sed -i '/free_matrix/d' *.c
sed -i '/eq_main/d' *.c
sed -i '/sys_main/d' *.c


