# My Configurations:

# modules to be compiled:
modules_path += $(TOP)/Src/Main/Modules/Manifold
modules_path += $(TOP)/Src/Main/Modules/Fields
modules_path += $(TOP)/Src/Main/Modules/Error_Handlings
modules_path += $(TOP)/Src/Main/Modules/Text_and_File_Tools
#modules_path += $(TOP)/Src/Main/Modules/Includes
modules_path += $(TOP)/Src/Main/Modules/Memory_Managers
modules_path += $(TOP)/Src/Main/Modules/Prints
modules_path += $(TOP)/Src/Main/Modules/Utilities
modules_path += $(TOP)/Src/Main/Modules/Maths/Linear_Algebra
modules_path += $(TOP)/Src/Main/Modules/Maths/Field_Analysis
modules_path += $(TOP)/Src/Main/Modules/Maths/Equation_Solvings
modules_path += $(TOP)/Src/Main/Modules/Maths/Matrix_Solvers
modules_path += $(TOP)/Src/Main/Modules/Maths/General
modules_path += $(TOP)/Src/Main/Modules/Maths/Approximation
modules_path += $(TOP)/Src/Main/Modules/Maths/Analytic
modules_path += $(TOP)/Src/Main/Modules/Maths/Calculus
modules_path += $(TOP)/Src/Main/Modules/Maths/Complex_Numbers
modules_path += $(TOP)/Src/Main/Modules/Physics/StressEnergy_Tensor
modules_path += $(TOP)/Src/Main/Modules/Physics/EoS
modules_path += $(TOP)/Src/Main/Modules/Physics/Observables
modules_path += $(TOP)/Src/Main/Modules/Physics/Transformation

# Projects to be compiled:
projects_path += $(TOP)/Src/Projects/Projects_Setup
projects_path += $(TOP)/Src/Projects/Poisson/Laplace_Inhom
projects_path += $(TOP)/Src/Projects/Tests/Modules_Test
projects_path += $(TOP)/Src/Projects/TOV_Stars
projects_path += $(TOP)/Src/Projects/Binary_BH_NS_Initial_Data



# special includes
#SPECIAL_INCS +=

# special libraries
#SPECIAL_LIBS +=

# compiler
#CC = clang
CC = gcc

# compiler flags
OFLAGS = -fopenmp -g3
#OFLAGS = -fopenmp -O3
#OFLAGS = -fopenmp -g3 #-g3 #-O3 #-Wall -W # mild

# warninng flags for the compiler
WARN =# no warnings
#WARN = -std=c11 -pedantic -Wextra -Werror -Wall -W \
  -Wmissing-prototypes -Wstrict-prototypes\
  -Wconversion -Wshadow -Wpointer-arith\
  -Wcast-qual -Wcast-align\
  -Wwrite-strings -Wnested-externs\
  -fshort-enums -fno-common -Dinline= # -O3 #-fmax-errors=1 #-O2 # aggressive reliable C 11

# define flags for compiler
DEFFLAGS = -DPragma_OpenMP_2d \
	   -DPragma_OpenMP_1d -DPragma_OpenMP_Patch


