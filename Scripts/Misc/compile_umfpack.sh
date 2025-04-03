#!/usr/bin/env bash
## Alireza Rashti - Apr 2025 (C)
set -ue
#set -x

: <<'END_COMMENT'
usage:
# gnu compiler
./me -c gcc
# intel compilers, old versions
./me -c icc
# intel compilers, new versions
./me -c icx
END_COMMENT

function usage {
	echo "usage:"
	echo "./me --help"
	echo './me -c gcc'
	exit 1
}

function argparser {
  local OPTS

  OPTS="$(getopt -o c: -l help -- "$@")"
  if [[ $? != 0 ]]
  then
    echo "ERROR parsing arguments!" >&2
    exit 1
  fi
  
  if [[ "$#" -eq 0 ]]
  then
  	usage
  fi

  eval set -- "$OPTS" # set $1, $2, etc
  while true ; do
    case "$1" in
      "-h" | "--help")
        usage;;

      "-c")
        cc="$2";
        shift 2;
        break;;
      
      *)
      	usage;;
    esac
  done
}

argparser "$@"

## parallel job
pj=${pj:=4}

top="$(pwd)"
ff=""

# set fortran compiler
if [[ $cc = "gcc" ]]
then
	ff='gfortran'
	cpp='g++'
elif [[ $cc = "icc" ]]
then
	ff='ifort'
	cpp='icpc'
elif [[ $cc = "icx" ]]
then
	ff='ifx'
	cpp='icx'
else
	echo "no such option"
	exit 1
fi

############
## LAPACK ##
############
## ---------------------------------------------------------------------- ##
tar -zxvf v3.9.0.tar.gz || exit 1
mv lapack-3.9.0 "lapack-3.9.0_$cc"
cd "lapack-3.9.0_$cc" || exit 1
lapack_top="$(realpath .)"
cp make.inc.example make.inc || exit 1

echo "compiling LAPACK and BLAS using $cc:"
# compiler
sed -E -i "s/^CC = gcc/CC = $cc/" make.inc
sed -E -i "s/^FC = gfortran/FC = $ff/" make.inc

# optimization
sed -E -i "s/^FFLAGS = -O2/FFLAGS = -O3/" make.inc

# timer
sed -E -i "s/^TIMER = INT_ETIME/TIMER = NONE/" make.inc

# turn off testings
sed -i "s/all: lapack_install lib blas_testing lapack_testing/all: lapack_install lib blaslib cblaslib/" Makefile
sed -i "s/lib: lapacklib tmglib/lib: lapacklib/" Makefile

# intel specifics
if [[ $cc = "icc" || $cc = "icx" ]]
then
	sed -E -i "s/-frecursive//" make.inc
fi

# libs
sed -i 's:^BLASLIB      = $(TOPSRCDIR)/librefblas.a:BLASLIB      = $(TOPSRCDIR)/libblas.a:' make.inc
echo "CFLAGS += -fPIC" >> make.inc
echo "FFLAGS += -fPIC" >> make.inc

make -j"${pj}"

## ---------------------------------------------------------------------- ##

#############
## UMFPACK ##
#############
cd "$top"

tar -zxvf SuiteSparse-5.7.2.tar.gz
mv SuiteSparse-5.7.2 "SuiteSparse-5.7.2_$cc"
cd "SuiteSparse-5.7.2_$cc"
umfpack_top="$(realpath .)"

cd SuiteSparse_config

mv SuiteSparse_config.mk temp.mk
# set compilers
echo "CC = $cc" >> SuiteSparse_config.mk
echo "CXX = $cpp" >> SuiteSparse_config.mk
echo "F77 = $ff"  >> SuiteSparse_config.mk

if [[ $cc = "gcc" ]]
then
	echo "CFOPENMP = -fopenmp" >> SuiteSparse_config.mk
	echo "LDFLAGS += -fopenmp" >> SuiteSparse_config.mk
	echo "LDLIBS += -lgfortran -lm "  >> SuiteSparse_config.mk
elif [[ $cc = "icc" || $cc = "icx" ]]
then
	echo "CFLAGS += -D_GNU_SOURCE" >> SuiteSparse_config.mk
	echo "CFOPENMP = -qopenmp" >> SuiteSparse_config.mk
	echo "LDFLAGS += -qopenmp" >> SuiteSparse_config.mk
	echo "LDLIBS += -lm -lirc -lifcore -ldl"  >> SuiteSparse_config.mk
else
	echo "no option for compiler: $cc"
	exit 1
fi

# libs
echo "BLAS = -L${lapack_top} -lblas" >> SuiteSparse_config.mk
echo "LAPACK=-L${lapack_top} -llapack" >> SuiteSparse_config.mk
cat temp.mk >> SuiteSparse_config.mk
rm temp.mk

# turn off cuda
sed -i 's/CUDA = auto/CUDA = no/' SuiteSparse_config.mk

# turn off autoconf
sed -i 's/AUTOCC ?= no/AUTOCC = no/' SuiteSparse_config.mk

# makes the UMFPACK be dependent only on AMD
sed -E -i "s/# UMFPACK_CONFIG = -DNCHOLMOD/UMFPACK_CONFIG = -DNCHOLMOD/" SuiteSparse_config.mk

make
cd ../ && cd AMD && make
cd ../ && cd UMFPACK 
# turn off test
sed -i 's/( cd Demo   ; $(MAKE) )/#( cd Demo   ; $(MAKE) )/' Makefile
make -j"${pj}"
## ---------------------------------------------------------------------- ##

echo "umfpack lib: ${umfpack_top}/lib"
echo "umfpack inc: ${umfpack_top}/include"
