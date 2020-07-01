# gnu-without-mpi

TARGET = salmon
FC = mpif90
CC = mpicc
FFLAGS = -O3 -Wall -cpp
CFLAGS = -O3 -Wall -std=c99
OMPFLAG = -fopenmp
MPIFLAG = -DUSE_MPI
SIMD_SET = 

LIBLAPACK = -llapack -lblas
# LIBLAPACK = -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -lpthread -ldl -liomp5 -lm

LIBXC = 
# LIBXC = -DUSE_LIBXC -lxcf90 -Ilibxc_installed_dir/include -Llibxc_installed_dir/lib

CONFIG = \
    -DFORTRAN_COMPILER_HAS_2MB_ALIGNED_ALLOCATION \
    -DSYSTEM_HAS_POSIX \
    -DSYSTEM_HAS_POSIX_STAT \
    -DSYSTEM_HAS_POSIX_ACCESS \
    -DSYSTEM_HAS_POSIX_MKDIR \
    -DSYSTEM_HAS_STDIO_REMOVE \
    -DSYSTEM_HAS_POSIX_NFTW \
    -DSYSTEM_HAS_PATH_MAX_IN_LIMITS_H \
    -DSYSTEM_HAS_PATH_MAX_IN_LINUX_LIMITS_H \
    -DUSE_OPT_DOMAIN_IS_POW2 \
    -DUSE_OPT_ARRAY_PADDING \
    -DUSE_OPT_SOFTWARE_PREFETCH 
#    -DFORTRAN_COMPILER_HAS_MPI_VERSION3 \
#    -DUSE_OPT_EXPLICIT_VECTORIZATION \

MODULE_SWITCH = -J

ifneq (,$(wildcard make.body))
include make.body
else 
include gnu_makefiles/make.body
endif


