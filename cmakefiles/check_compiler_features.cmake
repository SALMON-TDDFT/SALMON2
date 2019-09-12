include(CheckFortranSourceCompiles)
include(CheckIncludeFile)
include(CheckSymbolExists)

if (USE_MPI)
  find_package(MPI REQUIRED)
  if ("${MPI_Fortran_VERSION_MAJOR}" VERSION_GREATER_EQUAL "3")
    set(FORTRAN_COMPILER_HAS_MPI_VERSION3 ON CACHE STRING "")
    message(STATUS "SALMON will uses the optimization by MPI version 3.")
  endif ()
endif ()


# Can the fortran compiler compiles 2MB aligned memory allocation?
## NOTE: this directive is specified in Fortran 2018 standard.
check_fortran_source_compiles([[
complex(8),allocatable :: zbuf(:,:)
!dir$ attributes align : 2097152 :: zbuf
allocate(zbuf(10,20))
end]]
FORTRAN_COMPILER_HAS_2MB_ALIGNED_ALLOCATION SRC_EXT F90)


set(CMAKE_REQUIRED_DEFINITIONS "-D_XOPEN_SOURCE=500")
check_include_file(unistd.h SYSTEM_HAS_POSIX)
if (SYSTEM_HAS_POSIX)
  check_symbol_exists(stat   "sys/stat.h;sys/types.h;unistd.h" SYSTEM_HAS_POSIX_STAT)
  check_symbol_exists(access "unistd.h"                        SYSTEM_HAS_POSIX_ACCESS)
  check_symbol_exists(mkdir  "sys/stat.h;sys/types.h;unistd.h" SYSTEM_HAS_POSIX_MKDIR)
  check_symbol_exists(nftw   "ftw.h"                           SYSTEM_HAS_POSIX_NFTW)
else ()
  message(FATAL_ERROR "SALMON requires POSIX API for IO access")
endif ()
unset(CMAKE_REQUIRED_DEFINITIONS)

check_symbol_exists(remove "stdio.h" SYSTEM_HAS_STDIO_REMOVE)

check_symbol_exists(PATH_MAX "limits.h"       SYSTEM_HAS_PATH_MAX_IN_LIMITS_H)
check_symbol_exists(PATH_MAX "linux/limits.h" SYSTEM_HAS_PATH_MAX_IN_LINUX_LIMITS_H)
