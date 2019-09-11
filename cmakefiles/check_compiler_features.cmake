include(CheckFortranSourceCompiles)
include(CheckFunctionExists)

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


check_function_exists("mkdir" SYSTEM_HAS_POSIX_MKDIR_ROUTINE)
check_function_exists("rmdir" SYSTEM_HAS_POSIX_RMDIR_ROUTINE)
