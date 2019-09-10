include(CheckFortranSourceCompiles)
include(CheckFunctionExists)


# Can the fortran compiler compiles 2MB aligned memory allocation?
## NOTE: this directive is specified in Fortran 2018 standard.
check_fortran_source_compiles([[
complex(8),allocatable :: zbuf(:,:)
!dir$ attributes align : 2097152 :: zbuf
allocate(zbuf(10,20))
end]]
FORTRAN_COMPILER_PASS_2MB_ALIGNED_ALLOCATION SRC_EXT F90)

add_compile_definitions_if(FORTRAN_COMPILER_PASS_2MB_ALIGNED_ALLOCATION  SALMON_ENABLE_2MB_ALIGNED_ALLOCATE)


check_function_exists("mkdir" SYSTEM_HAS_POSIX_MKDIR_ROUTINE)
check_function_exists("rmdir" SYSTEM_HAS_POSIX_RMDIR_ROUTINE)
