if (NOT PETSC_SOURCE_DIR)
  message(FATAL_ERROR "PETSC_SOURCE_DIR is required.")
endif ()

set(_petsc_macros "${PETSC_SOURCE_DIR}/include/petscmacros.h")
if (NOT EXISTS "${_petsc_macros}")
  message(FATAL_ERROR "PETSc header not found: ${_petsc_macros}")
endif ()

file(READ "${_petsc_macros}" _petsc_macros_content)

set(_old_condition [=[#elif defined(__GNUC__)
  /* GCC 4.8+, Clang, Intel and other compilers compatible with GCC (-std=c++0x or above) */
  #define PetscUnreachable() __builtin_unreachable()]=])
set(_new_condition [=[#elif defined(__GNUC__) && !defined(__FUJITSU)
  /* GCC 4.8+, Clang, Intel and other compilers compatible with GCC (-std=c++0x or above) */
  #define PetscUnreachable() __builtin_unreachable()]=])

string(FIND "${_petsc_macros_content}" "${_old_condition}" _old_position)
if (NOT _old_position EQUAL -1)
  string(REPLACE "${_old_condition}" "${_new_condition}"
         _petsc_macros_content "${_petsc_macros_content}")
  file(WRITE "${_petsc_macros}" "${_petsc_macros_content}")
  message(STATUS "Patched PetscUnreachable for the Fujitsu compiler")
else ()
  string(FIND "${_petsc_macros_content}" "${_new_condition}" _new_position)
  if (_new_position EQUAL -1)
    message(FATAL_ERROR
      "The expected PetscUnreachable definition was not found in ${_petsc_macros}.")
  endif ()
  message(STATUS "PetscUnreachable Fujitsu patch is already applied")
endif ()
