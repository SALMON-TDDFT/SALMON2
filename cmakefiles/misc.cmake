# add_compile_definitions_if(ENV_VAL, DEF_VAL)
#   If ENV_VAL is defined, this command calls add_definitions command.
macro (add_compile_definitions_if ENV_VAL DEF_VAL)
  if (${ENV_VAL})
    add_compile_definitions(${DEF_VAL} ${ARGN})
  endif ()
endmacro (add_compile_definitions_if)

# option_set(ENV_VAL, HELP_STR, DEFAULT_VAL)
#   If option is not declare, this command initializes the option with ${DEFAULT_VAL}
function (option_set ENV_VAL HELP_STR DEFAULT_VAL)
  if (${ENV_VAL}_DEFAULT)
    option(${ENV_VAL} ${HELP_STR} ${${ENV_VAL}_DEFAULT})
  else ()
    option(${ENV_VAL} ${HELP_STR} ${DEFAULT_VAL})
  endif()
endfunction (option_set)

# list_prepend(VAR, PREFIX)
macro (list_prepend VAR PREFIX)
  set(LISTVAR)
  foreach (STR ${${VAR}})
    list(APPEND LISTVAR "${PREFIX}/${STR}")
  endforeach (STR)
  set(${VAR} "${LISTVAR}")
endmacro (list_prepend)

# check_mpi_compiler(COMPILER_NAME, RESULT)
#   This command checks the prefix of ${COMPILER_NAME} is `mpi`.
#   In almost all of MPI compiler, they takes `mpi` prefix. 
#
#   in:  COMPILER_NAME
#   out: ${RESULT} (TRUE or FALSE)
function (check_mpi_compiler COMPILER_NAME RESULT)
  get_filename_component(COMPILER_NAME ${COMPILER_NAME} NAME)
  string(LENGTH ${COMPILER_NAME} NAME_LEN)
  if (${NAME_LEN} LESS 3)
    set(RET FALSE)
  else ()
    string(SUBSTRING ${COMPILER_NAME} 0 3 COMPILER_HEADER)
    string(TOLOWER ${COMPILER_HEADER} COMPILER_HEADER)
    string(COMPARE EQUAL ${COMPILER_HEADER} "mpi" RET)
  endif ()
  set(${RESULT} ${RET} PARENT_SCOPE)
endfunction (check_mpi_compiler)
