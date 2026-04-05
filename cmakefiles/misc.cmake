# Miscellaneous CMake utilities
# This file provides common utilities and macros used in SALMON build

# Define custom macros
macro(option_set name doc default)
  option(${name} "${doc}" ${default})
endmacro(option_set)

# Prepend a directory path to each item in a list
macro(list_prepend list_name prefix)
  set(temp_list "")
  foreach(item IN LISTS ${list_name})
    list(APPEND temp_list "${prefix}/${item}")
  endforeach()
  set(${list_name} ${temp_list})
endmacro(list_prepend)
