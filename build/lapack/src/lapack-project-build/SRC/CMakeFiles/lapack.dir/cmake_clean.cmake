file(REMOVE_RECURSE
  "../lib/liblapack.a"
  "../lib/liblapack.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/lapack.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
