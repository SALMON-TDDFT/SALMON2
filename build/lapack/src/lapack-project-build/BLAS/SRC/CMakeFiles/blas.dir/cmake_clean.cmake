file(REMOVE_RECURSE
  "../../lib/libblas.a"
  "../../lib/libblas.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/blas.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
