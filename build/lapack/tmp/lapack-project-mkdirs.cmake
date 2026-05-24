# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/workspace/build/lapack/src/lapack-project"
  "/workspace/build/lapack/src/lapack-project-build"
  "/workspace/build/lapack"
  "/workspace/build/lapack/tmp"
  "/workspace/build/lapack/src/lapack-project-stamp"
  "/workspace/build/lapack/src"
  "/workspace/build/lapack/src/lapack-project-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/workspace/build/lapack/src/lapack-project-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/workspace/build/lapack/src/lapack-project-stamp${cfgdir}") # cfgdir has leading slash
endif()
