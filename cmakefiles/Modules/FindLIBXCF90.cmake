#
#  Copyright 2018-2020 SALMON developers
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#-----------------------------------------------------------------------------------------

# FindLIBXCF90.cmake
# ----------
# This module finds a installed binary and Fortran90 includes of 
# Libxc (http://www.tddft.org/programs/libxc/) package.
#
# This module sets the following variables:
#  LIBXCF90_INCLUDE_DIRS
#  LIBXCF90_LIBRARIES
#  LIBXCF90_DEFINITIONS
#  LIBXCF90_FOUND
#  LIBXCF90_VERSION

find_package(PkgConfig)

# Detect libxc setting
pkg_check_modules(PC_LIBXCF90 QUIET libxcf90)

set(_libxc_pathes
  ~/program            # Uemoto's environment :P
  ~/.local
  /usr/local/Celler    # Homebrew
  /opt/local           # MacPorts
  ~/Library/Frameworks # MacOS
  /Library/Frameworks  #
  /sw                  # Flink
  /usr                 # General linux package manager
  /opt                 #
)


# Search libxc libxcf90 (fortran90 library)
find_library(
  LIBXCF90_LIBRARY_XCF90
  NAMES xcf90
  PATHS ${PC_LIBXCF90_LIBRARY_DIRS} ${_libxc_paths}
)

# Search libxc (general libxc library)
find_library(
  LIBXCF90_LIBRARY_XC
  NAMES xc 
  PATHS ${PC_LIBXCF90_LIBRARY_DIRS} ${_libxc_paths}
)

set(LIBXCF90_LIBRARIES ${LIBXCF90_LIBRARY_XCF90} ${LIBXCF90_LIBRARY_XC})



# Search 'include' directory containing header file
find_path(
  LIBXCF90_INCLUDE_DIRS
  NAMES xc_version.h
  PATHS ${PC_LIBXCF90_INCLUDE_DIRS} ${_libxc_paths} 
)

file(STRINGS "${LIBXCF90_INCLUDE_DIRS}/xc_version.h" XC_VERSION_H REGEX "^#define XC_VERSION \"[^\"]*\"$")
string(REGEX REPLACE "^.*XC_VERSION \"([0-9\.]+).*$" "\\1" LIBXCF90_VERSION "${XC_VERSION_H}")

# Show error messages
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  LIBXCF90 DEFAULT_MSG
  LIBXCF90_LIBRARIES LIBXCF90_INCLUDE_DIRS
)

# Store results
mark_as_advanced(LIBXCF90_INCLUDE_DIRS LIBXCF90_LIBRARIES LIBXCF90_VERSION)

