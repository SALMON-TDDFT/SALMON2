# Check build environment settings
# This file detects and configures the build environment

# Set default compiler flags if not already set
if (NOT CMAKE_Fortran_FLAGS)
  set(CMAKE_Fortran_FLAGS "")
endif()

if (NOT CMAKE_C_FLAGS)
  set(CMAKE_C_FLAGS "")
endif()

if (NOT CMAKE_EXE_LINKER_FLAGS)
  set(CMAKE_EXE_LINKER_FLAGS "")
endif()

# Initialize external flags
set(EXTERNAL_FLAGS "")
set(EXTERNAL_LIBS "")

# Check for MPI
if(USE_MPI)
  find_package(MPI REQUIRED Fortran C)
  list(APPEND EXTERNAL_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
  list(APPEND EXTERNAL_LIBS ${MPI_Fortran_LIBRARIES})
endif()

# Set default optimization level
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
