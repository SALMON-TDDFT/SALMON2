if (CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
  execute_process(
    COMMAND uname -m
    OUTPUT_VARIABLE SALMON_MACOS_MACHINE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(
    COMMAND brew --prefix
    OUTPUT_VARIABLE SALMON_BREW_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )

  if (NOT SALMON_BREW_PREFIX)
    if (SALMON_MACOS_MACHINE STREQUAL "arm64")
      set(SALMON_BREW_PREFIX "/opt/homebrew")
    elseif (SALMON_MACOS_MACHINE STREQUAL "x86_64")
      set(SALMON_BREW_PREFIX "/usr/local")
    else ()
      message(FATAL_ERROR "Unsupported macOS CPU architecture: ${SALMON_MACOS_MACHINE}")
    endif ()
  endif ()

  set(SALMON_LLVM_PREFIX "${SALMON_BREW_PREFIX}/opt/llvm")
  set(SALMON_LIBOMP_PREFIX "${SALMON_BREW_PREFIX}/opt/libomp")
  set(SALMON_LLVM_CC "${SALMON_LLVM_PREFIX}/bin/clang")
  set(SALMON_LLVM_CXX "${SALMON_LLVM_PREFIX}/bin/clang++")
  set(SALMON_LIBOMP_LIBRARY "${SALMON_LIBOMP_PREFIX}/lib/libomp.dylib")

  if (NOT EXISTS "${SALMON_LLVM_CC}" OR NOT EXISTS "${SALMON_LLVM_CXX}")
    message(FATAL_ERROR
      "Homebrew LLVM clang/clang++ not found under ${SALMON_LLVM_PREFIX}; "
      "install llvm and avoid Apple /usr/bin/clang for OpenMP builds."
    )
  endif ()
  if (NOT EXISTS "${SALMON_LIBOMP_LIBRARY}")
    message(FATAL_ERROR
      "Homebrew libomp not found under ${SALMON_LIBOMP_PREFIX}; install libomp for OpenMP builds."
    )
  endif ()

  if (DEFINED ENV{CC} AND "$ENV{CC}" STREQUAL "/usr/bin/clang")
    message(FATAL_ERROR "Apple /usr/bin/clang cannot be used for macOS OpenMP builds")
  endif ()
  if (DEFINED ENV{CXX} AND "$ENV{CXX}" STREQUAL "/usr/bin/clang++")
    message(FATAL_ERROR "Apple /usr/bin/clang++ cannot be used for macOS OpenMP builds")
  endif ()

  set(ENV{CC} "${SALMON_LLVM_CC}")
  set(ENV{CXX} "${SALMON_LLVM_CXX}")
  set(ENV{OMPI_CC} "${SALMON_LLVM_CC}")
  set(ENV{OMPI_CXX} "${SALMON_LLVM_CXX}")
  set(CMAKE_C_COMPILER "${SALMON_LLVM_CC}" CACHE FILEPATH "Homebrew LLVM C compiler" FORCE)
  set(CMAKE_CXX_COMPILER "${SALMON_LLVM_CXX}" CACHE FILEPATH "Homebrew LLVM C++ compiler" FORCE)

  set(OpenMP_ROOT "${SALMON_LIBOMP_PREFIX}" CACHE PATH "Homebrew libomp prefix" FORCE)
  set(OpenMP_C_FLAGS "-fopenmp" CACHE STRING "OpenMP C flags" FORCE)
  set(OpenMP_C_LIB_NAMES "omp" CACHE STRING "OpenMP C library names" FORCE)
  set(OpenMP_omp_LIBRARY "${SALMON_LIBOMP_LIBRARY}" CACHE FILEPATH "Homebrew libomp library" FORCE)

  if (NOT "$ENV{CPPFLAGS}" MATCHES "(^| )-I${SALMON_LIBOMP_PREFIX}/include( |$)")
    set(ENV{CPPFLAGS} "-I${SALMON_LIBOMP_PREFIX}/include $ENV{CPPFLAGS}")
  endif ()
  if (NOT "$ENV{LDFLAGS}" MATCHES "(^| )-L${SALMON_LIBOMP_PREFIX}/lib( |$)")
    set(ENV{LDFLAGS} "-L${SALMON_LIBOMP_PREFIX}/lib $ENV{LDFLAGS}")
  endif ()
  if (NOT "$ENV{LIBS}" MATCHES "(^| )-lomp( |$)")
    set(ENV{LIBS} "-lomp $ENV{LIBS}")
  endif ()
  if (NOT "$ENV{CFLAGS}" MATCHES "(^| )-fopenmp( |$)" AND
      NOT CMAKE_C_FLAGS_INIT MATCHES "(^| )-fopenmp( |$)")
    set(CMAKE_C_FLAGS_INIT "-fopenmp ${CMAKE_C_FLAGS_INIT}")
  endif ()
  if (NOT CMAKE_C_FLAGS_INIT MATCHES "(^| )-I${SALMON_LIBOMP_PREFIX}/include( |$)")
    set(CMAKE_C_FLAGS_INIT "-I${SALMON_LIBOMP_PREFIX}/include ${CMAKE_C_FLAGS_INIT}")
  endif ()
  if (NOT "$ENV{LDFLAGS}" MATCHES "(^| )-L${SALMON_LIBOMP_PREFIX}/lib( |$)" AND
      NOT CMAKE_EXE_LINKER_FLAGS_INIT MATCHES "(^| )-L${SALMON_LIBOMP_PREFIX}/lib( |$)")
    set(CMAKE_EXE_LINKER_FLAGS_INIT
        "-L${SALMON_LIBOMP_PREFIX}/lib ${CMAKE_EXE_LINKER_FLAGS_INIT}")
  endif ()
  if (NOT "$ENV{LIBS}" MATCHES "(^| )-lomp( |$)" AND
      NOT CMAKE_EXE_LINKER_FLAGS_INIT MATCHES "(^| )-lomp( |$)")
    set(CMAKE_EXE_LINKER_FLAGS_INIT "-lomp ${CMAKE_EXE_LINKER_FLAGS_INIT}")
  endif ()

  message(STATUS "macOS OpenMP: machine=${SALMON_MACOS_MACHINE}, brew_prefix=${SALMON_BREW_PREFIX}")
  message(STATUS "macOS OpenMP: CC=${SALMON_LLVM_CC}")
  message(STATUS "macOS OpenMP: OpenMP_ROOT=${SALMON_LIBOMP_PREFIX}")
endif ()
