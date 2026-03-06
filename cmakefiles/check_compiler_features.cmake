# Check compiler features and capabilities

# Fortran compiler checks
include(CheckFortranSourceCompiles)

# Check for Fortran features as needed
# (Add specific feature checks here if needed)

# C compiler checks  
include(CheckCSourceCompiles)
include(CheckSymbolExists)

# Check for C features as needed
# (Add specific feature checks here if needed)

# Check for POSIX API support
check_symbol_exists(stat "sys/stat.h" SYSTEM_HAS_POSIX_STAT)
check_symbol_exists(access "unistd.h" SYSTEM_HAS_POSIX_ACCESS)
check_symbol_exists(mkdir "sys/stat.h" SYSTEM_HAS_POSIX_MKDIR)
check_symbol_exists(nftw "ftw.h" SYSTEM_HAS_POSIX_NFTW)

if(SYSTEM_HAS_POSIX_STAT AND SYSTEM_HAS_POSIX_ACCESS AND SYSTEM_HAS_POSIX_MKDIR AND SYSTEM_HAS_POSIX_NFTW)
  set(SYSTEM_HAS_POSIX ON)
endif()

# Check for PATH_MAX location
check_symbol_exists(PATH_MAX "limits.h" SYSTEM_HAS_PATH_MAX_IN_LIMITS_H)

# Check for remove function
check_symbol_exists(remove "stdio.h" SYSTEM_HAS_STDIO_REMOVE)
