# SALMON sources
set(SALMON_SOURCES)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ARTED)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/GCEED)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/parallel)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/io)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/math)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gs)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/rt)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/xc)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/atom)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/poisson)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/maxwell)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/common)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/misc)

# Executable binary
add_executable(${TARGET_NAME} main.f90 ${SALMON_SOURCES})
target_link_libraries(${TARGET_NAME} ${EXTERNAL_LIBS})

# Install binary
install(TARGETS ${TARGET_NAME} DESTINATION "bin")
