# SALMON sources
set(SALMON_SOURCES)

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
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/symmetry)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/so)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/plusu)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ms)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ttm)

# Executable binary
add_executable(${TARGET_NAME} main.f90 ${SALMON_SOURCES})
target_link_libraries(${TARGET_NAME} ${EXTERNAL_LIBS})

# Install binary
install(TARGETS ${TARGET_NAME} DESTINATION "bin")
