find_package(Spglib 2.1 CONFIG REQUIRED)

if (NOT TARGET Spglib::symspg)
  message(FATAL_ERROR "Spglib package does not export required target Spglib::symspg")
endif ()

set(EXTERNAL_LIBS Spglib::symspg ${EXTERNAL_LIBS})
set(HAVE_SPGLIB TRUE)
