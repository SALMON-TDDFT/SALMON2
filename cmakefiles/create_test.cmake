function(create_test)
  set(MPI_REQUIRED ${ARGV0})
  if (MPI_REQUIRED)
    if (USE_MPI)
    else ()
      return()
    endif ()
  endif ()

  get_filename_component(TEST_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
  set(TEST_PREP  "prep_${TEST_NAME}.sh")
  set(TEST_EXEC  "mpiexec_${TEST_NAME}.sh")

  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/preparation  DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/inputfile    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/verification DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE)

  add_custom_command(OUTPUT  ${TEST_PREP}
                     COMMAND echo '\#! /bin/sh' > ${TEST_PREP}
                     COMMAND echo 'find ./ -mindepth 1 -maxdepth 1 -type l | xargs -r unlink' >> ${TEST_PREP}
                     COMMAND echo 'find ./ -mindepth 1 -maxdepth 1 ! -name inputfile ! -name preparation ! -name verification ! -name \"*.sh\" ! -name CMakeFiles ! -name \"*.cmake\" ! -name Makefile | xargs -r rm -rf' >> ${TEST_PREP}
                     COMMAND echo './preparation' >> ${TEST_PREP}
                     COMMAND chmod +x ${TEST_PREP})

  add_custom_command(OUTPUT  ${TEST_EXEC}
                     COMMAND echo '\#! /bin/sh'                                         >  ${TEST_EXEC}
                     COMMAND echo "${TEST_COMMAND}" '< ./inputfile |& tee ./outputfile' >> ${TEST_EXEC}
                     COMMAND chmod +x ${TEST_EXEC})

  add_custom_target("gen_${TEST_NAME}" ALL DEPENDS ${TEST_PREP} ${TEST_EXEC} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  add_test(NAME "prep_${TEST_NAME}"   COMMAND ${TEST_PREP})
  add_test(NAME "run_${TEST_NAME}"    COMMAND ${TEST_EXEC})
  add_test(NAME "verify_${TEST_NAME}" COMMAND verification)

  set_tests_properties("prep_${TEST_NAME}"   PROPERTIES FIXTURES_SETUP    "setup_${TEST_NAME}")
  set_tests_properties("run_${TEST_NAME}"    PROPERTIES RUN_SERIAL ENABLE)
  set_tests_properties("verify_${TEST_NAME}" PROPERTIES FIXTURES_CLEANUP  "cleanup_${TEST_NAME}")

  set_tests_properties("run_${TEST_NAME}"    PROPERTIES FIXTURES_REQUIRED "setup_${TEST_NAME}")
  set_tests_properties("verify_${TEST_NAME}" PROPERTIES FIXTURES_REQUIRED "setup_${TEST_NAME};run_${TEST_NAME}")
endfunction(create_test)

function(create_mpi_test)
  create_test(TRUE)
endfunction()
