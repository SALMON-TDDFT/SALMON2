set(source "${EIGENEXA_SOURCE_DIR}/src/eigen_prd_t8.F")
file(READ "${source}" contents)

set(init_guard "      if ( local_size > 1 ) then\n\n        allocate(u0_z")
set(init_fixed "        allocate(u0_z")
string(REPLACE "${init_guard}" "${init_fixed}" contents "${contents}")

set(init_end "        call CSTAB_round_offset(offset4)\n\n      end if\n\n      allocate (tsave")
set(init_end_fixed "        call CSTAB_round_offset(offset4)\n\n      allocate (tsave")
string(REPLACE "${init_end}" "${init_end_fixed}" contents "${contents}")

set(final_guard "      if ( local_size > 1 ) then\n\n        deallocate(u0_z, v0_z)")
set(final_fixed "        deallocate(u0_z, v0_z)")
string(REPLACE "${final_guard}" "${final_fixed}" contents "${contents}")

set(final_end "        deallocate(u1_z, v1_z)\n\n      end if\n\n      deallocate(tsave)")
set(final_end_fixed "        deallocate(u1_z, v1_z)\n\n      deallocate(tsave)")
string(REPLACE "${final_end}" "${final_end_fixed}" contents "${contents}")

if(contents MATCHES "if \\( local_size > 1 \\) then[\r\n]+[\r\n]+        allocate\\(u0_z")
  message(FATAL_ERROR "Failed to patch EigenExa single-thread initialization")
endif()
file(WRITE "${source}" "${contents}")
