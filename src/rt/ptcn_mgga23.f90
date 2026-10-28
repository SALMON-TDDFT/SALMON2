! Backward-compatible entry point for MGGA23 callers.
module ptcn_mgga23
 use ptcn_native,only:mgga23_ptcn_step=>native_ptcn_step
 implicit none
 public::mgga23_ptcn_step
end module
