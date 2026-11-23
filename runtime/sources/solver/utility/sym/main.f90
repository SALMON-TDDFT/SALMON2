program main

  use cif_format_module

  integer,parameter :: unit=5

  call read_atom_cif( unit )

end program main
