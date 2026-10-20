# SALMON: Scalable Ab-initio Light-Matter simulator for Optics and Nanoscience

SALMON is an open-source software based on first-principles time-dependent density functional theory
to describe optical responses and electron dynamics in matters induced by light electromagnetic fields.

SALMON has been tested and optimized to run in the following supercomputer platforms:

- K-computer (It is scheduled to be permanently shut down at the end of August, 2019)
- Fujitsu FX100 supercomputer system
- Fujitsu supercomputer (Fugaku, FX1000, FX700) with A64FX processor
- Linux PC Cluster with Intel Xeon Phi (Knights Landing architecture)
- Linux PC Cluster with x86-64 CPU

For more information, please visit our website.

http://salmon-tddft.jp/

## Fujitsu MPI topology mapping

The CMake option `USE_FJMPI` controls the Fujitsu MPI topology routines used for
Tofu network-oriented process mapping. It defaults to `ON` when MPI is enabled
and the Fujitsu compiler is detected, and to `OFF` otherwise.

When using the Fujitsu compiler with another MPI implementation, pass
`--disable-fjmpi` to `configure.py` (or `-DUSE_FJMPI=OFF` to CMake) to use the
generic process mapping. `--enable-fjmpi` (or `-DUSE_FJMPI=ON`) requires
`USE_MPI=ON` and an MPI installation providing the `mpi_ext` Fortran module and
`FJMPI_Topology_*` routines. For builds using the GNU makefiles, define
`USE_FJMPI` in the Fortran preprocessor flags only when these routines are
available.

## License

SALMON is available under Apache License version 2.0.

    Copyright 2017-2026 SALMON developers

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
