# Experimental EXX k-grid reduction

Hidden developer namelist: `exx_kgrid_reduction(3)` in `&functional`, default `1,1,1` (disabled). Not a production accuracy recommendation.

Each target retains sources whose relative k-grid indices are divisible by the factors. All target states remain. CPU HSE06, full uniform mesh, integer occupations, k-only MPI only; checkpoint/restart is rejected. Coset FFTs and periodic kernel folding implement the reduced quadrature. MPI traffic and density products are not reduced by the source-count factor.

Independent action tests on mixed and 4x4x4 meshes, 1--8 MPI ranks: maximum absolute error 3.64e-13. Si 4x4x4 smoke (160 steps), reduction 2,2,2: current relative L2 difference 1.92%, electron-count error 2.0e-8. Default-route short trajectory identical to reference. Long-time accuracy and performance remain unvalidated.
