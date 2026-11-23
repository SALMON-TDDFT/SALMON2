# Explicit cell resolution

Extract cell-vector/length/spacing precedence into a value-configured resolver.
Keep explicit al before diagonal vectors before point-count times spacing; retain
nonorthogonal vector norms and spacing arithmetic. The compatibility caller
captures input and preserves historical al normalization. Invalid configurations
return before changing the result. Test precedence, tilted cells, independent
results and invalid retention; native reference/regression and packaging remain
local. This is not full heterogeneous TDDFT: k/pseudo/XC globals remain.
