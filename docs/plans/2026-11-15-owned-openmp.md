# Owned OpenMP execution selection

Continue approved state extraction: store stencil/current/force scheduling flags
in s_stencil. Preserve legacy auto thresholds and default input capture; numerical
consumers use only the passed stencil flags. Add optional internal constructor
selection (auto/orbital/rgrid, uniform across the three kernels), validated before
allocation. This does not change process thread count or introduce concurrency.

Test explicit and auto flags on independent owners, then full contexts with A
orbital/B rgrid against isolated references, both orders and close/reopen. Keep
unchanged native full-driver comparison, audits, Python/package tests, and synced
docs. Common grid/k/XC/global input guards remain. No Actions.
