# WPW Candidate Routing Performance Design

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

## Problem

The Si64 2x2x2-fragment production run spends more than 22 hours in
`route_dg_wpw_candidate_halo`. A live stack sample places every sampled main
thread in the canonical-face candidate route. All eight ranks consume one CPU
core continuously and use about 22 GiB RSS in total.

The route materializes millions of face candidates and insertion-sorts the
derived-type records four times: twice for duplicate validation and twice for
publication. The resulting O(M^2) work makes the production case effectively
non-terminating.

## Selected Design

Replace record-moving insertion sort with a stable O(M log M) index mergesort.
The candidate records remain immutable; integer indices are ordered by the
existing WP or PP key. Duplicate validation and aggregation consume the same
precomputed order instead of sorting copies repeatedly.

The ordering keys and numerical accumulation order remain identical to the
current implementation:

- WP: `(wp_p, wp_w, source_fragment, image_id)`
- PP: `(pp_r, pp_c, source_fragment, image_id)`

WP candidates precede PP candidates in WP order and PP candidates precede WP
candidates in PP order, preserving the existing mixed-kind behavior.

## Scope

- Preserve candidate ownership and MPI payload formats.
- Preserve duplicate-image fail-closed validation.
- Preserve deterministic aggregation order and floating-point summation order.
- Add a large reverse-order regression that is impractical for insertion sort.
- Add route-stage candidate counts and wall-clock timing for production logs.
- Do not change window, basis, metric, quadrature, face, or checkpoint physics.

## Alternatives Deferred

A specialized dense face router and hash aggregation can provide O(M) behavior,
but both change more routing logic and require broader physics validation. They
remain follow-up options if O(M log M) sorting is not sufficient.

## Verification

Run the focused serial and MPI candidate-halo tests, canonical-face scanner and
production-consumer contracts, matrix-free SCF tests, the full MPI build, and
`git diff --check`. Then run a reduced Si64 preflight and require the candidate
route to finish with bounded memory before starting another full production GS.
