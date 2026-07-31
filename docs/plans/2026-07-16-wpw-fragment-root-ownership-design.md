# WPW Fragment-Root Ownership Remediation Design

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

## Scope

This design repairs the production ownership mismatch discovered while
connecting Task 5B rank-local quadrature to `lcfo_flux`.  It supplements the
approved Task 5-10 remediation design and plan; it does not replace or modify
those intentionally untracked documents.  Tasks 6-10 remain unchanged except
that their distributed layout fingerprint must identify the ownership scheme
defined here.

The production basis remains `windowed_kg`.  Stable P-column identity remains

```text
column_id(K,G_id) = (K-1)*n_G + G_id.
```

The mathematical basis and operator conventions are unchanged:

```text
P_(K,G)(r) = chi_K(r) exp(i G.r) / sqrt(Omega_cell).
```

SIPG remains symmetric with `sigma=10/h`, periodic-H1 PP has no face term,
and every canonical face is integrated exactly once.

## Problem

The provisional production context partitions stable P columns arithmetically
over every rank in `dc%icomm_tot`.  `lcfo_flux`, however, owns volume geometry,
Wannier data, and potential data by fragment.  A rank can therefore evaluate a
fragment-local volume contribution whose P row is owned by a different
arithmetic rank.  Rejecting that candidate loses an operator contribution;
recomputing every fragment on the arithmetic owner would introduce forbidden
global fragment scans and global recomputation.

Production ownership must follow the fragment-local data decomposition rather
than an unrelated all-rank arithmetic partition.

## Communicators and ownership

### Fragment communicator

`dc%icomm_frag` remains the communicator over all ranks assigned to one
fragment.  Its rank `dc%id_frag==0` is the fragment root.  Non-root ranks may
compute their existing distributed local contribution, but they never publish
a production context, sparse block, callback, or checkpoint partition.

### Production communicator

All fragment roots collectively create one production communicator.  The
split key is `dc%i_frag`, not `dc%id_tot`, so production rank order is exactly
fragment ID order.  Initialization validates:

- exactly one root participates for each fragment;
- production communicator size equals `dc%n_frag`;
- production rank is `dc%i_frag-1`;
- fragment IDs are positive and unique;
- `n_frag*n_G` and every stable column ID fit the storage integer kind.

Failure publishes no context and terminates collectively on the containing
SALMON communicator.  The canonical-face scanner itself remains
communication-free.

### P-column ownership

Fragment root K owns all and only

```text
(K-1)*n_G+1 : K*n_G.
```

No global owner array is stored.  Ownership is computed in O(1) from the
stable column ID:

```text
K = (column_id-1)/n_G + 1
owner_rank = K-1
```

The context stores only its contiguous owned P IDs and bounded support P IDs.
Support P IDs are generated from the local fragment and its rank-local window
support records; they are sorted and deduplicated without scanning all
fragments.  Because `chi_K` is a product of three one-dimensional windows,
edge and corner overlap can require a bounded 26-neighbor stencil; window
support is therefore not identified with canonical-face adjacency.  Support
IDs do not imply local operator-row ownership: they permit local quadrature
contributions to be staged for bounded neighbor delivery to the owning root.

### W ownership

Fragment root K owns the W rows whose global IDs occur in fragment K's local
basis.  Support W rows are the sorted union of owned W rows and W rows supplied
by bounded window-support records.  No O(N) W-owner table is introduced.

### Neighbor schedules

The canonical-face schedule contains only physical face records and retains
separate periodic minus/plus records for a two-fragment axis.  A distinct
window-support schedule contains every fragment whose buffered support
intersects the local core, including edge and corner neighbors.

Window support is constructed by periodic structured-fragment coordinate
arithmetic or preexisting decomposition-local records, never by scanning
`1:dc%n_frag`.  Its degree is bounded by the window extent and is at most 26
for the accepted one-neighbor-per-axis support.  Initialization rejects a
buffer requiring a larger unimplemented stencil.  Duplicate periodic images
remain distinct overlap records but map to one fragment owner with explicit
image/shift IDs.  Window normalization deduplicates by fragment K: multiple
periodic images select the applicable periodic shift of one `q_K` and never
enter `sum_K q_K^2` as separate fragments.  Image IDs remain distinct only for
geometry regions, face identity, and message routing.

## Volume quadrature data flow

1. Each rank evaluates only its existing locally assigned grid/orbital work
   for its fragment.  It does not traverse another fragment's core box.
2. Normalized windows are evaluated from the bounded support-fragment set
   formed by the local fragment plus its window-support records.  The evaluator
   receives explicit buffer and transition widths and performs no MPI.
3. Before quadrature, a volume-support halo exchanges remote-owned W values
   and gradients only on receiver core points where the sender's buffered W
   support is nonzero.  Messages carry fragment/image, W IDs, grid-box bounds,
   epoch, and payload checksum.  The receiver never reconstructs remote W by
   scanning or globally recomputing fragments.
4. At a uniquely owned core grid point, the spatial fragment evaluates every
   P row in its bounded support set.  This is required because a normalized
   window owned by neighbor K can overlap the spatial fragment's core.  WP
   uses local plus prepared volume-halo W data; PP uses bounded support P IDs.
5. Contributions are accumulated in deterministic stable-ID order.  Dense
   global WP/PP matrices are never allocated.
6. Existing fragment-distributed partial values are reduced on
   `dc%icomm_frag` to `dc%id_frag==0`.  Reduction buffers are bounded by the
   fragment's owned/support sparse block sizes, not by total fragment count.
7. The spatial fragment root partitions the reduced candidates by their O(1)
   P-row owner.  Locally owned rows remain local.  Non-owned rows are sent
   only to roots in the precomputed window-support set by a dedicated
   candidate-halo exchange.  No global all-to-all, global owner metadata, or
   scanner communication is used.
8. The receiving P owner merges contributions in deterministic
   `(P_row,W_row)` or `(P_row,P_column)` order and rejects an unexpected
   sender, ID outside the declared support relation, duplicate contribution
   record, missing epoch, or incomplete neighbor set.
9. The P owner validates finiteness, ordering, dimensions, ownership, complete
   support-neighbor receipt, and the fragment epoch before publishing owned
   candidates to the sparse builder.
10. A failed rank, incomplete reduction, or incomplete neighbor exchange leaves
   the prior valid context
   untouched and publishes no partial block.

Every halo phase uses a prevalidated fixed neighbor schedule.  All expected
receives are posted before sends, all posted requests are completed, and only
then are payload/status failures reduced.  A rank detecting invalid local data
marks its outbound header invalid but does not return early while a peer can
still be waiting.  Header counts are bounded and validated before payload
allocation.  Epoch/image-derived tags are checked against `MPI_TAG_UB`; tag
reuse is forbidden while an earlier epoch remains active.

The production path must not obtain a missing contribution by traversing all
fragments or by reconstructing a remote fragment globally.  If the existing
fragment-local buffers do not cover a required owned-P/support-W pair, route
initialization fails with the detecting fragment, IDs, and required support
extent.

## Canonical faces and trace ownership

The fragment-root layout does not change the canonical-face convention.
Neighbor discovery is completed before scanning.  The provider/halo layer
prepares the W and P traces for one face epoch; the scanner reads those traces
without communication and accumulates by `(W_id,P_id)` into a temporary face
block.  If its P row is not locally owned, the same bounded candidate-halo
layer delivers it to the P owner after the complete face scan.  The scanner
never sends, receives, or publishes a partial face.

For a periodic axis containing two fragments, the physical minus and plus
interfaces remain distinct face records even though they connect the same
fragment pair.  A face identity therefore includes the periodic image/side,
not only the unordered fragment pair.  Each face is integrated once by its
canonical owner and never mirrored into a second scan.

PP remains periodic H1 and receives volume terms only.  A PP face candidate is
an invariant violation.

## Sparse builder and callbacks

`s_dg_wpw_column_layout` gains an explicit fragment-root ownership kind.  The
arithmetic all-rank layout remains available only to mathematical fixtures if
needed; production initialization rejects it.

The sparse builder validates candidate row ownership using the O(1)
fragment-root rule.  WP candidates are stored on the owner of their P row;
PP candidates are stored on the owner of their first P row.  Candidate IDs
must be strictly ordered and unique.  A non-owned candidate is legal only in
the unpublished staging state before candidate-halo delivery; a non-owned
candidate reaching the sparse builder is rejected rather than silently
filtered.

H/S callbacks operate on the production communicator layout and use the same
layout in DG-DC, checkpoint serialization, and RT.  Before applying local
blocks, the callback/halo layer exchanges only the support P and W coefficient
slices declared by the rank-local sparse blocks.  The neighbor schedule is
derived from bounded fragment support and stable IDs; it contains no global
vector owner table.  Callback output is reduced back to the owning fragment
roots through the reverse bounded schedule before it is exposed.  Quadrature
and canonical-face scanners perform none of this communication.

## Lifecycle and failure behavior

The production context lifecycle is:

```text
uninitialized
  -> communicator_valid
  -> quadrature_epoch_valid
  -> trace_epoch_valid
  -> operator_valid
  -> callbacks_bound
```

Every transition validates its complete candidate state before replacing the
previous state.  Errors move the candidate state to invalid but do not expose
partially filled arrays.  Communicator creation and validation occur in two
phases: all ranks first report fragment-local validity on `dc%icomm_tot`, then
validated roots enter root-only collectives.  This prevents a non-root failure
from leaving roots blocked in a production-communicator collective.
Communicator mismatch, stale halo epoch, incomplete neighbor receipt, missing
or duplicate IDs, nonfinite data, unsupported buffer extent, and callback use
after invalidation are terminal for the explicit production route.

The layout fingerprint includes at least:

- ownership kind and schema version;
- `n_frag`, `n_G`, stable-column convention;
- production rank and owned fragment ID;
- owned/support W and P ID hashes;
- fragment/core/buffer geometry hash;
- window, metric, SIPG, and operator-convention hashes.

DG-DC and RT must load the same fingerprint.  RT neither redistributes to a
different ownership layout nor reselects the metric.

## Testing strategy

All production changes follow RED-GREEN-REFACTOR.

### Communicator and layout fixtures

- two fragments with two ranks per fragment: only two roots publish contexts;
- `dc%id_tot` order different from fragment ID order: production rank still
  equals `K-1`;
- missing root, duplicate fragment ID, wrong communicator size, and overflow
  fail closed;
- owned P IDs equal exactly the stable-ID interval for the owned fragment;
- per-rank metadata remains bounded and contains no global owner array.
- edge/corner overlap builds the bounded window-support schedule without a
  global fragment loop;
- two-fragment periodic images contribute one `q_K` per fragment to window
  normalization while preserving distinct overlap/image records;
- a buffer requiring support beyond the implemented stencil fails closed.

### Volume fixtures

- distributed fragment partials reduce to the dense mathematical oracle;
- remote-owned W values and gradients on a local core are read only from a
  complete current volume-support halo epoch;
- missing, stale, truncated, or nonfinite volume-support payloads fail before
  quadrature publication;
- a neighbor-window contribution integrated on the spatial fragment is routed
  exactly once to its P-row owner and matches the dense oracle;
- a required support contribution is neither dropped nor evaluated globally;
- non-owned output, missing support, duplicate pair, nonfinite value, and
  incomplete epoch publish no block;
- unexpected neighbor sender, missing expected neighbor message, and duplicate
  delivery fail before publication;
- one rank reporting an invalid payload still completes the fixed exchange and
  causes collective fail-closed invalidation without a peer hang;
- fragment count scaling keeps per-root storage dependent only on owned and
  neighbor support sizes.

### Face fixtures

- canonical face contribution matches the dense oracle;
- two-fragment periodic plus/minus faces remain distinct and are each counted
  once;
- stale/missing halo epochs fail before scanner publication;
- a face block with a remote-owned P row is delivered after the complete scan
  and counted once by its owner;
- scanner source contains no MPI or global fragment loop;
- PP face candidates remain structurally absent.

### Integration gates

- `lcfo_flux` initializes the root communicator, builds volume candidates,
  prepares trace halos, scans canonical faces, builds sparse H/S blocks, binds
  callbacks, and invokes bounded matrix-free DG-DC only on the explicit route;
- non-root ranks participate only in required existing fragment collectives;
- the Task 5 operator MPI oracle, matrix-free SCF fixture, full MPI/EigenExa
  build, all WPW contracts, and `git diff --check` pass;
- later checkpoint and RT tests verify the same ownership fingerprint.

## Out of scope

- repartitioning G modes among multiple ranks of one fragment;
- arbitrary dynamic fragment migration;
- global candidate all-to-all as an alternative ownership scheme;
- long Si64 HHG production or spectrum interpretation;
- changes to stable column identity, basis normalization, SIPG convention,
  PP face convention, or the accepted midpoint/observable contracts.
