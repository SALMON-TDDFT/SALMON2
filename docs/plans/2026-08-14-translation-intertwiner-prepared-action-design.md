# Prepared Translation Intertwiner Design

## Problem

`build_dg_translation_character_intertwining_phase` rebuilds and revalidates every finite-translation element map for every target character.  For the Si64 case this repeats the same character-independent global map collectives for 32 characters and dominates the post-Wannier runtime.

## Selected design

Introduce a prepared translation-action object containing the validated full element maps, canonical catalog metadata, provenance fingerprint, workspace receipt, and construction collective count.  Production constructs it once before the character loop.  Per-character phase construction consumes the prepared action and performs only character-ratio propagation, covariance validation, and phase fingerprinting.

The existing one-shot API remains available and is implemented as a compatibility wrapper around prepare plus apply.  This preserves focused fixtures while ensuring production pays preparation cost once.

## Contracts

- Preparation collectively agrees dimensions, generator orders, element words, product table, character-independent catalog metadata, row ownership, and generator maps.
- Preparation validates identity, permutation, generator order, free action, and full composition exactly once.
- Apply collectively agrees the prepared-action provenance and the reference/target character payload.
- Apply does not reconstruct element maps or gather generator maps.
- Prepared maps are row-independent replicated integer metadata.  Their checked memory receipt is published.
- Numerical phase values and phase fingerprints remain equal to the legacy one-shot result.
- Production releases the prepared action after the final character pair.

## Performance target

Character-independent global generator-map reductions scale with the generator count, not with the number of characters.  For Si64, the repeated `32 x 32 x 32` map construction path is replaced by one preparation plus inexpensive indexed map reads.

## Testing

- RED: prepared APIs are missing and the production route prepares inside no character loop.
- Compare prepared apply against the legacy one-shot API for phase payload, provenance, and workspace validity on MPI 1/2/4/8.
- Assert the preparation collective receipt is independent of the number of target characters.
- Retain malformed catalog, ownership, rank-disagreement, and phase covariance adverse fixtures.
- Run the construction MPI fixture, route checker, full configured build, and `git diff --check`.

## Production scope

Only the intertwining-phase construction changes.  Translation-sector splitting, periodic-phase alignment, Gamma sewing, inverse accumulation, and the currently running Si64 process are not modified.
