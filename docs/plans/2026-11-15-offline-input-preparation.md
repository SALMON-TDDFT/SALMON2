# Offline input preparation implementation plan

Goal: support preparing material configuration and explicitly selected local assets
on a permitted preparation host, then verify the transferred set offline.
Architecture: reuse MaterialRegistry/resolve_materials, stage explicitly named
relative files, record original recipe and resolved settings plus SHA256/size/
provenance manifest. No fetching, installation, compilation or job submission.
The registry resolver validates backend-specific references against the staged
root. Opaque custom parameters are never scanned for guessed file names.
Tech stack: Python standard library and existing SALMONpy config API.

Approved by the user's “進めましょう” following the preparation-script proposal.
This first preparation milestone leaves ASE structure adapters and grid-spacing
resolution tracked separately as M2/M3; it does not claim scientific preparation
or a runnable solver job when only configuration has been checked.

1. Add failing tests in salmonpy/tests/test_preparation.py: transfer/relocation,
   config validation, hashes/provenance, missing files, corruption, duplicate paths,
   symlinks/traversal, unexpected inventory, no overwrite and failure cleanup.
2. Implement salmonpy.preparation prepare_inputs/verify_inputs. Explicit assets
   mapping: portable relative destination → local source path. Optional provenance
   keyed by asset, unknown entries rejected. Missing provenance remains unknown.
   Build in temporary sibling; reserve final directory without replacing it.
   Verification never loads backends; optional manifest hash anchors transfer.
3. Add executable example and verify CLI, paired guides, AI catalog; register APIs
   and install new module/docs in wheel. Example uses a tiny synthetic local input,
   not an actual Si calculation or production material.
4. Local focused/full tests; sdist→wheel install in fresh venv; execute sample from
   installed package and relocate/verify with sockets disabled. No Actions.
5. Record evidence and mark preparation milestone only. Do not commit/push unless asked.
