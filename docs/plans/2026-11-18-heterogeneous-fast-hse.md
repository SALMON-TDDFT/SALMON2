# Heterogeneous fast-PR HSE model registry

Goal: explicit host-MPI registry for independent different-input HSE models in one
process, retaining existing default/same-model calls and arithmetic. No new GS/SCF.

Each model owns value-only parsed global settings/unit tables and the symmetry
workspace; numerical/Libxc/MPI owners stay in their existing contexts. Selecting a
state selects its model before native step/info/checkpoint/transaction/close.
Generation-checked model tokens separately own model creation, and close prevents
new states while retaining existing children. Shared model Python leases reuse a
single native attachment. Thread concurrency/general native failure recovery are
outside scope. No conventional HSE or automatic model/schema guessing.

1. Capture/activate value inputs and complete symmetry workspace; fresh parser
   coordinate allocation and local status, preserving old single initialization.
2. Model register/spawn/close C ABI and per-state model selection; retain raw old
   handle entry and checkpoints. Scalar IEEE inquiries for affected native adapter.
3. Explicit Python backend.open_registry/context, registered shared-model leases.
4. Existing full/reduced and Taylor/PT-CN GS fixtures: mixed vs independent controls,
   both current/energy/history, independent close/spawn and direct fresh restore;
   MPI1/2, plus old same-model regression. No newly computed GS or references.
5. Native source and complete application-object audits; local builds only, no CI.
   Source/wheel checks, JP/EN/AI/sample/checklist and selected canonical evidence.
