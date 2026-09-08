# DG Density Owner-Distributed Communication Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** DG 密度再構築の world 通信を fragment subgroup 内へ分散し、`rank 0` 集中の `recv_wait` / `unpack` を削減する

**Architecture:** subgroup 内 reduce は維持し、その後の world 通信だけを subgroup root 固定から subgroup handler 分散へ置き換える。handler はまず `target_world_rank mod isize_frag` の静的規則で決め、sender/receiver が同じ規則を再計算して一致するようにする。

**Tech Stack:** Fortran, MPI wrapper (`comm_isend`, `comm_irecv`, `comm_wait_all`), existing DG RT density maps

---

### Task 1: Isolate Handler Selection

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`
- Test: build via `make -C build -j2 salmon`

**Step 1: Add a small helper for handler ownership**

Add a local function or inline helper logic that computes `handler_id_frag = modulo(target_rank, dg_frag%isize_frag)` and tests `dg_frag%id_frag == handler_id_frag`.

**Step 2: Use the helper only for logging first**

Update trace output so a handler-enabled run can show which subgroup rank owns each world target without changing communication yet.

**Step 3: Run build to verify it compiles**

Run: `make -C build -j2 salmon`  
Expected: `Built target salmon`

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90
git commit -m "refactor: prepare handler selection for density comm"
```

### Task 2: Reassign World Send Buffer Ownership

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`
- Test: build via `make -C build -j2 salmon`

**Step 1: Change `rho_send` allocation gating**

Replace subgroup-root-only `rho_send(irank)` allocation with handler-rank-only allocation based on the new helper.

**Step 2: Change `send_sum -> rho_send` copy ownership**

Only the handler rank for `irank` should extract the packed buffer from `send_sum` and populate `rho_send(irank)`.

**Step 3: Preserve self-owner local accumulation**

Keep the existing local `rho/rho_s` accumulation path for subgroup targets that map back to the handler’s own world rank.

**Step 4: Run build to verify it compiles**

Run: `make -C build -j2 salmon`  
Expected: `Built target salmon`

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90
git commit -m "refactor: distribute density send buffers across subgroup ranks"
```

### Task 3: Reassign World Receive Ownership

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`
- Test: build via `make -C build -j2 salmon`

**Step 1: Change `rho_recv` allocation gating**

Allocate `rho_recv(source_rank)` only on the handler rank responsible for receiving from that source into this world owner.

**Step 2: Change `comm_irecv` / `comm_isend` posting conditions**

Post world receive/send requests only on handler ranks, not on subgroup root.

**Step 3: Keep unpack local to receiving handler**

Unpack using existing `density_recv_map` on the handler rank that actually received the payload.

**Step 4: Run build to verify it compiles**

Run: `make -C build -j2 salmon`  
Expected: `Built target salmon`

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90
git commit -m "refactor: distribute density recv buffers across subgroup ranks"
```

### Task 4: Extend Trace For Load Distribution

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`
- Test: build via `make -C build -j2 salmon`

**Step 1: Add handler-aware communication trace**

Print `rank`, `id_frag`, and summary counts for send/recv ownership so a run can show whether traffic is still concentrated.

**Step 2: Keep existing timing breakdown intact**

Do not remove current `comm breakdown`; extend it only as needed for interpretation.

**Step 3: Run build to verify it compiles**

Run: `make -C build -j2 salmon`  
Expected: `Built target salmon`

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90
git commit -m "chore: add handler-aware density comm traces"
```

### Task 5: Verify Runtime Behavior

**Files:**
- Runtime validation only

**Step 1: Run the target MPI case**

Run the same reproducer used for the captured `output.47962447.zip`.

**Step 2: Extract root traces**

Confirm that:
- `rank 0` no longer dominates `recv_wait`
- multiple handler ranks now own nonzero `recv_wait` / `unpack`
- former sender-only ranks show reduced `send_wait`

**Step 3: Check correctness signals**

Confirm:
- normalization logs remain sane
- all ranks still reach `before-hartree`
- no communicator mismatch or deadlock appears

**Step 4: Commit**

```bash
git add .
git commit -m "test: validate distributed density communication ownership"
```
