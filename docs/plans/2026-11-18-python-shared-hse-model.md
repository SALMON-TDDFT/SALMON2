# Python shared fast-PR model

Use backend.open_model(initial,dt,context) to prepare one verified work bundle and
attach one host communicator. Return NativeFastHSEModel with collective
new_state(point_id=None,dt=None,capacity=None). Each _State uses a token-bound
library facade and its own lease, accepted endpoint, transaction and output label.
Reuse existing numerical/checkpoint/error semantics. Closing model prevents new
states; live states keep their lease and can finish. Detach after model closed and
last state closed. No hidden numerical template or per-state GS directory copies.

Add a separate native model_new capability so zero-live-state models can allocate
fresh contexts without retaining a full hidden state. Preserve live-token new
semantics and old ABI. One immutable non-symmetric model/process remains; same
communicator/main thread. New API does not automatically change MPI point plans.
Output journal inventory is restricted to each state's two declared streams.

Verify missing entry/lease/filter tests first; then unit tests, real MPI1/2
Taylor/learned PT-CN single vs interleaved, close/reopen, save/direct restore and
candidate publication; independent resumed futures, unchanged sibling state and
normal old routes. JP/EN/AI/sample and sdist-to-wheel/installed checks. No Actions.
