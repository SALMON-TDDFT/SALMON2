# DG PW Energy Cutoff Design

## Goal

Allow DG mixed-basis plane waves to be specified by energy cutoff in the input file while keeping the existing wavevector cutoff input backward compatible.

## Scope

- Add a new `&propagation` input: `pw_energy_cutoff`
- Keep existing `k_cutoff_plane_wave` unchanged as a wavevector cutoff
- Prefer `pw_energy_cutoff` when it is positive
- Convert internally from energy cutoff to wavevector cutoff using
  `k_cutoff = sqrt(2 * Ecut_au)`
- Log both:
  - user-provided `pw_energy_cutoff`
  - effective internal `k_cutoff_plane_wave`

## Behavior

1. Existing inputs that only set `k_cutoff_plane_wave` continue to behave exactly as before.
2. If `pw_energy_cutoff > 0`, the code uses it as the authoritative cutoff.
3. `pw_energy_cutoff` is interpreted in the current input energy unit and converted to a.u. in input processing.
4. The internal plane-wave initialization works only with the effective `k` cutoff.

## Approach Options

### Option 1: Reinterpret `k_cutoff_plane_wave`

Treat the old field as energy.

- Pros: smallest code change
- Cons: breaks old inputs and silently changes meaning

### Option 2: Add `pw_energy_cutoff` and preserve `k_cutoff_plane_wave`

Recommended.

- Pros: backward compatible, explicit, low risk
- Cons: two related inputs need precedence rules

### Option 3: Add a mode flag

Use one numeric field plus a separate mode switch.

- Pros: explicit
- Cons: more input complexity than needed

## Recommended Design

Use Option 2.

- Add `pw_energy_cutoff` to `salmon_global` and `&propagation`
- Default to `-1.0d0` meaning "unused"
- During input broadcast/unit conversion:
  - convert `pw_energy_cutoff` with `uenergy_to_au`
- In `init_plane_wave_basis`:
  - if `pw_energy_cutoff > 0`, set `energy_cutoff_pw = pw_energy_cutoff`
  - else use legacy path derived from `k_cutoff_plane_wave`
- Print both values in startup logs
- Record both in `variables.log`

## Error Handling

- If both are unset/nonpositive, use existing legacy default path
- If `pw_energy_cutoff > 0` and `k_cutoff_plane_wave > 0`, use `pw_energy_cutoff` and log that it overrides the legacy k cutoff

## Verification

- Build succeeds
- Old input with only `k_cutoff_plane_wave` selects same PW count as before
- New input with only `pw_energy_cutoff` produces expected converted `k` cutoff
- Logs show both the energy input and effective internal `k` cutoff
