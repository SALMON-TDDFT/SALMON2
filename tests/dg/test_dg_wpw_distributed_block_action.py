#!/usr/bin/env python3
"""Independent oracle for the distributed sparse WPW block action."""

import numpy as np


def main():
    rng = np.random.default_rng(20260714)
    nw, npw, nvec, nrank = 7, 8, 3, 3
    w_owner = np.arange(nw) % nrank
    p_owner = np.arange(npw) % nrank

    ww = rng.normal(size=(nw, nw)) + 1j * rng.normal(size=(nw, nw))
    ww = 0.5 * (ww + ww.conj().T)
    wp = rng.normal(size=(nw, npw)) + 1j * rng.normal(size=(nw, npw))
    pp = rng.normal(size=(npw, npw)) + 1j * rng.normal(size=(npw, npw))
    pp = 0.5 * (pp + pp.conj().T)

    # Window support makes WP and PP sparse.  The masks are symmetric where
    # required and deliberately include cross-rank support neighbors.
    wp *= np.fromfunction(lambda i, j: ((2 * i + j) % 4) != 0, (nw, npw))
    pp_mask = np.fromfunction(lambda i, j: (abs(i - j) <= 2), (npw, npw))
    pp *= pp_mask

    xw = rng.normal(size=(nw, nvec)) + 1j * rng.normal(size=(nw, nvec))
    xp = rng.normal(size=(npw, nvec)) + 1j * rng.normal(size=(npw, nvec))
    dense_w = ww @ xw + wp @ xp
    dense_p = wp.conj().T @ xw + pp @ xp

    # WW and PP rows live on their output-row owner.  WP entries live on the
    # PW-column owner: it computes owned Y_P directly and sends only partial
    # Y_W sums back to the W-row owners.
    dist_w = np.zeros_like(dense_w)
    dist_p = np.zeros_like(dense_p)
    for rank in range(nrank):
        wrows = np.flatnonzero(w_owner == rank)
        prows = np.flatnonzero(p_owner == rank)
        dist_w[wrows] += ww[wrows, :] @ xw
        dist_p[prows] += pp[prows, :] @ xp
        dist_p[prows] += wp[:, prows].conj().T @ xw
        dist_w += wp[:, prows] @ xp[prows, :]

    np.testing.assert_allclose(dist_w, dense_w, rtol=2e-14, atol=2e-14)
    np.testing.assert_allclose(dist_p, dense_p, rtol=2e-14, atol=2e-14)

    # The same ownership algebra applies to S.  Exercise the cleaned-basis
    # WW identity route with nontrivial sparse WP and PP overlap blocks.
    swp = 0.03 * wp
    spp = np.eye(npw, dtype=complex) + 0.01 * pp
    dense_sw = xw + swp @ xp
    dense_sp = swp.conj().T @ xw + spp @ xp
    dist_sw = xw.copy()
    dist_sp = np.zeros_like(dense_sp)
    for rank in range(nrank):
        prows = np.flatnonzero(p_owner == rank)
        dist_sp[prows] += spp[prows, :] @ xp
        dist_sp[prows] += swp[:, prows].conj().T @ xw
        dist_sw += swp[:, prows] @ xp[prows, :]
    np.testing.assert_allclose(dist_sw, dense_sw, rtol=2e-14, atol=2e-14)
    np.testing.assert_allclose(dist_sp, dense_sp, rtol=2e-14, atol=2e-14)

    # No PW owner stores all PW columns, and PP touches only support neighbors.
    assert max(np.count_nonzero(p_owner == rank) for rank in range(nrank)) < npw
    assert np.count_nonzero(pp) < npw * npw
    print("PASS distributed sparse WW/WP/PW block-action oracle")


if __name__ == "__main__":
    main()
