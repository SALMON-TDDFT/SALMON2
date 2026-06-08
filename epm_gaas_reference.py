#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone single-machine Python/NumPy/SciPy reference implementation of the
local Empirical Pseudopotential Method (EPM) ground-state solver for GaAs
(Cohen-Bergstresser local pseudopotential, zincblende structure).

This is a *debugging reference* for the Fortran `theory='epm'` module
(src/epm/epm_solver.f90, src/epm/epm_cohen_bergstresser.f90) -- it reproduces
the same lattice/basis/Hamiltonian/momentum-matrix construction, with no MPI
and no OpenMP (plain NumPy/SciPy on a single machine), and writes exactly the
same SYSNAME_k.data / SYSNAME_eigen.data / SYSNAME_tm.data files consumed by
gs_info_ssbe -- so its output can be diffed against the Fortran output, or fed
directly into an SBE real-time run for cross-checking.

Monolithic, self-contained: every parameter (lattice constant, plane-wave
cutoff, k-grid, number of bands/electrons, sysname, output directory) is a
hardcoded constant below -- edit them directly to change the run. No external
input files, no command-line arguments.

Run:
    python3 epm_gaas_reference.py
"""

import numpy as np
from numpy.linalg import norm
from scipy.linalg import eigh

# =============================================================================
# Hardcoded run parameters (edit here -- mirrors a GaAs &epm / &system input)
# =============================================================================
MATERIAL            = 'GaAs'
SYSNAME             = 'GaAs'
OUTPUT_DIR          = './'          # written as SYSNAME_k.data etc. in this dir

A_LATTICE_AU        = 10.68                          # lattice constant a [Bohr]
PW_CUTOFF_RY        = 11.1                          # |G|^2 cutoff [(2*pi/a)^2 units == Ry-style]

NUM_KGRID           = (4, 4, 4)     # Monkhorst-Pack grid n1 x n2 x n3 (uniform weights, no symmetry)

NSTATE              = 8             # number of bands to keep (must match the &system/&sbe input)
NELEC               = 8             # number of valence electrons (closed shell -> NELEC/2 occupied bands)

# =============================================================================
# Cohen-Bergstresser (1966) local pseudopotential form factors for GaAs,
# tabulated at G^2 = 3, 4, 8, 11 in units of (2*pi/a)^2, given in Rydberg and
# converted here to Hartree (SALMON's internal convention) by dividing by 2.
# =============================================================================
RY_TO_HA = 0.5

_CB_FORM_FACTORS_RY = {
    # G^2 : (V^S [Ry], V^A [Ry])
    3:  (-0.23,  0.07),
    4:  ( 0.00,  0.05),   # V^S(4) = 0: structure factor vanishes for zincblende
    8:  ( 0.01,  0.00),
    11: ( 0.06,  0.01),
}


def form_factors(material, G2):
    """Return (V^S, V^A) in Hartree for integer |G|^2 in units of (2*pi/a)^2."""
    if material != 'GaAs':
        raise ValueError("Only 'GaAs' Cohen-Bergstresser form factors are tabulated")
    vs_ry, va_ry = _CB_FORM_FACTORS_RY.get(G2, (0.0, 0.0))
    return vs_ry * RY_TO_HA, va_ry * RY_TO_HA


def lattice_vectors_fcc(a):
    """Conventional fcc primitive vectors for zincblende (Cohen-Bergstresser convention)."""
    h = 0.5 * a
    a1 = np.array([0.0, h,   h])
    a2 = np.array([h,   0.0, h])
    a3 = np.array([h,   h,   0.0])
    return a1, a2, a3


def tau_zincblende(a):
    """Internal two-atom-basis displacement, origin midway between cation/anion."""
    return np.array([a / 8.0, a / 8.0, a / 8.0])


def reciprocal_lattice(a1, a2, a3):
    """b_i = 2*pi (a_j x a_k)/V,  V = a1.(a2 x a3). Returns (b1,b2,b3,volume)."""
    a23 = np.cross(a2, a3)
    a31 = np.cross(a3, a1)
    a12 = np.cross(a1, a2)
    volume = np.dot(a1, a23)
    b1 = (2.0 * np.pi / volume) * a23
    b2 = (2.0 * np.pi / volume) * a31
    b3 = (2.0 * np.pi / volume) * a12
    return b1, b2, b3, volume


def build_plane_wave_basis(a_lattice, b_matrix, cutoff_ry):
    """
    Fixed (k-independent) set of reciprocal lattice vectors
    G = m1*b1 + m2*b2 + m3*b3 with |G|^2 <= cutoff [a.u.^2] (same convention
    as the Fortran build_plane_wave_basis: the cutoff bounds |G|^2 directly).
    Returns (Gcart[npw,3], G2_units[npw] integer |G|^2 in (2*pi/a)^2 units).
    """
    bnorm = norm(b_matrix, axis=1)
    gcut2 = cutoff_ry
    nmax = np.ceil(np.sqrt(gcut2) / bnorm).astype(int) + 1

    g2cart_to_units = (a_lattice / (2.0 * np.pi)) ** 2

    Gcart_list = []
    G2_list = []
    for m1 in range(-nmax[0], nmax[0] + 1):
        for m2 in range(-nmax[1], nmax[1] + 1):
            for m3 in range(-nmax[2], nmax[2] + 1):
                Gtmp = m1 * b_matrix[0] + m2 * b_matrix[1] + m3 * b_matrix[2]
                g2_units = np.dot(Gtmp, Gtmp)
                if g2_units <= gcut2 + 1.0e-8:
                    Gcart_list.append(Gtmp)
                    G2_list.append(int(round(g2_units * g2cart_to_units)))

    Gcart = np.array(Gcart_list)
    G2 = np.array(G2_list, dtype=int)
    return Gcart, G2


def monkhorst_pack_grid(b_matrix, num_kgrid):
    """Uniform MP grid, equal weights 1/nk, no symmetry reduction (matches Fortran ordering)."""
    n1, n2, n3 = num_kgrid
    nk = n1 * n2 * n3
    kpoint = np.zeros((nk, 3))
    kweight = np.full(nk, 1.0 / nk)
    ik = 0
    for i1 in range(1, n1 + 1):
        for i2 in range(1, n2 + 1):
            for i3 in range(1, n3 + 1):
                f1 = (2 * i1 - n1 - 1) / (2.0 * n1)
                f2 = (2 * i2 - n2 - 1) / (2.0 * n2)
                f3 = (2 * i3 - n3 - 1) / (2.0 * n3)
                kpoint[ik] = f1 * b_matrix[0] + f2 * b_matrix[1] + f3 * b_matrix[2]
                ik += 1
    return kpoint, kweight


def build_hamiltonian(material, kvec, Gcart, a_lattice, tau):
    """
    H_{G,G'}(k) = (1/2)|k+G|^2 delta_{G,G'}
                + V^S(|G-G'|^2) cos((G-G').tau) + i V^A(|G-G'|^2) sin((G-G').tau)
    """
    npw = Gcart.shape[0]
    H = np.zeros((npw, npw), dtype=complex)
    units = (a_lattice / (2.0 * np.pi)) ** 2

    kpg = kvec[None, :] + Gcart                      # (npw,3)
    diag = 0.5 * np.einsum('ij,ij->i', kpg, kpg)     # (1/2)|k+G|^2
    np.fill_diagonal(H, diag)

    for i in range(npw):
        for j in range(i + 1, npw):
            dG = Gcart[i] - Gcart[j]
            dG2 = int(round(np.dot(dG, dG) * units))
            VS, VA = form_factors(material, dG2)
            if VS == 0.0 and VA == 0.0:
                continue
            phase = np.dot(dG, tau)
            val = complex(VS * np.cos(phase), VA * np.sin(phase))
            H[i, j] = val
            H[j, i] = np.conj(val)   # Hermitian: H(j,i) uses dG' = -dG -> phase -> -phase, VS even/VA odd
    return H


def momentum_matrix(kvec, Gcart, evec):
    """
    p_mn(k) = sum_G conjg(c_m(G)) (k+G) c_n(G)   (diagonal in the plane-wave basis)
    evec: (npw, nb) matrix of the lowest `nb` eigenvectors (columns).
    Returns p_mn[nb,nb,3].
    """
    npw, nb = evec.shape
    kpg = kvec[None, :] + Gcart                      # (npw,3)
    p_mn = np.zeros((nb, nb, 3), dtype=complex)
    for idir in range(3):
        Dc = kpg[:, idir][:, None] * evec             # (npw, nb)
        p_mn[:, :, idir] = evec.conj().T @ Dc
    return p_mn


def main():
    a1, a2, a3 = lattice_vectors_fcc(A_LATTICE_AU)
    tau = tau_zincblende(A_LATTICE_AU)
    b1, b2, b3, volume = reciprocal_lattice(a1, a2, a3)
    b_matrix = np.array([b1, b2, b3])

    Gcart, G2 = build_plane_wave_basis(A_LATTICE_AU, b_matrix, PW_CUTOFF_RY)
    npw = Gcart.shape[0]

    kpoint, kweight = monkhorst_pack_grid(b_matrix, NUM_KGRID)
    nk = kpoint.shape[0]

    nb = NSTATE
    nocc = NELEC // 2

    print('# EPM (local Cohen-Bergstresser pseudopotential) -- Python/NumPy/SciPy reference')
    print(f'#   material           = {MATERIAL}')
    print(f'#   lattice constant a = {A_LATTICE_AU:.5e} a.u.')
    print(f'#   plane waves        = {npw}')
    print(f'#   k-points           = {nk}')
    print(f'#   bands requested    = {nb} / valence electrons = {NELEC}')

    eigen = np.zeros((nb, nk))
    occup = np.zeros((nb, nk))
    p_tm = np.zeros((nb, nb, 3, nk), dtype=complex)

    for ik in range(nk):
        H = build_hamiltonian(MATERIAL, kpoint[ik], Gcart, A_LATTICE_AU, tau)
        evals, evecs = eigh(H)              # ascending order, like ZHEEV

        eigen[:, ik] = evals[:nb]
        occup[:nocc, ik] = 2.0
        occup[nocc:nb, ik] = 0.0

        p_mn = momentum_matrix(kpoint[ik], Gcart, evecs[:, :nb])
        p_tm[:, :, :, ik] = p_mn

        if (ik + 1) % max(1, nk // 10) == 0 or ik == nk - 1:
            print(f'#   ... diagonalized k-point {ik + 1}/{nk}')

    write_epm_files(SYSNAME, OUTPUT_DIR, MATERIAL, kpoint, kweight, eigen, occup, p_tm)


# =============================================================================
# Output writers -- byte-for-byte compatible with gs_info_ssbe::read_k_data /
# read_eigen_data / read_tm_data (free-format read(fh,*), fixed header-line
# counts and fixed numeric-field counts per data line; see write_epm_files in
# src/epm/epm_solver.f90, which this mirrors line-for-line).
# =============================================================================
def write_epm_files(sysname, outdir, material, kpoint, kweight, eigen, occup, p_tm):
    nk = kpoint.shape[0]
    nb = eigen.shape[0]

    # --- SYSNAME_k.data --- (5 header lines, then "ik kx ky kz weight")
    with open(f'{outdir}{sysname}_k.data', 'w') as f:
        f.write('# k-point data\n')
        f.write('# generated by EPM (Cohen-Bergstresser local pseudopotential) -- Python reference\n')
        f.write(f'# material = {material}, nk = {nk}\n')
        f.write('# units: kx,ky,kz [a.u.], weight (sums to 1)\n')
        f.write('# ik, kx, ky, kz, weight\n')
        for ik in range(nk):
            f.write('{:6d}{:18.10E}{:18.10E}{:18.10E}{:18.10E}\n'.format(
                ik + 1, kpoint[ik, 0], kpoint[ik, 1], kpoint[ik, 2], kweight[ik]))

    # --- SYSNAME_eigen.data --- (3 header lines; per ik: 1 header + nb lines)
    with open(f'{outdir}{sysname}_eigen.data', 'w') as f:
        f.write('# eigenvalue data\n')
        f.write('# generated by EPM (Cohen-Bergstresser local pseudopotential) -- Python reference\n')
        f.write(f'# nk = {nk:6d}, nb = {nb:6d}\n')
        for ik in range(nk):
            f.write(f'# ik = {ik + 1:6d}\n')
            for ib in range(nb):
                f.write('{:6d}{:18.10E}{:18.10E}\n'.format(ib + 1, eigen[ib, ik], occup[ib, ik]))

    # --- SYSNAME_tm.data --- (3 header lines; block 1 p_tm; 1 header; block 2 rvnl_tm = 0)
    with open(f'{outdir}{sysname}_tm.data', 'w') as f:
        f.write('# transition matrix data\n')
        f.write('# generated by EPM (Cohen-Bergstresser local pseudopotential) -- Python reference\n')
        f.write('# block 1: p_tm = <u_m|p|u_n>  (ik, ib, jb, Re px, Im px, Re py, Im py, Re pz, Im pz)\n')
        for ik in range(nk):
            for ib in range(nb):
                for jb in range(nb):
                    px, py, pz = p_tm[ib, jb, 0, ik], p_tm[ib, jb, 1, ik], p_tm[ib, jb, 2, ik]
                    f.write('{:6d}{:6d}{:6d}{:18.10E}{:18.10E}{:18.10E}{:18.10E}{:18.10E}{:18.10E}\n'.format(
                        ik + 1, ib + 1, jb + 1,
                        px.real, px.imag, py.real, py.imag, pz.real, pz.imag))
        f.write('# block 2: rvnl_tm = -i[r,Vnl]  (all zero: local pseudopotential, no nonlocal correction)\n')
        zeros9 = '{:6d}{:6d}{:6d}' + '{:18.10E}' * 6 + '\n'
        for ik in range(nk):
            for ib in range(nb):
                for jb in range(nb):
                    f.write(zeros9.format(ik + 1, ib + 1, jb + 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

    print(f'# EPM (Python reference): wrote ground-state data files for sysname = {sysname}')


if __name__ == '__main__':
    main()
