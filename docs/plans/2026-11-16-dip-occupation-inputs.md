# TDDFT electron-count diagnostic ownership

Continue the approved input-isolation design. subdip reads nelec/nelec_spin from the existing immutable system occupation record, with legacy fallback. Keep strict tenfold bounds, spin-count precedence, diagnostics, reduction arithmetic and STOP behavior. Do not remove TDDFT configuration guards or expand molecular work.

Test against original subdip for unpolarized, polarized positive-spin-count and fallback-count cases; poison globals before current calls. Include uninitialized-system fallback and rejection at both strict bounds. Repeat complete TDDFT interleaving/reference regression, source/installed Python tests and sdist-to-wheel validation. No remote CI.
