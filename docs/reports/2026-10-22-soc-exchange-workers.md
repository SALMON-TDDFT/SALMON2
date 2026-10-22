# SOC exchange workers integrated on latest code

Base: c2b85dff; integration retains the latest KY reader, symmetry and forbidden-eigensolver checks.

- Static CPU SOC HSE06, full q/full support, ordinary ACE, orthogonal symmetry-reduced k layout.
- Explicit nproc_k keeps SCF owners; extra launched ranks compute bounded r/rprime blocks.
- Total exchange ranks may exceed representative and full k counts. No automatic owner-count inference; default nproc_k retains ordinary decomposition.
- Namelist tile divisor16 / rprime groups1 defaults; no historical environment selectors.
- Eligible SOC auto selects spinor pair; explicit reference/legacy overrides remain.
- No production executable/job was replaced. RT helper service is unsupported.

## Verification

GNU/Homebrew MPI and serial HSE builds passed; MPI build uses bounds checking.
Application object and source PZHEEV prohibition audits passed.
Registered tests445/455: 7/7 PASS. Four SCF owners, 4/8/12 total MPI,
nominal8 full k points; extra helpers all execute exchange and stop normally.
Two SCF iterations on tiny six-atom WSe2 are an implementation regression,
not physical SCF convergence or material band validation.

Energy output differences0 eV; eigenvalue differences 3.89e-12 and6.98e-12 Ha.
Max rank RSS 113.75/112.80/117.70 MiB (individual process peaks, not simultaneous node peak).
No robust speedup claim: short single measurements and oversubscription are unsuitable for performance certification.
Explicit reference+helpers, invalid group count, and band unfolding reject before service entry.
Review found/fixed legacy selection conflict and unfolding/Finalize helper lifecycle collision.
No remaining critical/important review findings.

Fujitsu production acceptance remains unverified; the historical compiler/layout-dependent discrepancy is not declared resolved.
SALMON-DOCS update is supplied as docs/SALMON-DOCS-soc-exchange-workers.patch;
no separate manual checkout was available here, so that upstream patch has not been applied/pushed.
