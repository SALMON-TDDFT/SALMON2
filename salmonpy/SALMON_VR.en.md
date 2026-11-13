# Public SALMON VR desktop bridge

[日本語](SALMON_VR.md)

`salmonpy.salmon_vr` exports paired density-change and Ac=A/c CUBE sequences for
the public Mac application. No solver, ASE or plotting dependency is required.
This is a file bridge, not recording control or live streaming.

```python
from salmonpy.salmon_vr import export_cubes
result = export_cubes(
    ['drho_step0.cube', 'drho_step8.cube'], 'new-view',
    steps=[0, 8], dt_au=0.02, missing_ac='zero',
    # ac_files=['Ac_step0.cube', 'Ac_step8.cube'],
    # application='/path/to/SALMON VR.app',
)
# result.launch()  # explicit macOS launch, only after copying the app
```

Use atomic-unit inputs, density changes relative to the initial state, and Ac=A/c.
No unit conversion, density subtraction or E-to-A integration is performed.
Provide strictly increasing native step IDs in 0..999999. `dt_au` is the native
time step, not the interval between saved frames. Files are paired by list order;
no physical meaning or time is inferred from names. Grid, origin and fixed atoms
must agree. The viewer linearly interpolates gaps; the exporter does not.

Outputs are config.json, Ac_Je_drho_CUBE/itNNNNNN_drho.cube and corresponding
itNNNNNN_Ac_x.cube (component selectable), plus salmonpy-vr.json provenance.
The intentionally misspelled public key `InfomationText` is preserved.
Existing destinations are rejected; failed exports remove only their new folder.

Missing A is an error by default. `missing_ac='zero'` explicitly creates display
placeholders, labeled in the GUI configuration and provenance. These are not
computed zero fields. The public Mac app loaded zero placeholders and displayed density isosurfaces
and timeline updates in a smoke test; this is not a guarantee for arbitrary data.

The public Mac app reads configuration/data beside the .app bundle. Pass
`application` to copy your downloaded app into the new export directory, then
explicitly call `launch()`, or place it there manually and open it with Finder.
OS security decisions remain with the user; no quarantine or security changes
are performed. The manual specifies Apple Silicon and macOS 14.5+. Export works
on other Python platforms; Windows app launch is not implemented or verified.

For an OutputArchive, choose indexed CUBE paths and supply their native steps.
A run without saved density distributions cannot provide a density animation.
Budget memory and disk for CUBE data and the copied app; the viewer loads all
frames. CLI sample: `python -m salmonpy.examples.salmon_vr_export --help`.

References: [public Mac manual](https://salmon-tddft.jp/download/macos_manual.pdf),
[application](https://github.com/SALMON-TDDFT/SALMON-VR),
[Otobe et al., Virtual Reality 30, 52 (2026)](https://doi.org/10.1007/s10055-026-01319-0).

Initial scope is positive, axis-aligned orthogonal grids; skew/rotated axes are rejected. Single-frame display and a nonzero first step are unverified in the viewer. If the public tar lacks executable permissions, follow the official Mac manual.

Before launch, the bridge checks Info.plist, CFBundleExecutable, file existence
and execute permission. If permission is missing, it does not launch and reports
a shell-quoted `chmod u+x` command for the user to review and run. No permissions
or OS security settings are changed automatically.
