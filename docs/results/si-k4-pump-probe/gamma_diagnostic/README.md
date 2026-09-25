# Equilibrium negative-response diagnostic (2026-09-25)

Existing unpumped trajectories only; no new propagation. Fixed alpha=.2, beta=0, gamma=.001/.004. Run `analyze.py` with NumPy and Matplotlib to regenerate metrics.json and diagnostic.png.

For gamma=.004 the negative Im epsilon minimum remains near2.04–2.06 eV with cubic observation windows440/600/720/880 au, with minima−85.81/−118.54/−141.98/−173.38. Exponential diagnostic windows eta=.1/.2/.4 eV also retain negative minima−98.37/−42.33/−9.18. This disfavors a simple cubic-window sidelobe explanation. Growth of peak height with observation duration alone is not evidence of dynamical instability: an undamped spectral line also sharpens this way.

Gamma=.001 has no comparable negative structure in1.7–2.4 eV; exponential-window responses there are positive. This limited band test does not establish global passivity for gamma=.001.

The implemented equation is Axc''+beta Axc'+gamma Axc=alpha J. Gamma restores, beta damps; beta=0 does not relax Axc monotonically to zero. The bare auxiliary frequencies are0.8605 eV and1.7210 eV respectively. The measured discrete recurrence relative residual is3.413e-12 and4.131e-12; thus output traces satisfy that recurrence. This does not validate the complete physical model or all current conventions.

A three-frequency least-squares fit to gamma=.004 gives a low mode near2.068 eV, Axc/J magnitude113.20 and phase177.57 degrees. The equation predicts ratio−112.77 at that frequency. This is consistent with an auxiliary-field-coupled response. The fit leaves25% relative current residual and is diagnostic rather than a precision mode decomposition. Agreement with the equation is not independent proof of physical validity, nor does opposite phase alone explain a negative oscillator strength (other modes also have opposite phase).

The evidence supports investigating a gamma-dependent, nonpassive equilibrium response before interpreting pump-induced changes as exciton dynamics. It does not uniquely isolate gamma from the coupled current feedback, underlying response, or extraction conventions. The previous gamma=.004 pump comparison remains a numerical comparison; physical bleaching/shift claims based on that equilibrium baseline are suspended. No parameter changes, additional trajectories, or automatic restarts are authorized by this diagnostic.
