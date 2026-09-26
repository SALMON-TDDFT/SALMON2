# Three-dimensional MLWF support implementation plan

User approved replacing the axial cutoff with the previously proposed 3D periodic distance.
Keep initial centers fixed and use the existing radius input (0 means full).
For the supported orthorhombic grid, wrap each coordinate to its minimum image
and cut by Euclidean distance. Full support requires radius >= half box diagonal.
Any unreliable periodic center component protects the whole WF from truncation.
Save xyz centers with diagnostic format version 2; links format remains version 1.

1. Add a real-source regression with tails in all axes, periodic wrapping,
   radius above half x length, full support, and unreliable transverse centers.
   Run against old code and confirm failure before implementation.
2. Extend moments and centers to xyz and apply the same mask to source and loss.
   Update transport test reader and expected masks for version 2.
3. Build and run serial probes and native MPI regression one job at a time.
   Pause and resume the existing spectrum job only during numerical validation.
4. Document new radius semantics and distinguish all historical axial results.
