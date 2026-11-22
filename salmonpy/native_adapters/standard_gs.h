/* SALMONpy GS protocol1. Runtime capabilities are profile-specific. Host MPI, collective main-thread calls. */
#ifndef SALMONPY_STANDARD_GS_H
#define SALMONPY_STANDARD_GS_H
#include <stdint.h>
#include "native_context.h"
#ifdef __cplusplus
extern "C" {
#endif
/* State operations are collective and main-thread/sequential, including info
 * and array reads (token agreement communicates). resources and pure protocol/
 * format queries are local. Caller owns buffers and authenticates disk records;
 * no pointers retained. Output buffers are valid only on successful status0.
 * Metadata declaration is not blanket numerical/allocator/MPI-failure recovery.
 */
int sp_gs_protocol(void);
/* Optional1 confirms Cartesian force output after symmetry averaging. */
int sp_gs_cartesian_force_protocol(void);
/* Immutable relative input and dependencies must be visible on every subgroup rank.
 * Standard profile: create computes SCF; fast import-only profile returns2. Unconverged SCF still succeeds and retains the last iterate.
 * Handles retain numerical/XC/MPI owners, with unique tokens invalidated at close. */
int sp_gs_create(const char *filename,int64_t *token);
/* Direct wavefunction/occupation seed, no mixing/HSE resource copy or GS I/O.
 * Same live attachment/model/layout (positions may change). Caller verifies PP
 * contents; Python authenticates model metadata. steps0 full SCF,positive budget.
 * Source token remains unchanged; target owns independent arrays/resources. */
int sp_gs_guess_protocol(void);
/* Synchronous rank-local reader: bounds[10] one-based inclusive x/y/z starts,
 * ends, band start/end,k start/end; shape[5]=nx,ny,nz,global bands,kpoints.
 * Fill complex128 waves[k,band,z,y,x], float64 occupations[k,band],
 * reduced reciprocal points[k,3], weights[k]. No pointers retained.
 * Caller authenticates declared model/pseudopotentials; Fortran validates
 * extents, finite values and sampling, never truncates complex real-Gamma input.
 */
typedef int (*sp_gs_array_reader)(void*,const int64_t*,const int64_t*,void*,double*,double*,double*);
int sp_gs_array_guess_protocol(void);
int sp_gs_create_from_arrays(const char*,int,sp_gs_array_reader,void*,int64_t*);
int sp_gs_create_from_gs(const char *filename,int64_t source,int steps,int64_t *token);
/* Additive active-SCF1: CPU periodic nonsymmetric scalar PZ/Libxc-PZ/PBE, fast-PR HSE06.
 * HSE needs its explicit compute profile; Gamma orbital GS additionally requires sp_gs_hse_orbital_protocol()==1. Imports, symmetry,
 * SOC, +U, DC and ionic evolution are outside this entry.
 * All calls are collective in the same order. Start reads/initializes once;
 * continue skips initialization/no-mix warmup and retains the original total cap.
 * Snapshot tokens retain complete numerical copies in this live GS owner only;
 * no file/process/topology portability or allocation-failure atomicity claimed.
 * scf_info: paused,iterations,total cap,initializations,executed main iterations,
 * live snapshots. Restore changes no lifetime execution/init counters.
 * Snapshots exclude live MPI/Libxc/FFT/backend/output handles; text output is not rewound.
 * HSE numerical ACE/cache/history/phase and true mixing age are retained.
 * Status11 rejects bad/mixed budgets/snapshot tokens before numerical writes;
 * status15 denotes an unsupported SCF profile. Unconverged final results warn.
 */
int sp_gs_scf_protocol(void);
int sp_gs_start(const char *filename,int steps,int64_t *token);
int sp_gs_continue(int64_t token,int steps);
int sp_gs_scf_info(int64_t token,int values[6]);
int sp_gs_scf_snapshot(int64_t token,int64_t *snapshot);
int sp_gs_scf_restore(int64_t token,int64_t snapshot);
int sp_gs_scf_snapshot_close(int64_t token,int64_t snapshot);
/* Optional candidate1: stage detached numerical payload (mode2/snapshot), then
 * host votes globally before publish (mode3/snapshot0), or discard (mode4/0).
 * Same live owner, one pending payload; pending owners reject evolution/direct
 * restore/snapshot/close. Metadata reads remain available. No resource creation,
 * SCF replay, exchange evaluation or file access. Publication assignments can
 * allocate: fatal OOM/process failure is not a handled rollback guarantee.
 * All ranks agree token/mode/snapshot before writes; mode4 is idempotent. */
/* Rank-local SCF disk v1: same binary/input/topology. mode0 writes a NEW file;
 * mode2 stages a numerical read collectively, then candidate(token,0,3) publishes.
 * Host verifies manifest/SHA256 first. Bootstrap prepares static runtime resources
 * but skips initial orbitals/density, mixing and every SCF iteration. */
/* Optional changed-layout SCF restart. Every rank enters in matching order.
 * k protocol1: ordinary XC; protocol2 also dense fast-PR HSE.
 * orbital protocol1: ordinary XC; protocol2 also dense Gamma orbital HSE. Bounds are old_size consecutive
 * (first,last) pairs in rank order; no spatial/mixed migration. */
int sp_gs_scf_k_repartition_protocol(void);
int sp_gs_scf_k_repartition(int64_t token,const char *directory,int old_size,const int *kbounds);
int sp_gs_scf_o_repartition_protocol(void);
int sp_gs_scf_o_repartition(int64_t token,const char *directory,int old_size,const int *obounds);
/* Optional fast-PR HSE Gamma orbital GS + same-layout SCF disk restart.
 * Runtime/input and full-support CPU requirements are in NATIVE_GS. No orbital
 * GS-to-RT handoff is implied. Protocol2 adds explicit one-rank orbital mode;
 * changed orbital rank count additionally requires orbital-repartition2. */
int sp_gs_hse_orbital_protocol(void);
/* Optional HSE spatial compute1.
 * Explicit full-support CPU Gamma profile, selected on both hosts including1 rank.
 * Spatial-repartition2 additionally maps dense ACE/workspace physical grid rows. */
int sp_gs_hse_spatial_protocol(void);
/* CPU SCF bridge: local rows, replicated FFT scratch. 1=full k, 2=full/reduced k, 3=orbital groups. */
int sp_gs_hse_kspace_protocol(void);
/* Spatial repartition1: ordinary XC;2 also dense full-support Gamma HSE.
 * boxes is 6*old_size integers:
 * (xlo,xhi,ylo,yhi,zlo,zhi) for each old rank's physical interior. Pure r-grid
 * layouts only. Copies accepted interiors/mixing and periodic halos; static
 * new-layout FFT/communication resources and dead scratch remain live. */
int sp_gs_scf_r_repartition_protocol(void);
int sp_gs_scf_r_repartition(int64_t token,const char *directory,int old_size,const int *boxes);
/* Optional mixed ordinary-XC repartition1. boxes has10*old_size integers:
 * (xlo,xhi,ylo,yhi,zlo,zhi,olo,ohi,klo,khi) per old rank. Wavefunction/CG
 * five-dimensional ownership must cover the target exactly. Scalar density,
 * potentials and mixing use one replica per spatial box, never summed twice.
 * Same runtime/physical input/assets; stage then publish. Protocol2 maps dense Gamma orbital/spatial HSE; protocol3 also maps full-k/spatial HSE. */
int sp_gs_scf_m_repartition_protocol(void);
int sp_gs_scf_m_repartition(int64_t token,const char *directory,int old_size,const int *boxes);
/* Read-only native stream representation; kinds[8], bytes[128], no GS needed. */
int sp_gs_scf_state_revision(void); /* semantic revision; not an arbitrary source version */
int sp_gs_scf_format_protocol(void);
int sp_gs_scf_format_signature(int *kinds,int8_t *bytes);
int sp_gs_scf_disk_protocol(void);
int sp_gs_scf_bootstrap(const char *filename,int64_t *token);
int sp_gs_scf_disk(int64_t token,const char *filename,int mode);
int sp_gs_scf_candidate_protocol(void);
int sp_gs_scf_candidate(int64_t token,int64_t snapshot,int mode);
/* Import existing wfn with yn_restart=y/read_gs_restart_data=wfn; no SCF or GS writes.
 * Imported forces are unavailable (array_info/read return11). */
int sp_gs_import(const char *filename,int64_t *token);
/* Owned GS→RT seed. Same subgroup/decomposition; inputs A/c and returned charge j.
 * Model owns a copy and survives source GS close. RT registry manages returned handles. */
int sp_tddft_memory_gs_protocol(void);
int sp_tddft_model_from_gs(int64_t gs,double dt,int capacity,const int decomposition[5],
                          int64_t *model,int64_t *state,double j[3]);
/* stats: iterations, native convergence flag (-1 unknown on import), atoms, bands.
 * Standard values: total Hartree energy, residual, tolerance.
 * Fast import-only: values are undefined placeholders; Python returns None.
 * Fast import does not expose eigenvalues/forces or the optional GS-to-RT symbols. */
int sp_gs_info(int64_t token,int stats[4],double values[3]);
/* q:1 density,2 occupations,3 eigenvalues,4 force,5 orbital,6 SCF history (4 rows: iteration/residual/threshold/PBE phase). Fortran-order arrays;
 * shape/lower use7 slots; ndim selects meaningful slots. Local density/orbitals
 * may include halos; eigen/occupation/force are native replicated arrays. */
int sp_gs_array_info(int64_t token,int quantity,int *ndim,int64_t shape[7],int64_t lower[7]);
int sp_gs_read(int64_t token,int quantity,double *real_part,double *imag_part,int64_t length);
/* close destroys only this owner. Native stale/mixed tokens return10; Python
 * balanced close is idempotent. Caller closes every owner before detaching MPI. */
int sp_gs_close(int64_t token);
/* states,live child communicators,management slot high-water. */
int sp_gs_resources(int counts[3]);
/* Rank-local completed timer intervals; diagnostics, not numerical state. */
int sp_gs_timer_protocol(void);
int sp_gs_timers(int64_t token, double seconds[4]);

/* Current directory protocol3: main-thread enter/leave around one call,
 * no host chdir; not concurrent solver support. No old-version fallback. */
int sp_io_directory_protocol(void);
int sp_io_directory_enter(const char *absolute_directory);
void sp_io_directory_leave(void);
#ifdef __cplusplus
}
#endif
#endif
