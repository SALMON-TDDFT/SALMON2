/* Copyright 2026 SALMONpy contributors. Apache-2.0.
 * Explicitly link the current host-MPI fast-PR runtime and its versioned header.
 * Run in a new directory containing a full-k RT/GS fixture for this MPI size.
 * This is a composition example, not a converged material calculation.
 */
#include <mpi.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "fast_hse_handles.h"

static void check(int code) {
    if (code) {
        fprintf(stderr, "SALMONpy returned status %d\n", code);
        MPI_Abort(MPI_COMM_WORLD, code);
    }
}
int main(int argc, char **argv) {
    int rank, step, accepted;
    int64_t a, b, fresh;
    double old_a[3] = {0}, old_b[3] = {0};
    double next_a[3], next_b[3], j[3], energy, queried_a[3];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (argc != 2) {
        if (!rank) fprintf(stderr, "Usage: native_fast_hse_handles rt.inp\n");
        MPI_Finalize();
        return 2;
    }
    check(sp_attach((int64_t)MPI_Comm_c2f(MPI_COMM_WORLD)));
    if (sp_fast_hse_handle_protocol() != 1) MPI_Abort(MPI_COMM_WORLD, 2);
    check(sp_fast_hse_handle_init(argv[1], .001, 4, 0, &a, j));
    check(sp_fast_hse_handle_new(a, .002, 4, &b, j));
    for (step = 1; step <= 3; ++step) {
        next_a[0] = next_a[1] = 0;
        next_a[2] = .0001 * sin(.3 * step);
        next_b[0] = .0001 * (cos(.2 * step) - 1);
        next_b[1] = .0001 * sin(.2 * step);
        next_b[2] = 0;
        check(sp_fast_hse_handle_step(a, old_a, next_a, j, &energy));
        check(sp_fast_hse_handle_step(b, old_b, next_b, j, &energy));
        memcpy(old_a, next_a, sizeof old_a);
        memcpy(old_b, next_b, sizeof old_b);
    }
    check(sp_fast_hse_handle_close(a));
    check(sp_fast_hse_handle_new(b, .001, 4, &fresh, j));
    /* A stale token cannot acquire the fresh state occupying the same slot. */
    if (sp_fast_hse_handle_info(a, &accepted, queried_a, j) == 0)
        MPI_Abort(MPI_COMM_WORLD, 3);
    check(sp_fast_hse_handle_info(b, &accepted, queried_a, j));
    if (!rank) printf("B continues at step %d: j = %.12g %.12g %.12g\n",
                      accepted, j[0], j[1], j[2]);
    check(sp_fast_hse_handle_close(fresh));
    check(sp_fast_hse_handle_close(b));
    check(sp_detach());
    MPI_Finalize();
    return 0;
}
