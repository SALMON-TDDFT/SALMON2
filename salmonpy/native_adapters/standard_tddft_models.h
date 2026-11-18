/* Copyright 2026 SALMONpy contributors. Apache-2.0. */
#ifndef SALMONPY_STANDARD_TDDFT_MODELS_H
#define SALMONPY_STANDARD_TDDFT_MODELS_H
#include "standard_tddft_handles.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Registry1: one host-MPI attachment, distinct immutable full-k standard models.
 * All ranks enter on their main thread in matching order; same subgroup.
 * Model tokens are opaque NEGATIVE generation IDs; state tokens stay positive.
 * Register initializes and returns a first state. Each model owns dt/capacity/A0,
 * decomposition, parser inputs and derived/XC scratch. Only existing scalar
 * periodic semilocal profile; no other native initializer/thread concurrency.
 * Host supplies correct model workdir for register/spawn/state calls and immutable
 * GS/PP/input. Caller authenticates checkpoints/inputs. j is charge-current density.
 * Do not mix legacy init/model_new with registry mode; use explicit model_spawn.
 */
int sp_tddft_registry_protocol(void);
int sp_tddft_model_register(const char *input,double dt,int capacity,const double A0[3],
                           const int decomposition[5],int64_t *model_token,
                           int64_t *state_token,double j[3]);
int sp_tddft_model_spawn(int64_t model_token,int64_t *state_token,double j[3]);
/* Retire creation token. Existing states remain usable and release owner at last close. */
int sp_tddft_model_close(int64_t model_token);
#ifdef __cplusplus
}
#endif
#endif
