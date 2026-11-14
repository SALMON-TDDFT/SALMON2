/* Synthetic reference plugin, not a fitted material or a performance kernel.
 * Configuration: two native C doubles {gain, decay}, copied at model_create.
 * All ABI calls are serialized by the host. No fixed model/state-count cap. */
#include "salmonpy_material.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>

typedef struct Model {
    salmonpy_handle id; double gain, decay; size_t live; struct Model *next;
} Model;
typedef struct State {
    salmonpy_handle id; Model *model; double dt,t0,time,A[3],j[3];
    int64_t steps; struct State *next;
} State;
static Model *models=NULL;
static State *states=NULL;
static int64_t next_id=1;
static int32_t fail(int32_t code,char *error,size_t size,const char *message) {
    if(error && size)snprintf(error,size,"%s",message);
    return code;
}
static Model *get_model(salmonpy_handle id) {
    Model *m;for(m=models;m;m=m->next)if(m->id==id)return m;return NULL;
}
static State *get_state(salmonpy_handle id) {
    State *s;for(s=states;s;s=s->next)if(s->id==id)return s;return NULL;
}
static int32_t model_create(const void *config,size_t size,salmonpy_handle *out,char *error,size_t cap) {
    double c[2];Model *m;
    if(!out)return fail(SALMONPY_ARGUMENT,error,cap,"missing output handle");
    *out=0;
    if(!config || size!=sizeof(c))return fail(SALMONPY_ARGUMENT,error,cap,"expected two native doubles: gain,decay");
    memcpy(c,config,sizeof(c));
    if(!isfinite(c[0]) || !isfinite(c[1]))return fail(SALMONPY_NUMERICAL,error,cap,"nonfinite model configuration");
    if(next_id==INT64_MAX)return fail(SALMONPY_MEMORY,error,cap,"handle IDs exhausted");
    m=(Model*)calloc(1,sizeof(*m));
    if(!m)return fail(SALMONPY_MEMORY,error,cap,"model allocation failed");
    m->id=next_id++;m->gain=c[0];m->decay=c[1];m->next=models;models=m;*out=m->id;
    return SALMONPY_OK;
}
static int32_t model_destroy(salmonpy_handle id,char *error,size_t cap) {
    Model **link=&models,*m;
    while(*link && (*link)->id!=id)link=&(*link)->next;
    m=*link;
    if(!m)return fail(SALMONPY_HANDLE,error,cap,"invalid model handle");
    if(m->live)return fail(SALMONPY_BUSY,error,cap,"model has live states");
    *link=m->next;free(m);return SALMONPY_OK;
}
static int32_t state_create(salmonpy_handle id,double dt,double t0,const double A[3],salmonpy_handle *out,char *error,size_t cap) {
    Model *m=get_model(id);State *s;int k;
    if(!out)return fail(SALMONPY_ARGUMENT,error,cap,"missing output handle");
    *out=0;
    if(!m)return fail(SALMONPY_HANDLE,error,cap,"invalid model handle");
    if(!A || !isfinite(dt) || dt<=0 || !isfinite(t0))return fail(SALMONPY_ARGUMENT,error,cap,"invalid initial clock or A");
    for(k=0;k<3;k++)if(!isfinite(A[k]))return fail(SALMONPY_NUMERICAL,error,cap,"nonfinite initial A");
    if(next_id==INT64_MAX || m->live==SIZE_MAX)return fail(SALMONPY_MEMORY,error,cap,"state counters exhausted");
    s=(State*)calloc(1,sizeof(*s));
    if(!s)return fail(SALMONPY_MEMORY,error,cap,"state allocation failed");
    s->id=next_id++;s->model=m;s->dt=dt;s->t0=t0;s->time=t0;
    memcpy(s->A,A,sizeof(s->A));s->next=states;states=s;m->live++;*out=s->id;
    return SALMONPY_OK;
}
static int32_t state_observe(salmonpy_handle id,double *time,double current[3],char *error,size_t cap) {
    State *s=get_state(id);
    if(!s)return fail(SALMONPY_HANDLE,error,cap,"invalid state handle");
    if(!time || !current)return fail(SALMONPY_ARGUMENT,error,cap,"missing output buffers");
    *time=s->time;memcpy(current,s->j,sizeof(s->j));return SALMONPY_OK;
}
static int32_t state_step(salmonpy_handle id,const double A[3],double *time,double current[3],char *error,size_t cap) {
    State *s=get_state(id);double a[3],j[3],clock;int k;
    if(!s)return fail(SALMONPY_HANDLE,error,cap,"invalid state handle");
    if(!A || !time || !current)return fail(SALMONPY_ARGUMENT,error,cap,"missing state buffers");
    if(s->steps==INT64_MAX)return fail(SALMONPY_NUMERICAL,error,cap,"step counter exhausted");
    for(k=0;k<3;k++) {
        if(!isfinite(A[k]))return fail(SALMONPY_NUMERICAL,error,cap,"nonfinite A");
        a[k]=A[k];j[k]=s->model->decay*s->j[k]+s->model->gain*(a[k]-s->A[k]);
        if(!isfinite(j[k]))return fail(SALMONPY_NUMERICAL,error,cap,"current overflow");
    }
    clock=s->t0+(double)(s->steps+1)*s->dt;
    if(!isfinite(clock))return fail(SALMONPY_NUMERICAL,error,cap,"clock overflow");
    /* Every validation precedes state mutation. */
    s->steps++;s->time=clock;memcpy(s->A,a,sizeof(a));memcpy(s->j,j,sizeof(j));
    *time=clock;memcpy(current,j,sizeof(j));return SALMONPY_OK;
}
static int32_t state_destroy(salmonpy_handle id,char *error,size_t cap) {
    State **link=&states,*s;
    while(*link && (*link)->id!=id)link=&(*link)->next;
    s=*link;
    if(!s)return fail(SALMONPY_HANDLE,error,cap,"invalid state handle");
    *link=s->next;s->model->live--;free(s);return SALMONPY_OK;
}
int32_t salmonpy_material_query(int32_t version,size_t size,salmonpy_material_api_v1 *out) {
    salmonpy_material_api_v1 api={SALMONPY_MATERIAL_ABI,sizeof(api),"example-c-memory/1",0,1,
        model_create,model_destroy,state_create,state_observe,state_step,state_destroy};
    if(version!=SALMONPY_MATERIAL_ABI || size<sizeof(api) || !out)return SALMONPY_ABI;
    *out=api;return SALMONPY_OK;
}
