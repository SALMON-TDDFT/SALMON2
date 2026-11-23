#include "salmonpy_material.h"
#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <limits>
int main() {
    salmonpy_material_api_v1 a;
    assert(salmonpy_material_query(99,sizeof(a),&a)==SALMONPY_ABI);
    assert(salmonpy_material_query(1,1,&a)==SALMONPY_ABI);
    assert(salmonpy_material_query(1,sizeof(a),&a)==SALMONPY_OK);
    assert(a.abi_version==1 && a.struct_size==sizeof(a));
    char error[256];double cfg[2]={2.,.8},zero[3]={0.,0.,0.};
    salmonpy_handle model=0,one=0,two=0;
    assert(a.model_create(cfg,sizeof(cfg),&model,error,sizeof(error))==0);
    cfg[0]=999.; /* Configuration is copied. */
    assert(a.state_create(model,.02,0.,zero,&one,error,sizeof(error))==0);
    assert(a.state_create(model,.02,0.,zero,&two,error,sizeof(error))==0 && one!=two);
    assert(a.model_destroy(model,error,sizeof(error))==SALMONPY_BUSY);
    double time,j[3],last[3]={0.,0.,0.},reference[3]={0.,0.,0.};
    for(int n=1;n<=8;n++) {
        double A[3]={n*.01,n*.02,-n*.005};
        assert(a.state_step(one,A,&time,j,error,sizeof(error))==0);
        for(int k=0;k<3;k++) {
            reference[k]=.8*reference[k]+2.*(A[k]-last[k]);last[k]=A[k];
            assert(std::isfinite(j[k]) && std::abs(j[k]-reference[k])<1e-14);
        }
        assert(time==n*.02);
        assert(a.state_observe(two,&time,j,error,sizeof(error))==0 && time==0.);
        for(double value:j)assert(value==0.);
    }
    double bad[3]={std::numeric_limits<double>::quiet_NaN(),0.,0.};
    assert(a.state_step(one,bad,&time,j,error,sizeof(error))==SALMONPY_NUMERICAL);
    assert(a.state_observe(one,&time,j,error,sizeof(error))==0 && time==8*.02);
    for(int k=0;k<3;k++)assert(j[k]==reference[k]);
    assert(a.state_destroy(one,error,sizeof(error))==0);
    assert(a.state_observe(one,&time,j,error,sizeof(error))==SALMONPY_HANDLE);
    salmonpy_handle three=0;
    assert(a.state_create(model,.02,0.,zero,&three,error,sizeof(error))==0 && three!=one);
    assert(a.state_destroy(two,error,sizeof(error))==0);
    assert(a.state_destroy(three,error,sizeof(error))==0);
    assert(a.model_destroy(model,error,sizeof(error))==0);
    assert(a.model_destroy(model,error,sizeof(error))==SALMONPY_HANDLE);
    std::puts("C++ material ABI, independent handles, ownership and atomic failure PASS");
}
