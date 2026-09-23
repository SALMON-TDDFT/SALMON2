"""Check native HSE snapshot against independent full-kernel reference."""
import argparse,json
import numpy as np
from model import NativeModel,hartree,semilocal_potential
from semilocal import Semilocal
from distance_exchange import DistanceExchange

def verify(path):
 m=NativeModel(path);u=m.psi
 with Semilocal() as xc:v,es=semilocal_potential(m.rho,m.nab,xc,m.dv)
 vh,eh=hartree(m.rho,m.h)
 v_error=float(np.max(abs(v-m.read('vxc',m.shape))))
 # Libxc WPBEH shows ~1e-9 derivative noise under one-ulp density changes.
 # This tolerance is for cross-implementation potential parity, not RT residual acceptance.
 assert v_error<2e-9,('GGA divergence potential',v_error)
 w,_=DistanceExchange(m.shape,m.h,m.k).apply(u,u)
 expected=m.core(u)+(vh+v)*u+.25*w
 action_error=float(np.linalg.norm(expected-m.native_hpsi)/np.linalg.norm(expected))
 ex=.125*m.expectation(u,w)
 total=m.expectation(u,m.core(u))+eh+es+ex+m.native_energies[4]
 energy_error=float(total-m.native_energies[0])
 result=dict(hpsi_relative_error=action_error,vxc_max_error=v_error,total_energy_error_Ha=energy_error,xc_energy_error_Ha=float(es+ex-m.native_energies[3]),hse_exchange_Ha=float(ex),total_energy_Ha=float(total))
 assert abs(result["xc_energy_error_Ha"])<1e-9,result
 assert action_error<2e-9,result
 assert abs(energy_error)<1e-9,result
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('export');a=p.parse_args();print(json.dumps(verify(a.export),indent=2))
