"""python -m salmonpy.examples.fdtd_dielectric EXECUTABLE SHA256 NEW_DIRECTORY"""
import argparse
from salmonpy.fdtd import FDTD, Dielectric, PlaneWave, SalmonExecutable
from salmonpy.geometry import Grid, Geometry, Shape
from salmonpy.fields import AcosPulse

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('executable');p.add_argument('sha256');p.add_argument('directory')
a=p.parse_args()
geometry=Geometry(Grid((32,32,32),(2,2,2)),(Shape('ellipsoid',(8,8,8)),))
pulse=AcosPulse(amplitude=.001,omega=10.,duration=.05)
calculation=FDTD(geometry,(Dielectric(epsilon=2.),),dt=.001,steps=200,
                 source=PlaneWave(pulse,axis=0,position=-12.))
print(calculation.input_text())
result=SalmonExecutable(a.executable,a.sha256).run(calculation,a.directory)
print(result.observations[0].columns)
print(result.observations[0].values[-1])
