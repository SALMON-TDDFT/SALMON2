"""Synthetic broadband charge-current response; no solver, download or GS."""
import numpy as np
from salmonpy.analysis import (optical_conductivity, dielectric_function,
                               optical_linearity, optical_sum_rule)

def main():
    time=np.arange(128)*0.1
    field=np.zeros(128);field[0]=1
    options=dict(current_convention='charge',field_kind='total',
                 units='atomic_gaussian',window='none',
                 metadata={'model':'synthetic sigma=2','polarization':'x'})
    sigma=optical_conductivity(time,2*field,field,**options)
    epsilon=dielectric_function(sigma)
    scaled=optical_conductivity(time,4*field,2*field,**options)
    print('sigma at first nonzero frequency:',sigma.values[1])
    print('epsilon:',epsilon.values[1])
    print('linearity:',optical_linearity(sigma,scaled,rtol=1e-12))
    print('finite-band diagnostic:',optical_sum_rule(sigma,
        expected_weight=2*sigma.omega[-1],normalization='synthetic constant conductivity over sampled band'))

if __name__=='__main__': main()
