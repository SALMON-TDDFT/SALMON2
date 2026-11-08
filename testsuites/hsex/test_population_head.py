"""Integer electron/hole pockets must not select an insulating head."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/hsex'))
from metal_head import population_head


def run():
    energies=np.array([[-2.,-1.,1.,2.],[-2.1,-.9,1.1,2.1]])
    filled=np.array([[2.,2.,0.,0.]]*2)
    result=population_head(filled,energies,2)
    assert result['insulating'] and result['indirect_gap_Ha'] == 1.9
    pockets=filled.copy();pockets[1]=[2.,0.,2.,0.]
    assert not population_head(pockets,energies,2)['insulating']
    partial=filled.copy();partial[1]=[2.,1.9,.1,0.]
    assert not population_head(partial,energies,2)['insulating']
    semimetal=energies.copy();semimetal[1,1:]=[.9,1.1,2.1];semimetal[0,2]=.8
    assert not population_head(filled,semimetal,2)['insulating']
    assert not population_head(filled[:,:2],energies[:,:2],2)['insulating']
    for bad in (np.nan,np.inf):
        invalid=energies.copy();invalid[-1,-1]=bad
        try:population_head(filled,invalid,2)
        except ValueError:pass
        else:raise AssertionError('nonfinite spectrum must reject')
    print('common gapped occupied manifold / integer pockets / missing gap PASS')


if __name__=='__main__':run()
