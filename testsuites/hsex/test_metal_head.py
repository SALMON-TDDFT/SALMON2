"""Finite-q cell average: anisotropy, zero-carrier singularity and guards."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools/hsex'))
from metal_head import metallic_cell, check_eta, metallic_sample_indices
from q_mesh import signed_indices


def run():
    mesh = np.array([12]*3)
    rec = np.array([[.08,0,0],[.01,.07,0],[0,.02,.09]])
    q = signed_indices(mesh)/mesh @ rec
    a = np.array([[2.,.2,0],[.2,3.,.1],[0,.1,4.]])
    scalar = 4*np.pi/(.8 + np.einsum('ni,ij,nj->n',q,a,q))
    matrix = np.array([[1.,.2],[.2,.7]])
    values = scalar[:,None,None]*matrix
    values[0] = 0 # the existing optical q=0 body is NOT a metallic limit
    point, average, report = metallic_cell(values,mesh,rec,mesh,tolerance=.01)
    subset=metallic_sample_indices(mesh,rec)
    sparse_point,sparse_average,_=metallic_cell(values[subset],mesh,rec,mesh,tolerance=.01,indices=subset)
    np.testing.assert_allclose(sparse_point,point,rtol=1e-13,atol=1e-13)
    np.testing.assert_allclose(sparse_average,average,rtol=1e-13,atol=1e-13)
    np.testing.assert_allclose(point,4*np.pi/.8*matrix,rtol=2e-6)
    # Independent product Gauss quadrature of an anisotropic screened Coulomb law.
    from numpy.polynomial.legendre import leggauss
    x,w=leggauss(10)
    pts=np.array(np.meshgrid(x,x,x,indexing='ij')).reshape(3,-1).T/2/mesh @ rec
    weights=np.einsum('i,j,k->ijk',w,w,w).ravel()/8
    ref=np.sum(weights*4*np.pi/(.8+np.einsum('ni,ij,nj->n',pts,a,pts)))
    np.testing.assert_allclose(average,ref*matrix,rtol=2e-6)
    # A coarser screening mesh must integrate the final electronic q cell.
    _, fine_average, _ = metallic_cell(values,mesh,rec,2*mesh,tolerance=.01)
    fine_pts = pts/2
    fine_ref = np.sum(weights*4*np.pi/(.8+np.einsum('ni,ij,nj->n',fine_pts,a,fine_pts)))
    np.testing.assert_allclose(fine_average,fine_ref*matrix,rtol=2e-6)
    assert np.linalg.norm(fine_average-average) > 1e-5
    assert report['nested_relative_error'] < .01
    check_eta(values,values*(1+1e-4),.01)
    def reject(f):
        try:f()
        except ValueError:return
        raise AssertionError('expected rejection')
    reject(lambda:check_eta(values,values*2,.01))
    singular=4*np.pi/np.maximum(np.sum(q*q,axis=1),1e-30)
    reject(lambda:metallic_cell(singular[:,None,None],mesh,rec,mesh,tolerance=.01))
    bad=values.copy();bad[3,0,0]=np.nan
    reject(lambda:metallic_cell(bad,mesh,rec,mesh))
    reject(lambda:metallic_cell(np.ones((8,1,1)),np.array([2]*3),rec,mesh))
    # Full matrix, including odd imaginary wings, must survive averaging.
    z=values.astype(complex); z[:,0,1]+=1j*q[:,0]; z[:,1,0]-=1j*q[:,0]
    _,av,_=metallic_cell(z,mesh,rec,mesh,tolerance=.01)
    np.testing.assert_allclose(av,average,atol=1e-11)
    # Large neutral/local-field modes must not hide the uniform charge head.
    small_mesh=np.array([8]*3)
    small_q=signed_indices(small_mesh)/small_mesh
    q2=np.sum(small_q*small_q,axis=1)
    uniform=np.ones(27)/np.sqrt(27)
    charge=np.outer(uniform,uniform)
    neutral=np.eye(27)-charge
    divergent=np.divide(1.,q2,out=np.zeros_like(q2),where=q2>0)
    masked=100*neutral[None]+divergent[:,None,None]*charge
    reject(lambda:metallic_cell(masked,small_mesh,np.eye(3),small_mesh))
    # A large neutral background alone must not trigger numerical wing failures.
    smooth=10000*neutral[None]+(1/(1+q2))[:,None,None]*charge
    _,_,smooth_report=metallic_cell(smooth,small_mesh,np.eye(3),small_mesh)
    assert smooth_report['nested_channel_relative_errors']['charge'] < .005
    first=10000*neutral[None]+np.ones((2,1,1))*charge
    second=10000*neutral[None]+np.ones((2,1,1))*2*charge
    reject(lambda:check_eta(first,second,.02))
    # Neutral-charge mixing must also be checked independently.
    orthogonal=np.zeros(27);orthogonal[0]=1;orthogonal[1]=-1;orthogonal/=np.sqrt(2)
    wing=np.outer(uniform,orthogonal)+np.outer(orthogonal,uniform)
    first=10000*np.eye(27)[None]+np.ones((2,1,1))*wing
    second=10000*np.eye(27)[None]+np.ones((2,1,1))*2*wing
    reject(lambda:check_eta(first,second,.02))
    # Roundoff in nominally zero wings must not fail a tighter physical tolerance.
    uniform216=np.ones(216)/np.sqrt(216)
    charge216=np.outer(uniform216,uniform216)
    neutral216=np.eye(216)-charge216
    smooth216=10000*neutral216[None]+(1/(1+q2))[:,None,None]*charge216
    _,_,roundoff_report=metallic_cell(smooth216,small_mesh,np.eye(3),small_mesh,tolerance=.001)
    assert roundoff_report['nested_channel_relative_errors']['charge'] < .001
    print('metallic cell quadrature / anisotropy / rejection PASS')

if __name__=='__main__':run()
