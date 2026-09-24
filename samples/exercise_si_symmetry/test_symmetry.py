"""Small, dependency-free geometry and runtime field-guard checks."""
import itertools
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def operations():
    rows = [[float(x) for x in line.split()] for line in (HERE / 'sym.dat').read_text().splitlines()]
    return [rows[i:i+3] for i in range(0, len(rows), 3)]


def irreducible_count(n):
    # SALMON's default half-shifted mesh, represented exactly in units 1/(2*n).
    grid = set(itertools.product(range(1-n, n, 2), repeat=3))
    count = 0
    while grid:
        k = min(grid)
        orbit = {tuple(int(sum(row[j]*k[j] for j in range(3))) for row in op) for op in operations()}
        assert orbit <= set(itertools.product(range(1-n, n, 2), repeat=3))
        grid -= orbit
        count += 1
    return count


class SymmetryTests(unittest.TestCase):
    def test_diamond_and_field(self):
        atoms = {(0,0,0),(0,2,2),(2,0,2),(2,2,0),(1,1,1),(1,3,3),(3,1,3),(3,3,1)}
        for op in operations():
            self.assertEqual([row[2] for row in op], [0,0,1])
            transformed = {tuple(int(sum(row[j]*a[j] for j in range(3))+4*row[3]) % 4 for row in op) for a in atoms}
            self.assertEqual(transformed, atoms)

    def test_affine_group_closure(self):
        ops = operations()
        self.assertEqual(len(ops), 32)
        for a in ops:
            for b in ops:
                composed = [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] +
                            [(sum(a[i][k]*b[k][3] for k in range(3))+a[i][3]) % 1] for i in range(3)]
                self.assertIn(composed, ops)

    def test_exact_counts(self):
        self.assertEqual(irreducible_count(4), 12)
        self.assertEqual(irreducible_count(16), 576)

    def test_runtime_rejects_wrong_group(self):
        compiler = shutil.which('gfortran-15') or shutil.which('gfortran')
        if compiler is None:
            self.skipTest('Fortran compiler unavailable')
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d/'probe.f90').write_text('''module communication
contains
subroutine comm_get_globalinfo(g,p,n)
integer,intent(out)::g,p,n
g=0;p=0;n=1
end subroutine
logical function comm_is_root(p)
integer,intent(in)::p
comm_is_root=p==0
end function
end module
module salmon_global
character(32)::xc='pz',tdcdft='lrc',theory='tddft_pulse'
end module
''')
            (d/'main.f90').write_text('''program probe
use sym_sub
implicit none
real(8)::a(3,3),b(3,3),atoms(3,8)
integer::i,species(8)
character(12)::mode
a=0d0;b=0d0
do i=1,3
a(i,i)=1d0;b(i,i)=2d0*acos(-1d0)
end do
call read_sw_symmetry('yyy')
call init_sym_sub(a,b)
atoms=reshape([0d0,0d0,0d0,0d0,.5d0,.5d0,.5d0,0d0,.5d0,.5d0,.5d0,0d0, &
.25d0,.25d0,.25d0,.25d0,.75d0,.75d0,.75d0,.25d0,.75d0,.75d0,.75d0,.25d0],[3,8])
species=1
call get_command_argument(1,mode)
if(mode=='wrong_atom')atoms(1,1)=.01d0
if(mode=='wrong_kind')species(1)=2
call symmetry_validate_atoms_cartesian(atoms,species)
end program
''')
            subprocess.run([compiler,'-J',str(d),str(d/'probe.f90'),str(ROOT/'src/symmetry/symmetry.f90'),str(d/'main.f90'),'-o',str(d/'probe')],check=True,cwd=d)
            shutil.copy(HERE/'sym.dat',d/'sym.dat')
            self.assertEqual(subprocess.run([str(d/'probe')],cwd=d,capture_output=True).returncode,0)
            for mode in ('wrong_atom','wrong_kind'):
                result = subprocess.run([str(d/'probe'),mode],cwd=d,capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                self.assertIn('does not preserve atomic positions and species',result.stderr)
            # Same yn_symmetry, but actual inversion group would reverse the field.
            (d/'sym.dat').write_text('1 0 0 0\n0 1 0 0\n0 0 1 0\n-1 0 0 0\n0 -1 0 0\n0 0 -1 0\n')
            result = subprocess.run([str(d/'probe')],cwd=d,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('changes the RT field direction',result.stderr)
            # Closed, integer, unimodular and z-preserving is still insufficient:
            # this shear-reflection is not a kinetic-energy symmetry.
            (d/'sym.dat').write_text('1 0 0 0\n0 1 0 0\n0 0 1 0\n1 2 0 0\n0 -1 0 0\n0 0 1 0\n')
            result = subprocess.run([str(d/'probe')],cwd=d,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('not a Cartesian isometry',result.stderr)
            # Eight coset representatives omit FCC centering translations.
            rows = (HERE/'sym.dat').read_text().splitlines()
            (d/'sym.dat').write_text('\n'.join(rows[i+j] for i in range(0,96,12) for j in range(3))+'\n')
            result = subprocess.run([str(d/'probe')],cwd=d,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('do not form a closed group',result.stderr)


if __name__ == '__main__':
    unittest.main()
