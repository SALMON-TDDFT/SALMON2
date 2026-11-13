"""材料ID格子を作成・保存。長さはbohr、0は真空。"""
import argparse
import numpy as np
from salmonpy.geometry import Grid, Geometry, Shape


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output', help='新規.npzファイル')
    args = p.parse_args()
    geometry = Geometry(Grid((32, 32, 32), (2., 2., 2.)),
                        (Shape('ellipsoid', (8., 8., 8.)),))
    print('PMLと接触する面:', geometry.boundary_contacts(('pml',)*6))
    with open(args.output, 'xb') as f:
        np.savez(f, material_ids=geometry.material_ids(),
                 midplane=geometry.slice(2, 16), spacing_bohr=geometry.grid.spacing)


if __name__ == '__main__':
    main()
