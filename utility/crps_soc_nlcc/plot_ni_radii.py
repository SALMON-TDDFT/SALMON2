#!/usr/bin/env python3
"""Plot existing radial candidate tables; no calculation or executable changes."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read(path):
    lines=path.read_text().splitlines();nc=int(lines[3].split()[2]);ch=np.array([list(map(float,x.split())) for x in lines[4:4+nc]])
    rows=np.array([list(map(float,x.split())) for x in lines[4+nc:-1]])
    return ch,rows

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('baseline',type=Path);p.add_argument('selected',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    fig,axes=plt.subplots(5,2,figsize=(11,14),layout='constrained')
    for path,label,color in [(a.baseline,'Initial 2.2 / 2.0 / 2.0 bohr','#e1812c'),(a.selected,'Larger 3.2 / 2.8 / 3.2 bohr','#3274a1')]:
        ch,rows=read(path);r=rows[:,0];h=r[1]
        for k,c in enumerate(ch):
            n=int(round(c[2]/h));v=rows[:n+1,5+2*k]
            axes[k,0].plot(r[:n+1],v,color=color,label=label,lw=2)
            axes[k,0].scatter([r[n]],[v[-1]],color=color,s=22,zorder=3)
            dd=np.diff(v,2)/h**2
            axes[k,1].plot(r[3:n-1],dd[2:-1],color=color,lw=1.6)
            for ax in axes[k]:ax.set_title(f'l={int(c[0])}, j={int(c[1])}/2');ax.grid(alpha=.2);ax.set_xlabel('r (Angstrom)')
    for k in range(5):
        axes[k,0].set_ylabel('Ionic potential (eV)');axes[k,1].set_ylabel('Second derivative (eV / Angstrom²)')
    axes[0,0].legend(fontsize=9)
    fig.suptitle('Ni SOC + NLCC: potential shape inside each channel Rps\nDots mark Rps; s-local in both candidates',fontsize=15)
    fig.savefig(a.output,dpi=150)
