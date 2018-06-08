### CLASSIC_EBM
### Jake Aylmer
###
### Generate a plot showing a typical temperature profile for a specified ice-
### edge (xi) (does not necessarily produce a steady-state solution since it
### uses Q specified in parameters.py).
### ---------------------------------------------------------------------------

from __future__ import division

import sys, os, numpy as np, matplotlib.pyplot as plt
import scipy.special as spec
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from src import parameters as pm, analytics as an
from src import fileIO, plotting as pl


def main(xi=np.sin(70*np.pi/180), smooth_coalbedo=False):
    
    x = np.arange(0.0, 1.001, 0.01)
    Tn = np.zeros(1+pm.nmax//2)
    T = np.zeros(len(x))
    Q = an.Q(xi, smooth_coalbedo=smooth_coalbedo)
    
    print "Using %s-coalbedo..." % ('smoothed' if smooth_coalbedo else 
        'step-function')
    
    for n in xrange(0, pm.nmax+2, 2):
        Tn[n//2] = an.Tn(n, xi, Q=Q, smooth_coalbedo=smooth_coalbedo)
    
    for k in xrange(len(x)):
        for n in xrange(0, pm.nmax+2, 2):
            T[k] += Tn[n//2]*spec.legendre(n)(x[k]) 
    
    fig, ax = pl.PlotTemperature(x, T, xi)
    subdir_name = ('Tprof_xi=%.2f'%xi) + ('_SmoothedCoalbedo'*smooth_coalbedo)
    fileIO.SaveFigures([fig], subdir_name)
    fileIO.SaveFigures([fig], subdir_name, '.svg')
    fig.show()
    
    pass


if __name__ == '__main__':
    pl.SetRCParams()
    main(smooth_coalbedo=('smooth_coalbedo' in sys.argv))
