### CLASSIC_EBM
### Jake Aylmer
###
### ---------------------------------------------------------------------------

from __future__ import division

import sys, os, numpy as np, matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from src import parameters as pm, analytics as an
from src import fileIO, plotting as pl


def main(smooth_coalbedo=False):
    
    xi = np.arange(0.0, 1.001, 0.01)
    Q = np.zeros(len(xi))
    
    print "Using %s-coalbedo..." % ('smoothed' if smooth_coalbedo else 
        'step-function')
    
    for k in xrange(len(Q)):
        Q[k] = an.Q(xi[k], smooth_coalbedo=smooth_coalbedo)
    
    fig, ax = pl.StabilityPlot(xi, np.array([Q])/pm.Q, np.array([1]))
    subdir_name = 'StandardCase' + ('_SmoothedCoalbedo'*smooth_coalbedo)
    fileIO.SaveFigures([fig], subdir_name)
    fileIO.SaveFigures([fig], subdir_name, '.svg')
    fig.show()
    
    pass


if __name__ == '__main__':
    pl.SetRCParams()
    main(smooth_coalbedo=('smooth_coalbedo' in sys.argv))
