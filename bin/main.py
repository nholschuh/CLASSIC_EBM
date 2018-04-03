### CLASSIC_EBM
### Jake Aylmer
###
### ---------------------------------------------------------------------------

from __future__ import division

import sys, os, numpy as np, matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from src import parameters as pm, analytics as an
from src import plotting as pl


def main(xi=0.9, smooth_coalbedo=False):
    
    x = np.arange(0.0, 1.001, 0.01)
    HT = np.zeros(len(x))
    
    print "Using %s-coalbedo..." % ('smoothed' if smooth_coalbedo else 
        'step-function')
    
    for k in xrange(len(x)):
        HT[k] = an.HeatTransport(x[k], xi, smooth_coalbedo=smooth_coalbedo)
    
    fig, ax = pl.PlotHeatTransport(x, HT, xi)
    fig.show()
    
    pass


if __name__ == '__main__':
    pl.SetRCParams()
    main(smooth_coalbedo=('smooth_coalbedo' in sys.argv))
