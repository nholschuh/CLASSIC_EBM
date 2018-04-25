### CLASSIC_EBM
### Jake Aylmer
###
### Generate standard case plots of heat transports and heat flux convergences
### across the globe in steady state.
### ---------------------------------------------------------------------------

from __future__ import division

import sys, os, numpy as np, matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from src import parameters as pm, analytics as an
from src import plotting as pl, fileIO


def main(xi=0.95, smooth_coalbedo=False, latitude_axis=False):
    
    x = np.arange(0.0, 1.001, 0.001)
    HT = np.zeros(len(x))
    HFC = np.zeros(len(x))
    Q = an.Q(xi, smooth_coalbedo=smooth_coalbedo)
    
    print "Using %s-coalbedo..." % ('smoothed' if smooth_coalbedo else 
        'step-function')
    
    for k in xrange(len(x)):
        HT[k] =an.HeatTransport(x[k], xi, Q=Q, smooth_coalbedo=smooth_coalbedo)
        HFC[k] = an.HeatFluxConvergence(
            x[k], xi, Q=Q, smooth_coalbedo=smooth_coalbedo)
    
    fig1, ax1 = pl.PlotHeatTransport(x, HT, xi, latitude_axis)
    fig2, ax2 = pl.PlotHeatFluxConvergence(x, HFC, xi, latitude_axis)
    subdir_name = 'StandardHeatTransports'+('_SmoothCoalbedo'*smooth_coalbedo)
    fileIO.SaveFigures([fig1, fig2], subdir_name)
    fileIO.SaveFigures([fig1, fig2], subdir_name, '.svg')
    fig1.show()
    fig2.show()
    
    pass


if __name__ == '__main__':
    pl.SetRCParams()
    main(smooth_coalbedo=('smooth_coalbedo' in sys.argv),
        latitude_axis=('latitude_axis' in sys.argv))
