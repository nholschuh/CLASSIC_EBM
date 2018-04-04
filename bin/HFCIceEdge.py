### CLASSIC_EBM
### Jake Aylmer
###
### Plot the heat flux convergence at the ice edge.
### ---------------------------------------------------------------------------

from __future__ import division

import sys, os, numpy as np, matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from src import plotting as pl, fileIO


def main(smooth_coalbedo=False):
    
    fig1, ax1 = pl.PlotHFCIceEdge(smooth_coalbedo=smooth_coalbedo)
    fig2, ax2 = pl.PlotHFCIceEdge(relative_D=np.array([1.0]),
        smooth_coalbedo=smooth_coalbedo, add_linear_fit=True)
    subdir_name = 'HFCIceEdge' + ('_SmoothedCoalbedo'*smooth_coalbedo)
    fileIO.SaveFigures([fig1, fig2], subdir_name)
    fileIO.SaveFigures([fig1, fig2], subdir_name, '.svg')
    fig1.show(); fig2.show()
    
    pass


if __name__ == '__main__':
    pl.SetRCParams()
    main('smooth_coalbedo' in sys.argv)
