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
    
    fig, ax = pl.PlotHFCIceEdge()
    fileIO.SaveFigures([fig], 'HFCIceEdgeStandardD')
    fileIO.SaveFigures([fig], 'HFCIceEdgeStandardD', '.svg')
    fig.show()
    
    pass


if __name__ == '__main__':
    pl.SetRCParams()
    main('smooth_coalbedo' in sys.argv)
