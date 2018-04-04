### CLASSIC_EBM
### Jake Aylmer
###
### Generate a single plot containing two Q(xi) profiles, one using the 
### original value of D = D0, and the other using a different value D = f*D0.
### ---------------------------------------------------------------------------

from __future__ import division

import sys, os, numpy as np, matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from src import parameters as pm, analytics as an, fileIO
from src import plotting as pl


def main(f=0.7, smooth_coalbedo=False):
    
    xi = np.arange(0.0, 1.0001, 0.001)
    relative_D = np.array([1.0, f])
    Q_arrays = np.zeros( (len(relative_D), len(xi)) )
    
    print "Using %s-coalbedo..." % ('smoothed' if smooth_coalbedo else 
        'step-function')
    
    for k in xrange(len(Q_arrays)):
        print "Calculating with D/D0 = %.2f" % relative_D[k]
        for j in xrange(len(xi)):
            Q_arrays[k][j] = an.Q(xi[j], D=relative_D[k]*pm.D, 
                smooth_coalbedo=smooth_coalbedo)
    
    fig, ax = pl.StabilityPlot(xi, Q_arrays/pm.Q, relative_D, ['grey', 'k'])
    ax.legend(loc='upper right')
    subdir_name =('change_D_to_%.2fD0'%f)+('_SmoothedCoalbedo'*smooth_coalbedo)
    fileIO.SaveFigures([fig], subdir_name)
    fileIO.SaveFigures([fig], subdir_name, '.svg')
    fig.show()
    
    pass


if __name__ == '__main__':
    pl.SetRCParams()
    if len(sys.argv)==3:
        try:
            main(float(sys.argv[1]), (sys.argv[2]=='smooth_coalbedo'))
        except:
            main()
    elif len(sys.argv)==2 and sys.argv[1]=='smooth_coalbedo':
        main(smooth_coalbedo=True)
    elif len(sys.argv)==2:
        main(float(sys.argv))
    else:
        main()
