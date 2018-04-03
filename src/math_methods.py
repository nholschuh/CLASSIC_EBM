### CLASSIC_EBM
### Jake Aylmer
###
### Mathematical methods.
### ---------------------------------------------------------------------------

from __future__ import division
import numpy as np


def SplitByGradient(x, y):
    """Split an array into multiple separate arrays based on where the gradient
    changes sign. For example, for the array [1,3,5,4,2,3,4], this function
    returns a NumPy array of arrays [ [1,3,5], [4,2], [3,4] ].
    
    If arrays containing consecutive equal values, e.g. [4,4,4,3,1], returns
    [ [4,4,4], [3,1] ].
    
    --Args--
    x, y : two (NumPy) arrays of equal length. Both are split according to the
           change in sign of the y array.
    """
    
    asc = (y[1] >= y[0]) #flag: True => in ascending/increasing part of array
    ind = [0] #indices of turning points.
    
    for j in xrange(2, len(y)):
        if (asc and y[j]<y[j-1]) or ( (not asc) and y[j]>y[j-1]):
            ind.append(j)
            asc = not asc #reverse flag
    
    if ind[-1] != len(y):
        ind.append(len(y))
    
    # Split original arrays:
    x_out = np.array( [x[ind[i-1]:ind[i]] for i in xrange(1,len(ind))] )
    y_out = np.array( [y[ind[i-1]:ind[i]] for i in xrange(1,len(ind))] )
    
    return x_out, y_out


def LinInt(xfind, xknown1, yknown1, xknown2, yknown2):
    """Estimate the value of y(xfind) given two known points (xknown1, yknown1)
    and (xknown2, yknown2) assuming a linear function for y."""
    return yknown1 + (yknown2 - yknown1)*(xfind - xknown1)/(xknown2 - xknown1)
