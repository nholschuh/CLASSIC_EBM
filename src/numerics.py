### CLASSIC EBM
### Jake Aylmer
###
### Functions for solving the classic EBM numerically.
### ---------------------------------------------------------------------------

from __future__ import division
import parameters as pm, math_methods as math, diffusion_scheme as dif
import numpy as np, scipy.special as spec


def SmoothedCoalbedo(x, xi):
    """"""
    return 0.5*(pm.af+pm.ai) + (pm.ai-pm.af)*spec.erf((x-xi)/pm.delta_x)


def PseudoDiffusivity(x, D=pm.D):
    """"""
    return (1-x**2)*D/pm.C


def CalculateIceEdge(x, T):
    """Determine the ice edge location from the temperature profile."""
    sign = T > pm.T_ice_edge
    if all(sign):
        xi = 1.0 # all T > T_ice_edge
    elif all(np.invert(sign)):
        xi = 0.0 # all T < T_ice_edge
    else:
        for Tj in xrange(len(sign)-1):
            if sign[Tj+1] != sign[Tj]:
                xi = math.LinInt(pm.T_ice_edge, T[Tj], x[Tj], T[Tj+1], x[Tj+1])
                break
    return xi


def InitialConditions():
    """"""
    x = np.linspace(0.5/pm.nbox, 1-0.5/pm.nbox, pm.nbox)
    T_init = pm.T_equator_0*(1.0 - (x/pm.xi_0)**2) + (x/pm.xi_0)*pm.T_ice_edge   
    return x, T_init


def FindQofxi(Q=np.arange(0.9,1.21,0.1)):
    
    xi_of_Q = np.zeros(len(Q))
    
    for j in xrange(len(Q)):
        print "Calculating Q=%.2f*Q0 (%i of %i)" % (Q[j], j+1, len(Q))
        x, T_old = InitialConditions()
        xi = CalculateIceEdge(x, T_old)
        
        integrating = True
        while integrating:
            RHS = (pm.Q*Q[j]*(1+pm.S2*spec.legendre(2)(x))*SmoothedCoalbedo(x,xi) -
                pm.A - pm.B*T_old) / pm.C
            T_new = dif.SolveDiffusionEquation(T_old, RHS, RHS, PseudoDiffusivity,
                pm.dt, theta=pm.theta)
            xi = CalculateIceEdge(x, T_new)
            
            if all(abs(T_new-T_old)<=pm.tolerance):
                integrating = False
                xi_of_Q[j] = xi
            else:
                T_old = T_new.copy()
    
    return xi_of_Q
