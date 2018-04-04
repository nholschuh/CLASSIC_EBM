### CLASSIC EBM
### Jake Aylmer
###
### Functions for using the analytic solutions to the classic EBM.
### ---------------------------------------------------------------------------

from __future__ import division
import parameters as pm
import numpy as np
import scipy.special as spec, scipy.integrate as integrate


def Hn_step_coalbedo(n, xi):
    """The term H_n(x_i) appearing in the T_n coefficient of the analytic
    solution to the classical EBM (see North et. al. 1981 eq (29)). It uses
    the step-function coalbedo: 
    
            { ai    x >= xi
        a = {
            { af    x < xi
    
    --Args--
    n  : integer determining which term in the expansion is being calculated.
    xi : float between 0 and 1, sine of ice-edge latitude.
    """
    integrand = spec.legendre(n)*(1.0 + pm.S2*spec.legendre(2))
    integrand1 = integrand * pm.af
    integrand2 = integrand * pm.ai
    integral = np.polyint(integrand1)(xi) - np.polyint(integrand1)(0)
    integral += np.polyint(integrand2)(1) - np.polyint(integrand2)(xi)
    return (2*n + 1) * integral


def Hn_smooth_coalbedo(n, xi):
    """The term H_n(x_i) appearing in the T_n coefficient of the analytic
    solution to the classical EBM (see North et. al. 1981 eq (29)). It uses
    the smoothed coalbedo:
    
        a = a1 + a2*erf[ (x-xi)/delta_x ]
        
    where a1=(ai+af)/2, a2=(ai-af)/2, delta_x specifies the degree of 
    smoothing about the ice-edge location and erf() is the error-function.
    
    This function uses the SciPy general numerical integration method
    (scipy.integrate()) and SciPy Special module for the error function
    (scipy.special.erf()).
    
    --Args--
    n  : integer determining which term in the expansion is being calculated.
    xi : float between 0 and 1, sine of ice-edge latitude.
    """
    a1 = 0.5*(pm.ai+pm.af); a2 = 0.5*(pm.ai-pm.af)
    integrand = lambda x: ( (2*n+1)*spec.legendre(n)(x)*(
        1+pm.S2*spec.legendre(2)(x))*(a1+a2*spec.erf((x-xi)/pm.delta_x)) )
    return integrate.quad(integrand, 0.0, 1.0)[0]


def Ln(n, D=pm.D):
    """The term denoted L_n = n(n+1)D + B in the solution to the classical EBM
    (see North et. al. 1981 equation (28)).
    
    --Args--
    n   : int, identifies the term in the spectral expansion.
    (D) : float, large-scale constant diffusivity [W m^-2 degC^-1]
    """
    return n*(n+1)*D + pm.B


def Tn(n, xi, Q=pm.Q, D=pm.D, smooth_coalbedo=False):
    """
    Returns the coefficient T_n [degC] in the expansion of T(x) for the
    solution of the classical EBM model (North et. al. 1981 equation (30)).
    
    --Args--
    n                 : int, identifies the term in the spectral expansion.
    xi                : float, sine of ice-edge latitude [dimensionless].
    (Q)               : float, solar constant divided by 4 [W m^-2].
    (D)               : float, large-scale constant diffusivity
                        [W m^-2 degC^-1].
    (smooth_coalbedo) : bool, whether to use smoothed coalbedo function.
    """
    Hn= Hn_smooth_coalbedo(n,xi) if smooth_coalbedo else Hn_step_coalbedo(n,xi)
    T_n = Q*Hn/Ln(n, D) - (n==0)*(pm.A/pm.B)
    return T_n


def Q(xi, D=pm.D, smooth_coalbedo=False):
    """Calculates analytically Q at ice edge position xi for the diffusive
    model including the ice albedo feedback effect, i.e. equation (37) in North
    et al.
    
    --Args--
    xi                : float, sine of ice-edge latitude [dimensionless].
    (D)               : float, the large-scale constant diffusivity
                        [W m^-2 degC^-1]
    (smooth_coalbedo) : bool, whether to use the smoothed coalbedo function.
    """
    sumterm = 0
    for n in xrange(0, pm.nmax+2, 2):
        Hn = Hn_smooth(n, xi) if smooth_coalbedo else Hn_step_coalbedo(n, xi)
        sumterm += Hn*spec.legendre(n)(xi)/Ln(n, D)
    
    return (pm.A + pm.B*pm.T_ice_edge) / (pm.B*sumterm)


def HeatFluxConvergence(x, xi, Q=pm.Q, D=pm.D, smooth_coalbedo=False):
    """Calculate the steady-state heat flux convergence (HFC) [W m^-2] at
    location x.
    
    --Args--
    x                 : float, sine of latitude at which to calculate HFC.
    xi                : float, sine of ice-edge latitude.
    (Q)               : float, solar constant divided by 4 [W m^-2].
    (D)               : float, large-scale constant diffusivity
                        [W m^-2 degC^-1].
    (smooth_coalbedo) : bool, whether to use the smoothed coalbedo function.
    """
    sumterm_ddx = 0; sumterm_ddx2 = 0
    for n in xrange(0, pm.nmax+2, 2):
        T_n = Tn(n, xi, Q, D, smooth_coalbedo)
        sumterm_ddx += T_n * np.polyder(spec.legendre(n), m=1)(x)
        sumterm_ddx2 += T_n * np.polyder(spec.legendre(n), m=2)(x)
    return D * ( (1-x**2)*sumterm_ddx2 - 2*x*sumterm_ddx )


def HeatTransport(x, xi, Q=pm.Q, D=pm.D, smooth_coalbedo=False):
    """Calculate the steady-state zonally-integrated heat transport [W] at
    location x.
    
    --Args--
    x                 : float, sine of latitude at which to calculate HFC.
    xi                : float, sine of ice-edge latitude.
    (Q)               : float, solar constant divided by 4 [W m^-2].
    (D)               : float, large-scale constant diffusivity
                        [W m^-2 degC^-1].
    (smooth_coalbedo) : bool, whether to use the smoothed coalbedo function.
    """
    sumterm = 0 # sum from n=0 to n=nmax (n even), of T_n*(d/dx)P_n
    for n in xrange(0, pm.nmax+2, 2):
        T_n = Tn(n, xi, Q, D, smooth_coalbedo)
        sumterm += T_n*np.polyder(spec.legendre(n), m=1)(x)
    return -2*np.pi*D*pm.RE**2*(1-x**2)*sumterm
