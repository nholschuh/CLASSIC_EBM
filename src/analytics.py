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
    Returns the coefficient T_n in the expansion of T(x) for the solution
    of the classical EBM model (North et. al. 1981 equation (30)).
    
    --Args--
    n                 : int, identifies the term in the spectral expansion.
    xi                : float, sine of ice-edge latitude [dimensionless].
    (Q)               : float, solar constant divided by 4 [W m^-2].
    (D)               : float, large-scale constant diffusivity
                        [W m^-2 degC^-1].
    (smooth_coalbedo) : bool, whether to use smoothed coalbedo function.
    """
    Hn = Hn_smooth(n, xi) if smooth_coalbedo else Hn_step_coalbedo(n, xi)
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
    (smooth_coalbedo) : bool, whether to use the smoothed albedo function.
    """
    sumterm = 0
    for n in xrange(0, pm.nmax+2, 2):
        Hn = Hn_smooth(n, xi) if smooth_coalbedo else Hn_step_coalbedo(n, xi)
        sumterm += Hn*spec.legendre(n)(xi)/Ln(n, D)
    
    return (pm.A + pm.B*pm.T_ice_edge) / (pm.B*sumterm)



def D(xs, A=203.3, B=2.09, Ts=-10.0, Q=335.0, ai=0.38, af=0.7, S2=-0.477):
    """
    Calculates D at ice edge position xs for the diffusive model including the
    ice albedo feedback. Experimental!
    """
    
    P2 = spec.legendre(2)(xs)
    P4 = spec.legendre(4)(xs)
    
    # Calculate J functions:
    integrand0 = 1.0 + S2*spec.legendre(2) # S(x)
    J = [0,0,0]
    
    for n in [0, 2, 4]:
        J[int(n/2)] = np.polyint( integrand0 * spec.legendre(n) * af)(xs
            ) - np.polyint( integrand0 * spec.legendre(n) * af)(0)
        J[int(n/2)] += np.polyint( integrand0 * spec.legendre(n) * ai
            )(1.0) - np.polyint( integrand0 * spec.legendre(n) * ai)(xs)

    J0 = J[0]
    J2 = J[1]
    J4 = J[2]
    
    alpha = 120.0*( ((A+B*Ts)/Q) - J0)
    beta = 2.0*B*( 13.0*(((A+B*Ts)/Q)-J0) - 50.0*J2*P2 - 27.0*J4*P4 )
    gamma = B*B*( ((A+B*Ts)/Q) - J0 - 5.0*J2*P2 - 9.0*J4*P4 )
    
    if beta**2 - 4.0*alpha*gamma < 0:
    
        return [np.nan, np.nan]
        
    else:
    
        Dplus = (-beta + np.sqrt(beta**2 - 4.0*alpha*gamma)) / (2.0*alpha)
        Dminus = (-beta - np.sqrt(beta**2 - 4.0*alpha*gamma)) / (2.0*alpha)
        
        if Dplus < 0:
            Dplus = np.nan
        if Dminus < 0:
            Dminus = np.nan
        
        return Dplus, Dminus



def heatfluxconvergence(x, xs, Q=335.0, D=0.649, A=203.3, B=2.09, af=0.7,
    ai=0.38, S2=-0.477, nmax=6):
    """
    Calculates the heat flux convergence at a given x for given parameters
    Q, D, A, B, xs, ai and af.
    """
    
    sumterm_ddx1 = 0
    sumterm_ddx2 = 0
    
    for n in xrange(0, nmax+2, 2):
    
        sumterm_ddx1 += tn_coef(n=n, Q=Q, D=D, A=A, B=B, amode='step', af=af,
            ai=ai, S2=S2, xs=xs)*np.polyder(spec.legendre(n), m=1)(x)
        sumterm_ddx2 += tn_coef(n=n, Q=Q, D=D, A=A, B=B, amode='step', af=af,
            ai=ai, S2=S2, xs=xs)*np.polyder(spec.legendre(n), m=2)(x)
    
    hfc = D*(1-x**2)*sumterm_ddx2 - 2*D*x*sumterm_ddx1
    
    return hfc



def heatflux(x, xs, Q=335.0, D=0.649, A=203.3, B=2.09, af=0.7, ai=0.38,
    S2=-0.477, R_E=6.371E6, nmax=6):
    """
    Calculates the heat flux at a given x for given parameters Q, D, A, B, xs,
    ai and af.
    """
    
    sumterm_ddx1 = 0
    for n in xrange(0, nmax+2, 2):
        sumterm_ddx1 += tn_coef(n=n, Q=Q, D=D, A=A, B=B, amode='step', af=af,
            ai=ai, S2=S2, xs=xs)*np.polyder(spec.legendre(n), m=1)(x)
    
    return -2* np.pi * R_E**2 * D * (1-x**2) * sumterm_ddx1



def heatfluxdensity(x, xs, Q=335.0, D=0.649, A=203.3, B=2.09, af=0.7, ai=0.38,
    S2=-0.477, nmax=6):
    """
    Calculates the heat flux density at a given x for given parameters Q, D, A,
    B, xs, ai and af - i.e. heatflux()/(area of latitude belt)
    """
    
    sumterm_ddx1 = 0
    for n in xrange(0, nmax+2, 2):
        sumterm_ddx1 += tn_coef(n=n, Q=Q, D=D, A=A, B=B, amode='step', af=af,
            ai=ai, S2=S2, xs=xs)*np.polyder(spec.legendre(n), m=1)(x)
    
    return - D * (1-x**2) * sumterm_ddx1
