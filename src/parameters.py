### CLASSIC_EBM
### Jake Aylmer
###
### Parameters for use in the Python 2.7 package for the classic energy balance
### model (EBM).
### ---------------------------------------------------------------------------

from __future__ import division
import numpy as np

RE = 6371000.0 # Mean radius of Earth [m]

A = 203.3 # OLR expansion constant term [W m^-2]
B = 2.09 # OLR expansion linear term [W m^-2 degC^-1]

# Co-albedo expansion parameters [dimensionless]:
a0 = 0.697
a2 = -0.0779
b0 = 0.38

# Co-albedo typical value for ice-covered and ice-free regions respecively:
ai = 0.38 # [dimensionless]
af = 0.70 # [dimensionless]

delta_x = 0.3 # width of smoothing of co-albedo [dimensionless]

S2 = -0.482 # Coefficient of degree-2 Legendre polynomial in spatial
            # distribution of solar radiation [dimensionless]

Q = 335.0 # Solar constant divided by 4 [W m^-2]
            
D = 0.649 # Large scale diffusivity (accounts for geometric factors too)
          # [W m^-2 degC^-1]

T_ice_edge = -10.0 # temperature contour definining ice-cap edge [degC]

C = 0.16*B # effective heat capacity [W yr m^-2 degC^-1]


### ANALYTIC SOLUTION PARAMETERS ###
nmax = 6 # expansion index to truncate (see North et. al. 1981 equation (25))


### NUMERICAL SOLUTION PARAMETERS ###
T_equator_0 = 25 # initial equatorial temperature [degC]
xi_0 = 0.9 # initial sine of ice-edge latitude [dimensionless]

nbox = 100 # number of grid boxes to use (in x-space, with T defined at centres)
theta = 0.5 # parameter to determine which scheme to use for the diffusion eqn.
dt = 1/365 # time step [years]
tolerance = 0.01 # T-profile change for which steady-state may be assumed [degC]


### PLOTTING PARAMETERS ###
Q_min = 0.8 # default minimum extent to plot Q [units of default Q value]
Q_max = 1.4 # default maximum extent to plot Q [units of default Q value]
x_min = -0.05 # default minimum extent to plot x [dimensionless]
x_max = 1.05 # default maximum extent to plot x [dimensionless]
