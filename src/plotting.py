### CLASSIC_EBM
### Jake Aylmer
###
### Plotting sub-routines.
### ---------------------------------------------------------------------------

from __future__ import division
import parameters as pm, analytics as an, math_methods as math
import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt


def StabilityPlot(xi, norm_Q_arrays, relative_D, cols=['k']):
    """Generate a plot of x_i against Q given multiple data sets for different
    values of D. Plots are added with labels so that a legend may be added to 
    the returned MatPlotLib axis object if desired. (Returns (fig, ax)).
    
    --Args--
    xi            : NumPy array, coordinates of sine of ice-edge latitude.
    norm_Q_arrays : NumPy array of shape (len(relative_D), len(xi)). Contains
                    a set of Q(x_i) data for each value of D = relative_D*D0.
    relative_D    : NumPy array containing the values of D relative to (i.e. in
                    units of) the standard value of D (D_0), in the order
                    corresponding to the data sets in norm_Q_arrays.
    (cols)        : array containing MatPlotLib color identifiers which are
                    selected sequentially as the plots are added, cycling back
                    to the beginning if len(cols) < len(relative_D).
    """
    fig, ax = plt.subplots()
    ax.axhline(0.0, color=[.5,.5,.5], linewidth=0.8)
    ax.axhline(1.0, color=[.5,.5,.5], linewidth=0.8)
    
    for k in xrange(len(norm_Q_arrays)):
        
        xi_split, Q_split = math.SplitByGradient(xi, norm_Q_arrays[k])
        
        for j in xrange(len(Q_split)):
           lnst = '--' if Q_split[j][1]<Q_split[j][0] else '-'
           ax.plot(Q_split[j], xi_split[j], color=cols[k%len(cols)],
               linestyle=lnst)
    
        # Add perennial ice free and snowball states:
        xi_cold = np.array([0, 0])
        Q_cold = np.array([pm.Q_min, norm_Q_arrays[k][0]])
        xi_warm = np.array([1, 1])
        Q_warm = np.array([norm_Q_arrays[k][-1], pm.Q_max])
        ax.plot(Q_cold, xi_cold, color=cols[k%len(cols)], linestyle='-')
        ax.plot(Q_warm, xi_warm, color=cols[k%len(cols)], linestyle='-',
            label=r'$D/D_0=%.2f$' % relative_D[k])
    
    ax.set_xlim([pm.Q_min, pm.Q_max])
    ax.set_ylim([pm.x_min, pm.x_max])
    ax.set_xlabel(r'Normalised solar constant, $Q/Q_0$')
    ax.set_ylabel(r'Ice-edge position, $x_\mathrm{i}=\sin \phi_\mathrm{i}$')
    fig.canvas.set_window_title(
        'StabilityPlot' + ('Multiple'*(len(norm_Q_arrays)>1)) )
    return FormatAxis(fig, ax, minorgrid=False)


def PlotTemperature(x, T, xi):
    """"""
    fig, ax1 = plt.subplots()
    ax1.axvline(xi, linestyle='--', color='b')
    ax1.plot(x, T, color='k')
    
    ax2 = ax1.twiny()
    ax2.set_xticks( np.sin((np.pi/180)*np.arange(0, 90.1, 10)) )
    ax2.set_xticklabels( np.arange(0, 90, 10) )
    
    ax1.set_xlabel(r'$x=\sin \phi$')
    ax1.set_ylabel(r'Surface temperature, $T$ ($^\circ$C)')
    ax2.set_xlabel(r'Latitude, $\phi$ (deg)', y=2)
    ax2.tick_params(axis='both', which='both', direction='out')
    ax2.tick_params(axis='both', which='major', labelsize=17, pad=0)
    return FormatAxis(fig, ax1, minorgrid=False)


def PlotHeatTransport(x, HT, xi, latitude_axis=False):
    """Plot the zonally integrated heat transport in PW over the hemisphere for
    a given solution to the EBM (note that heat transports are input to this
    function in W, which are then converted to PW automatically). Returns the
    MatPlotLib figure and axis objects (fig, ax).
    
    --Args--
    x               : (NumPy) array, containing x-coordinates between 0 and 1.
    HT              : (NumPy) array, containing Heat transport [W] at each x
                      coordinate.
    xi              : float, sine of ice-edge latitude.
    (latitude_axis) : bool, whether to convert to latitude (deg).
    """
    fig, ax = plt.subplots()
    if latitude_axis:
        x = np.degrees(np.arcsin(x))
        xi = np.degrees(np.arcsin(xi))
        ax.set_xlim([0,90])
        ax.set_xlabel(r'Latitude, $\phi$ ($^\circ$)')
    else:
        ax.set_xlim([0,1])
        ax.set_xlabel(r'$x=\sin \phi$')
    ax.axvline(xi, linestyle='--', label=r'Ice edge')
    ax.plot(x, HT/(1E15), color='k')
    ax.set_ylabel(r'Poleward Heat Transport (PW)')
    fig.canvas.set_window_title('HeatTransport')
    return FormatAxis(fig, ax, minorgrid=False)


def PlotHeatFluxConvergence(x, HFC, xi, latitude_axis=False):
    """Plot the heat flux convergence (HFC) [W m^-2] over the hemisphere for a
    given solution to the EBM. Returns the MatPlotLib figure and axis objects
    (fig, ax).
    
    --Args--
    x               : (NumPy) array, containing x-coordinates between 0 and 1.
    HFC             : (NumPy) array, containing heat flux convergences [W m^-2]
                      at each x.
    xi              : float, sine of ice-edge latitude.
    (latitude_axis) : bool, whether to convert to latitude (deg).
    """
    fig, ax = plt.subplots()
    ax.axhline(0, color=[.2,.2,.2], linewidth=0.8)
    if latitude_axis:
        x = np.degrees(np.arcsin(x))
        xi = np.degrees(np.arcsin(xi))
        ax.set_xlim([0,90])
        ax.set_xlabel(r'Latitude, $\phi$ ($^\circ$)')
    else:
        ax.set_xlim([0,1])
        ax.set_xlabel(r'$x=\sin \phi$')
    ax.axvline(xi, linestyle='--', label=r'Ice edge')
    ax.plot(x, HFC, color='k')
    ax.set_ylabel(r'Heat flux convergence (W m$^{-2}$)')
    fig.canvas.set_window_title('HeatFluxConvergence')
    return FormatAxis(fig, ax, minorgrid=False)


def PlotHFCIceEdge(relative_D=np.array([0.75,1.0,1.25]), smooth_coalbedo=False,
    add_linear_fit=False):
    """Plot the heat flux convergence (HFC) at the ice edge as the ice edge
    varies (i.e. HFC(x=xi) vs xi) for each value of D = relative_D * D0 where
    D0 is standard value (set in the parameters file). Calculations are done
    here and the MatPlotLib figure and axis objects (fig, ax) are returned.
    
    --Args--
    (relative_D)      : (NumPy) array of values of D to be used in units of D0.
    (smooth_coalbedo) : bool, whether to use the smoothed coalbedo function.
    (add_linear_fit)  : bool, if True, adds a linear fit to the first data set.
    """
    
    xi = np.arange(0.0, 1.001, 0.01)
    HFC = np.zeros( (len(relative_D), len(xi)) )
    
    for j in xrange(len(relative_D)):
        for k in xrange(len(xi)):
            Q = an.Q(xi[k], pm.D*relative_D[j], smooth_coalbedo)
            HFC[j][k] = an.HeatFluxConvergence(xi[k], xi[k], Q,
                pm.D*relative_D[j], smooth_coalbedo)
    
    fig, ax = plt.subplots()
    ax.axhline(0, color=[.2,.2,.2], linewidth=0.8)
    firstplot = ax.plot(xi, HFC[0], label=r'$D/D_0=%.2f$' % relative_D[0])
    
    if add_linear_fit:
        a = np.argmin(abs(xi-pm.xi_HFC_lim1)) # index of lowest xi to fit to
        b = np.argmin(abs(xi-pm.xi_HFC_lim2)) # index of upper xi to fit to
        fit = np.polyfit(xi[a:b], HFC[0][a:b], 1)
        ax.plot(xi[a:b], fit[1]+fit[0]*xi[a:b], linestyle='--',
            color=firstplot[0].get_color() )
    
    for j in xrange(1, len(relative_D)):
        ax.plot(xi, HFC[j], label=r'$D/D_0=%.2f$' % relative_D[j])
    ax.set_xlim([0,1])
    ax.set_xlabel(r'Ice edge position, $x_\mathrm{i}=\sin\phi_\mathrm{i}$')
    ax.set_ylabel(r'Heat flux convergence (W m$^{-2}$)')
    ax.legend(loc='upper left', fontsize=16)
    fig.canvas.set_window_title(
        'HeatFluxConvergenceIceEdge' + '_Multiple'*(len(relative_D)>1))
    return FormatAxis(fig, ax, minorgrid=False)


###############################################################################

def FormatAxis(fig, ax, ticksize=18, tickpad=8, gridon=True, minorgrid=True):
    """Set the layout and formatting of the plot on axis ax belonging to figure
    object fig.
    
    --Args--
    fig         : MatPlotLib figure object.
    ax          : MatPlotLib axis object associated with fig.
    (ticksize)  : int, font-size for axis tick labels.
    (tickpad)   : int, padding for the axis tick labels (see MatPlotLib).
    (gridon)    : boolean, whether to set the grid on.
    (minorgrid) : boolean, whether to show minor grid-lines.
    """
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='out')
    ax.tick_params(axis='both', which='major', labelsize=ticksize, pad=tickpad)
    if gridon:
        if minorgrid:
            ax.grid(which='minor', linestyle='-', color=[.92, .92, .92])
        ax.grid(which='major', linestyle='-', color=[.75, .75, .75])
    ax.set_axisbelow(True)
    fig.tight_layout()
    return fig, ax


def SetRCParams():
    """Set default MatPlotLib formatting styles (rcParams) which will be set
    automatically for any plotting method.
    """
    mpl.rcParams['font.sans-serif'] = 'Calibri' #set the font for sans-serif style
    mpl.rcParams['font.family'] = 'sans-serif' #choose the sans-serif font style
    mpl.rcParams['mathtext.fontset'] = 'custom' #allow customisation of maths font
    mpl.rcParams['mathtext.rm'] = 'sans' #maths roman font in sans-serif format
    mpl.rcParams['mathtext.it'] = 'sans:italic' #maths italic font
    mpl.rcParams['mathtext.default'] = 'it' #maths in italic by default
    mpl.rcParams['axes.titlesize'] = 20 #plt.title font size
    mpl.rcParams['axes.labelsize'] = 18 #x/y axis label font size
    mpl.rcParams['savefig.format'] = 'pdf' #default format to save to
    mpl.rcParams['lines.linewidth'] = 1.5 #default plot linewidth (thickness)
    pass
