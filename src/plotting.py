### CLASSIC_EBM
### Jake Aylmer
###
### Plotting sub-routines.
### ---------------------------------------------------------------------------

from __future__ import division
import parameters as pm, analytics as an, math_methods as math
import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt


def StandardPlot(xi, norm_Q_arrays, relative_D, cols=['k']):
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
    ax.set_xlabel(r'Relative solar constant, $Q/Q_0$')
    ax.set_ylabel(r'Ice-edge position, $y_\mathrm{i}=\sin \phi_\mathrm{i}$')
    
    return FormatAxis(fig, ax, minorgrid=False)


def PlotHeatTransport(x, HT, xi):
    """Plot the zonally integrated heat transport in PW over the hemisphere for
    a given solution to the EBM (note that heat transports are input to this
    function in W, which are then converted to PW automatically).
    
    --Args--
    x  : (NumPy) array, containing x-coordinates between 0 and 1.
    HT : (NumPy) array, containing Heat transports [W] at each x coordinate.
    xi : float, sine of ice-edge latitude.
    """
    fig, ax = plt.subplots()
    ax.axvline(xi, linestyle='--', label=r'$y_\mathrm{i}$')
    ax.plot(x, HT/(1E15), color='k')
    ax.set_xlim([0,1])
    ax.set_xlabel(r'$y=\sin \phi$')
    ax.set_ylabel(r'Poleward Heat Transport (PW)')
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
