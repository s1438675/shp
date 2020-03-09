"""
Matthew Christey-Reid
Senior Honours Project
s1438675@ed.ac.uk
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.optimize import curve_fit
from astropy.timeseries import LombScargle


def foldcurve(_band, _period):
    """
    Folds the magnitude measurements to a light curve using provided period
    :param _band: Observation band to be folded
    :param _period: Period of object
    :return: Array same size as _band, but with a phase instead of Julian date
    """
    # Set epoch to first date observed
    _epoch = _band[0][0]
    # Iterate through array, update date to phase
    for i in range(0, _band.shape[0]):
        _band[i, 0] = ((_band[i, 0] - _epoch) / _period) % 1
    # Return folded array
    return _band


def doublearrayphase(_inputarray):
    """
    Doubles a band array from phase 0 -> 1 to 0 -> 2 as convention
    :param _inputarray: Array to be doubled
    :return: Returns an array from phase 0 -> 2, size [n * 2, 3]
    """
    # Create a new array twice the size of the input
    _newarray = np.zeros((_inputarray.shape[0] * 2, _inputarray.shape[1]), dtype=float)
    # Iterate through the input array
    for i in range(0, _newarray.shape[0]):
        # Before phase 1 simply copy data into new array
        if i < _inputarray.shape[0]:
            _newarray[i] = _inputarray[i]
        # After phase 1, simply shift all phases by +1
        else:
            _newarray[i] = _inputarray[i - _inputarray.shape[0]]
            _newarray[i][0] = _newarray[i][0] + 1
    # Return the new doubled array
    return _newarray


def plotbands(_kband, _hband, _jband, _period):
    kx, kl = calclobf(_kband, _period)
    hx, hl = calclobf(_hband, _period)
    jx, jl = calclobf(_jband, _period)
    plt.plot(kx, kl, 'b.')
    plt.plot(hx, hl, 'r.')
    plt.plot(jx, jl, 'k.')
    plt.show()

def plotband(_band, _period):
    """
    Plots an observed band using pyplot
    :param _band: Array to be plotted
    :param _period: Period of object
    """
    _band = foldcurve(_band, _period)
    _band = doublearrayphase(_band)
    plt.plot(_band[:, 0], _band[:, 1], 'b.')
    plt.gca().invert_yaxis()
    plt.show()


def calclobf(_inputband, _period):
    """
    Creates a line of best fit using Lomb-Scargle methods
    :param _inputband: Band array to be fit
    :param _period: Period of object
    :return: Returns a linearly spaced x-axis, with y-axis values for line of best fit
    """
    # Create a model with 10 terms
    _ls = LombScargle(_inputband[:, 0], _inputband[:, 1], _inputband[:, 2], nterms=10)
    # Create n linearly spaced points between phase 0 and 1
    _xfit = np.linspace(0, 1, 1000)
    # Frequency = 1 / Period
    _freq = 1 / _period
    # Plot the line of best fit generated
    _lobf = _ls.model(_xfit / _freq, _freq)
    return _xfit, _lobf


def plotlobf(_inputband, _period):
    """
    Plots line of best fit to screen using pyplot
    :param _inputband: Band of array to be plotted
    :param _period: Period of object
    """
    # Frequency = 1 / Period
    _freq = 1 / _period
    _xfit, _lobf = calclobf(_inputband, _period)
    # Plot the data in the array to screen, lightly coloured and z rank behind the line of best fit
    plt.errorbar((_inputband[:, 0] * _freq) % 1, _inputband[:, 1], _inputband[:, 2], fmt='.', color='gray',
                 ecolor='lightgray', capsize=0, zorder=0)
    plt.plot(_xfit, _lobf, '-k', lw=2, zorder=2)
    # Get maxima and minima of the fit
    _max, _min = getcritpoints(_lobf)
    # Plot maxima and minima one the line of best fit
    for i in _min:
        plt.plot(_xfit[i], _lobf[i], '.r', zorder=3)
    for i in _max:
        plt.plot(_xfit[i], _lobf[i], '.b', zorder=3)
    # Invert y-axis as convention
    plt.gca().invert_yaxis()
    # Display to screen
    plt.show()


def getcritpoints(_inputarray):
    """
    Get the critical points of an array
    :param _inputarray: Array to get critical points [x, y]
    :return: Two arrays, value is the index of critical point in _inputarray
    """
    _max, _min = [], []
    if len(_inputarray) < 3:
        return _min, _max

    NEUTRAL, RISING, FALLING = range(3)

    def get_state(a, b):
        if a < b: return RISING
        if a > b: return FALLING
        return NEUTRAL

    ps = get_state(_inputarray[0], _inputarray[1])
    begin = 1
    for i in range(2, len(_inputarray)):
        s = get_state(_inputarray[i - 1], _inputarray[i])
        if s != NEUTRAL:
            if ps != NEUTRAL and ps != s:
                if s == FALLING:
                    _max.append((begin + i - 1) // 2)
                else:
                    _min.append((begin + i - 1) // 2)
            begin = i
            ps = s
    return _min, _max


def plotarray(_inputarray):
    plt.plot(_inputarray[:, 0], _inputarray[:, 1], 'b.')
    plt.xlim([0, 3])
    plt.ylim([0, 0.001])
    plt.show()


def plotblackbody(_inputarray):
    """
    mean, std = stats.norm.fit(_inputarray)
    x = np.linspace(0, 5, 100)
    y = stats.norm.pdf(x, mean, std)
    plt.plot(x, y)
    plt.show()
    return
    """
    params = stats.maxwell.fit(_inputarray, floc=0)
    x = np.linspace(0, 5, 100)
    plt.plot(x, stats.maxwell.pdf(x, *params))
    plt.show()
