import math
import string
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.optimize import curve_fit

import lightcurve

from astropy import modeling
from astropy.timeseries import LombScargle

VegaToAB = [0.528, 0.634, 0.938, 1.379, 1.9]

def plotallbands(_zband, _yband, _jband, _hband, _kband, _period):
    plt.style.use('seaborn-whitegrid')
    # Frequency = 1 / Period
    _freq = 1 / _period
    i = 0
    while i < 5:
        _colours = ['-b', '-g', '-r', '-c', '-m']
        _bands = [_zband, _yband, _jband, _hband, _kband]
        _legend = ['Z-band', 'Y-band', 'J-band', 'H-band', 'K-band']
        _xfit, _lobf = calclobf(_bands[i], _period)
        # Plot the data in the array to screen, lightly coloured and z rank behind the line of best fit
        plt.plot(_xfit, _lobf, _colours[i], lw=1, zorder=2, label=_legend[i])
        i += 1
    plt.xlim(0, 1)
    plt.xlabel("Phase")
    plt.ylabel("Magnitude")
    plt.title("Folded light curve")
    plt.legend()
    # Invert y-axis as convention
    plt.gca().invert_yaxis()
    # Display to screen
    plt.show()


def plotband(_band, _period):
    """
    Plots an observed band using pyplot
    :param _band: Array to be plotted
    :param _period: Period of object
    """
    # Frequency = 1 / Period
    _freq = 1 / _period
    _xfit, _lobf = calclobf(_band, _period)
    # Plot the data in the array to screen, lightly coloured and z rank behind the line of best fit
    plt.errorbar((_band[:, 0] * _freq) % 1, _band[:, 1], _band[:, 2], fmt='.', color='gray',
                 ecolor='lightgray', capsize=0, zorder=0)
    plt.plot(_xfit, _lobf, '-k', lw=2, zorder=2)
    plt.xlim(0, 1)
    plt.xlabel("Phase")
    plt.ylabel("Magnitude")
    plt.title("Folded light curve")
    # Invert y-axis as convention
    plt.gca().invert_yaxis()
    # Display to screen
    plt.show()


def calclobf(_band, _period):
    """
    Creates a line of best fit using Lomb-Scargle methods
    :param _inputband: Band array to be fit
    :param _period: Period of object
    :return: Returns a linearly spaced x-axis, with y-axis values for line of best fit
    """
    # Create a model with 10 terms
    _ls = LombScargle(_band[:, 0], _band[:, 1], _band[:, 2], nterms=10)
    # Create n linearly spaced points between phase 0 and 1
    _xfit = np.linspace(0, 1, 1000)
    # Frequency = 1 / Period
    _freq = 1 / _period
    # Plot the line of best fit generated
    _lobf = _ls.model(_xfit / _freq, _freq)
    return _xfit, _lobf


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


def plotblackbody(_zband, _yband, _jband, _hband, _kband):
    _bands = [_zband, _yband, _jband, _hband, _kband]
    _lambda = [0.9, 1.02, 1.22, 1.63, 2.2]
    _largestar = np.zeros((1, 2))
    _smallstar = np.zeros((1, 2))

    i = 0
    while i < 5:
        _max, _min = lightcurve.maxminvals(_bands[i])
        _largestar = np.append(_largestar, np.array([_lambda[i], (magtoflux(_min, i + 1)) * 1000], ndmin=2), axis=0)
        i += 1
    _largestar = np.delete(_largestar, 0, axis=0)

    i = 0
    while i < 5:
        _max, _min = lightcurve.maxminvals(_bands[i])
        _smallstar = np.append(_smallstar, np.array([_lambda[i], (magtoflux(_max, i + 1) -
                                                                  magtoflux(_min, i + 1)) * 1000], ndmin=2), axis=0)
        i += 1
    _smallstar = np.delete(_smallstar, 0, axis=0)

    getwientemp(_largestar)
    getwientemp(_smallstar)



def getwientemp(_inputdata):
    _xdata = _inputdata[:, 0]
    _ydata = _inputdata[:, 1]
    #_params, _paramscovariance = curve_fit(maxwelleqn, _xdata, _ydata, method='dogbox')
    _popt, _pcov = curve_fit(stats.maxwell.pdf, _xdata, _ydata, p0=[0, 1], method='dogbox')
    _x = np.linspace(0, 3, 100)
    plt.style.use('seaborn-whitegrid')
    plt.plot(_inputdata[:, 0], _inputdata[:, 1], 'k.')
    _yplot = stats.maxwell.pdf(_x, *_popt)
    print(_pcov)
    plt.plot(_x, _yplot)

    _maxval, _maxloc = 0, 0
    i = 0
    while i < len(_x):
        if _yplot[i] > _maxval:
            _maxval = _yplot[i]
            _maxloc = i
        i += 1

    plt.plot()
    plt.show()
    _wientemp = 2898 / (_maxloc * (3 / 100))
    print(_wientemp)

def magtoflux(_mag, _id):
    return math.pow(10, -0.4*(_mag + VegaToAB[_id - 1] - 8.9))