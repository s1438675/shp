"""
Matthew Christey-Reid
Senior Honours Project
s1438675@ed.ac.uk
"""

# Imports
from astropy.io import fits
import math
import numpy as np
import matplotlib.pyplot as plt

# Global Variables
_kbandarray = np.zeros((1, 3), dtype=float)
_hbandarray = np.zeros((1, 3), dtype=float)
_jbandarray = np.zeros((1, 3), dtype=float)
_ybandarray = np.zeros((1, 3), dtype=float)
_zbandarray = np.zeros((1, 3), dtype=float)
_epoch = 0.0


# Fold light curve
def foldcurve(_inputarray, _period):
    # Get global epoch variable
    global _epoch
    # Update x axis to phase
    for i in range(0, _inputarray.shape[0]):
        _inputarray[i, 0] = ((_inputarray[i, 0] - _epoch) / _period) % 1
    # Return updated array
    return _inputarray


# Apparent magnitude to absolute magnitude
def appmagnitudetoabsmagnitude(_inputarray, _pmas):
    # Distance in parsecs is 1 / parallax half angle
    _distance = 1 / (_pmas / 2)
    # Iterate through array
    for i in range(0, _inputarray.shape[0]):
        # Update magnitude from apparent magnitude to absolute magnitude
        _inputarray[i, 0] =_inputarray[i, 0] - (5 * (math.log10(_distance))) + 5
    # Return updated array
    return _inputarray

# Plot array using pyplot
def plotarray(_inputarray):
    # Create a new array twice the size of the input to plot from phase 0 -> 2
    _plotarray = np.zeros((_inputarray.shape[0] * 2, _inputarray.shape[1]), dtype=float)
    # Iterate through array
    for i in range(0, _plotarray.shape[0]):
        # Before phase 1 simply copy the data
        if i < _inputarray.shape[0]:
            _plotarray[i] = _inputarray[i]
        # After phase 1 shift all by +1
        else:
            _plotarray[i] = _inputarray[i - _inputarray.shape[0]]
            _plotarray[i][0] = _plotarray[i][0] + 1
    # Set pyplot type
    plt.style.use('seaborn-whitegrid')
    plt.errorbar(_plotarray[:, 0], _plotarray[:, 1], yerr=_plotarray[:, 2], fmt='.k')
    # Set axis values to the limits of the data array, considering error bar size
    axes = plt.gca()
    axes.set_xlim([0, 2])
    axes.set_ylim([np.amin(_plotarray[:, 1], axis=0) - np.amax(_plotarray[:, 2], axis=0),
                   np.amax(_plotarray[:, 1], axis=0) + np.amax(_plotarray[:, 2], axis=0)])
    # Determine line of best fit for data

    # Label x axis
    plt.xlabel("Phase")
    # Label y axis
    plt.ylabel("Magnitude")
    # Flip y axis as convention
    plt.gca().invert_yaxis()
    # Show plot
    plt.show()


# Plot all bands on the same graph
def plotbands(_inputarrayk, _inputarrayh, _inputarrayj, _inputarrayy, _inputarrayz):
    # Create a new array twice the size of the input to plot from phase 0 -> 2
    _plotarrayk = np.zeros((_inputarrayk.shape[0] * 2, _inputarrayk.shape[1]), dtype=float)
    _plotarrayh = np.zeros((_inputarrayh.shape[0] * 2, _inputarrayh.shape[1]), dtype=float)
    _plotarrayj = np.zeros((_inputarrayj.shape[0] * 2, _inputarrayj.shape[1]), dtype=float)
    _plotarrayy = np.zeros((_inputarrayy.shape[0] * 2, _inputarrayy.shape[1]), dtype=float)
    _plotarrayz = np.zeros((_inputarrayz.shape[0] * 2, _inputarrayz.shape[1]), dtype=float)
    # Iterate through array
    for i in range(0, _plotarrayh.shape[0]):
        # Before phase 1 simply copy the data
        if i < _inputarrayh.shape[0]:
            _plotarrayh[i] = _inputarrayh[i]
        # After phase 1 shift all by +1
        else:
            _plotarrayh[i] = _inputarrayh[i - _inputarrayh.shape[0]]
            _plotarrayh[i][0] = _plotarrayh[i][0] + 1
    # Set pyplot type
    plt.style.use('seaborn-whitegrid')
    plt.errorbar(_plotarrayk[:, 0], _plotarrayk[:, 1], yerr=_plotarrayk[:, 2], fmt='.k')
    plt.errorbar(_plotarrayh[:, 0], _plotarrayh[:, 1], yerr=_plotarrayh[:, 2], fmt='.k')
    plt.errorbar(_plotarrayj[:, 0], _plotarrayj[:, 1], yerr=_plotarrayj[:, 2], fmt='.k')
    plt.errorbar(_plotarrayy[:, 0], _plotarrayy[:, 1], yerr=_plotarrayy[:, 2], fmt='.k')
    plt.errorbar(_plotarrayz[:, 0], _plotarrayz[:, 1], yerr=_plotarrayz[:, 2], fmt='.k')
    # Set axis values to the limits of the data array, considering error bar size
    axes = plt.gca()
    axes.set_xlim([0, 2])
    # Determine line of best fit for data

    # Label x axis
    plt.xlabel("Phase")
    # Label y axis
    plt.ylabel("Magnitude")
    # Flip y axis as convention
    plt.gca().invert_yaxis()
    # Show plot
    plt.show()

# Main function
def main(_filename, _period):
    global _kbandarray
    global _hbandarray
    global _jbandarray
    global _ybandarray
    global _zbandarray
    # Open file with astropy fits library
    with fits.open(_filename) as _file:
        # Read in fits tables
        _data = _file[1].data
        # Set epoch to first data recorded
        _epoch = _data[0][0]
        for _row in _data:
            # If good k band reading, append to array
            if _row[1] > -100:
                _kbandarray = np.append(_kbandarray, np.array([_row[0], _row[1], _row[2]], ndmin=2), axis=0)
            # If a good h band reading, append to array
            if _row[3] > -100:
                _hbandarray = np.append(_hbandarray, np.array([_row[0], _row[3], _row[4]], ndmin=2), axis=0)
            # If a good j band reading, append to array
            if _row[5] > -100:
                _jbandarray = np.append(_jbandarray, np.array([_row[0], _row[5], _row[6]], ndmin=2), axis=0)
            # If a good y band reading, append to array
            if _row[7] > -100:
                _ybandarray = np.append(_ybandarray, np.array([_row[0], _row[7], _row[8]], ndmin=2), axis=0)
            # If a good z band reading, append to array
            if _row[9] > -100:
                _zbandarray = np.append(_zbandarray, np.array([_row[0], _row[9], _row[10]], ndmin=2), axis=0)
        # Delete blank first entry
        _kbandarray = np.delete(_kbandarray, 0, axis=0)
        _hbandarray = np.delete(_hbandarray, 0, axis=0)
        _jbandarray = np.delete(_jbandarray, 0, axis=0)
        _ybandarray = np.delete(_ybandarray, 0, axis=0)
        _zbandarray = np.delete(_zbandarray, 0, axis=0)
        # Fold curve
        _kbandarray = foldcurve(_kbandarray, _period)
        _hbandarray = foldcurve(_hbandarray, _period)
        _jbandarray = foldcurve(_jbandarray, _period)
        _ybandarray = foldcurve(_ybandarray, _period)
        _zbandarray = foldcurve(_zbandarray, _period)
        # Plot to screen
        plotbands(_kbandarray, _hbandarray, _jbandarray, _ybandarray, _zbandarray)


# Program entry point
main('Data/858994114024.fits', 0.558301)
