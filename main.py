"""
Matthew Christey-Reid
Senior Honours Project
s1438675@ed.ac.uk
"""

# Imports
from astropy.io import fits
from scipy.fftpack import fft, ifft
import math
import numpy as np
import matplotlib.pyplot as plt

# Global Variables
_kbandarray = np.zeros((1, 3), dtype=float)
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
    # Set pyplot stype
    plt.style.use('seaborn-whitegrid')
    plt.errorbar(_plotarray[:, 0], _plotarray[:, 1], yerr=_plotarray[:, 2], fmt='.k')
    # Set axis values to the limits of the data array, considering error bar size
    axes = plt.gca()
    axes.set_xlim([0, 2])
    axes.set_ylim([np.amin(_plotarray[:, 1], axis=0) - np.amax(_plotarray[:, 2], axis=0),
                   np.amax(_plotarray[:, 1], axis=0) + np.amax(_plotarray[:, 2], axis=0)])
    # Determine line of best fit for data
    _X = fft(_plotarray[:, 0])
    _x = ifft(_plotarray[:, 0])
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
    # Open file with astropy fits library
    with fits.open(_filename) as _file:
        # Read in fits tables
        _data = _file[1].data
        # Set epoch to first data recorded
        _epoch = _data[0][0]
        for _row in _data:
            # If a good reading, add to array
            if _row[1] > -100:
                _kbandarray = np.append(_kbandarray, np.array([_row[0], _row[1], _row[2]], ndmin=2), axis=0)
        # Delete blank first entry
        _kbandarray = np.delete(_kbandarray, 0, axis=0)
        # Fold curve
        _kbandarray = foldcurve(_kbandarray, _period)
        # Plot to screen
        plotarray(_kbandarray)


# Program entry point
main('Data/858994114024.fits', 0.558301)
