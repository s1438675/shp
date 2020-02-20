"""
Matthew Christey-Reid
Senior Honours Project
s1438675@ed.ac.uk
"""

# Imports
from astropy.io import fits
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
    # Return new array
    return _inputarray


# Plot array using pyplot
def plotarray(_inputarray):
    _plotarray = np.zeros((_inputarray.shape[0] * 2, _inputarray.shape[1]), dtype=float)

    for i in range(0, _plotarray.shape[0]):
        if i < _inputarray.shape[0]:
            _plotarray[i] = _inputarray[i]
        else:
            _plotarray[i] = _inputarray[i - _inputarray.shape[0]]
            _plotarray[i][0] = _plotarray[i][0] + 1

    plt.style.use('seaborn-whitegrid')
    plt.errorbar(_plotarray[:, 0], _plotarray[:, 1], yerr=_plotarray[:, 2], fmt='.k')
    plt.xlabel("Phase")
    plt.ylabel("Magnitude")
    plt.gca().invert_yaxis()
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
        _kbandarray = foldcurve(_kbandarray, _period)
        plotarray(_kbandarray)


# Program entry point
main('Data/858994112734.fits', 0.236334)
