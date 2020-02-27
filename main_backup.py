"""
Matthew Christey-Reid
Senior Honours Project
s1438675@ed.ac.uk
"""

# Imports
from astropy.io import fits
import numpy as np
import plotarray


# Fold light curve+
def foldcurve(_inputarray, _period):
    _epoch = _inputarray[0][0]
    # Update x axis to phase
    for i in range(0, _inputarray.shape[0]):
        _inputarray[i, 0] = ((_inputarray[i, 0] - _epoch) / _period) % 1
    # Return updated array
    return _inputarray


# Split fits array into k, h, j, y, z bands
def splitbands(_inputarray):
    _kband = np.zeros((1, 3), dtype=float)
    _hband = np.empty((1, 3), dtype=float)
    _jband = np.empty((1, 3), dtype=float)
    _yband = np.empty((1, 3), dtype=float)
    _zband = np.empty((1, 3), dtype=float)
    for _row in _inputarray:
        # If good k band reading, append to array
        if _row[1] > -100:
            _kband = np.append(_kband, np.array([_row[0], _row[1], _row[2]], ndmin=2), axis=0)
        # If a good h band reading, append to array
        if _row[3] > -100:
            _hband = np.append(_hband, np.array([_row[0], _row[3], _row[4]], ndmin=2), axis=0)
        # If a good j band reading, append to array
        if _row[5] > -100:
            _jband = np.append(_jband, np.array([_row[0], _row[5], _row[6]], ndmin=2), axis=0)
        # If a good y band reading, append to array
        if _row[7] > -100:
            _yband = np.append(_yband, np.array([_row[0], _row[7], _row[8]], ndmin=2), axis=0)
        # If a good z band reading, append to array
        if _row[9] > -100:
            _zband = np.append(_zband, np.array([_row[0], _row[9], _row[10]], ndmin=2), axis=0)
    _kband = delfirst(_kband)
    return _kband, _hband, _jband, _yband, _zband


def delfirst(_inputarray):
    _inputarray = np.delete(_inputarray, 0, axis=0)
    return _inputarray


# Main function
def main(_filename, _period):
    # Open file with astropy fits library
    with fits.open(_filename) as _file:
        # Read in fits tables
        _data = _file[1].data
        _kband, _hband, _jband, _yband, _zband = splitbands(_data)
        # _kband = foldcurve(_kband, _period)
        plotarray.plotlobf(_kband, _period)


# Program entry point
main('Data/858994114024.fits', 0.558301)
