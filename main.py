"""
Matthew Christey-Reid
Senior Honours Project
s1438675@ed.ac.uk
"""

# Imports
from astropy.io import fits
import numpy as np
import plotarray


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


def splitbands(_data):
    """
    Split observed bands in FITS file to seperate array for each band
    :param _data: Data from FITS file to be split
    :return: Returns k, h, j, y, z arrays - where each array is [date][mag][mag_err]
    """
    def createband(_data, _id):
        """
        Creates a new array for measurements in a band
        :param _data: FITS file to be split
        :param _id: ID of band, k = 1, h = 2, j = 3, y = 4, z = 5
        :return: Returns band data in a [n, 3] array
        """
        # Create empty array to store values
        _band = np.zeros((1, 3), dtype=float)
        # Go through each row of FITS data
        for _row in _data:
            # If a good reading, append to array
            if _row[(_id * 2) - 1] > -100:
                _band = np.append(_band, np.array([_row[0], _row[(_id * 2) - 1], _row[_id * 2]], ndmin=2), axis=0)
        # Delete the first (empty) row in the array
        _band = delfirst(_band)
        # Return band array
        return _band

    def delfirst(_band):
        """
        Deletes the first entry in an array
        :param _band: Band array to have first entry deleted
        :return: Array of size [n - 1, 3], empty first row deleted
        """
        # Delete the first row of array
        _band = np.delete(_band, 0, axis=0)
        # Return new array
        return _band
    # Create the band arrays from data and band ids
    _kband = createband(_data, 1)
    _hband = createband(_data, 2)
    # Return created band arrays
    return _kband, _hband


def main(_filename, _period):
    """
    Main method of the program, reads in a FITS table and determines stellar properties from spectroscopic data
    within the file
    :param _filename: FITS file to be opened
    :param _period: Period of object
    """
    # Open the file with astropy.io.fits library
    with fits.open(_filename) as _file:
        # Read data in file, file[0].data contains header information
        _data = _file[1].data
        # Split FITS file into seperate bands
        _kband, _hband = splitbands(_data)
        # Plot split bands
        plotarray.plotlobf(_kband, _period)
        plotarray.plotlobf(_hband, _period)


# Program entry point
main('Data/858994114024.fits', 0.558301)
