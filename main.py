"""
Matthew Christey-Reid
Senior Honours Project
s1438675@ed.ac.uk
"""

# Imports
from astropy.io import fits
import numpy as np
import plotarray
import math


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
    _jband = createband(_data, 3)
    _yband = createband(_data, 4)
    _zband = createband(_data, 5)
    # Return created band arrays
    return _kband, _hband, _jband, _yband, _zband


def calcfluxdensity(_mag, _id):
    """
    Calculates the flux density of an object given a magnitude in the vega photometric system
    fb = fa * 10 ^ (-fb / 2.5)
    :param _mag: Magnitude of the star being calculated
    :param _id: Band id of observation
    :return: Returns the flux density in Jy
    """
    # Set 0 point for flux density using provided band id
    _zpoint = 0
    if _id == 1:
        _zpoint = 667
    if _id == 2:
        _zpoint = 1075
    if _id == 3:
        _zpoint = 1603
    # Return calculated flux density
    return _zpoint * math.pow(10, -_mag / 2.5)


def bandmaxmag(_inputband, _period):
    """
    Calculate the maximum magnitude of the line of best fit for an observation band
    :param _inputband: Observation band to be worked on
    :param _period: Period of object
    :return: Maximum magnitude of the object
    """
    # Plot a line of best fit using Lomb-Scargle formulae
    _xfit, _lobf = plotarray.calclobf(_inputband, _period)
    # Get critical points of lobf
    _max, _min = plotarray.getcritpoints(_lobf)
    # Determine maximum magnitude
    if _lobf[_min[0]] > _lobf[_min[1]]:
        _magmax = _lobf[_min[0]]
    else:
        _magmax = _lobf[_min[1]]
    # Return found value
    return _magmax


def calcpeakfluxdensity(_inputband, _period, _id):
    """
    Calculate the peak flux density of a band
    :param _inputband: Band to be worked on
    :param _period: Period of object
    :param _id: Band id
    :return: Returns the flux density of the peak magnitude observed in Jy
    """
    # Determine maximum magnitude of the object
    _magmax = bandmaxmag(_inputband, _period)
    # Calculate flux density from max magnitude
    _fluxdensity = calcfluxdensity(_magmax, _id)
    # Return calculated flux density
    return _fluxdensity


def calcflux(_fluxdensity, _bandwidth):
    """
    Calculate the flux from flux density
    F = f * dlambda
    :param _fluxdensity: Flux density of object per unit wavelength in W/m^2/Hz
    :param _bandwidth: Width of the observation band in metres
    :return: Returns the flux of an object
    """
    return _fluxdensity * _bandwidth


def calcdistance(_parallax):
    """
    Calculate the distance to an object using parallax angle
    :param _parallax: Parallax angle to be converted to distnace
    :return: Distance to object in metres
    """
    return (1 / (math.tan(_parallax / 1000))) * 3.086 * (math.pow(10, 16))


def main(_filename, _period, _parallax):
    """
    Main method of the program, reads in a FITS table and determines stellar properties from spectroscopic data
    within the file
    :param _filename: FITS file to be opened
    :param _period: Period of object
    :param _parallax: Parallax angle of object
    """
    # Open the file with astropy.io.fits library
    with fits.open(_filename) as _file:
        # Read data in file, file[0].data contains header information
        _data = _file[1].data
        # Split FITS file into seperate bands
        _kband, _hband, _jband, _yband, _zband = splitbands(_data)
        _kfluxd = calcpeakfluxdensity(_kband, _period, 1)
        _hfluxd = calcpeakfluxdensity(_hband, _period, 2)
        _jfluxd = calcpeakfluxdensity(_jband, _period, 3)
        _bb = np.array([[1.2, _jfluxd], [1.6, _hfluxd], [2.2, _kfluxd]], dtype=float)
        plotarray.plotarray(_bb)
        plotarray.plotbands(_kband, _hband, _jband, _period)
        plotarray.plotlobf(_jband, _period)
        plotarray.plotblackbody(_bb)


# Program entry point
main('Data/858994114024.fits', 0.558301, 0.5187)
