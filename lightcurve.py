# Import numpy for data processing
import numpy as np


def splitbands(_data):
    """
    Split observed bands in FITS file to seperate array for each band
    :param _data: Data from FITS file to be split
    :return: Returns z, y, j, h, k arrays - where each array is [date][mag][mag_err]
    """
    def createband(_data, _id):
        """
        Creates a new array for measurements in a band
        :param _data: FITS file to be split
        :param _id: ID of band, z = 1, y = 2, j = 3, h = 4, k = 5
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
    _zband = createband(_data, 1)
    _yband = createband(_data, 2)
    _jband = createband(_data, 3)
    _hband = createband(_data, 4)
    _kband = createband(_data, 5)
    # Return created band arrays
    return _zband, _yband, _jband, _hband, _kband


def maxminvals(_band):
    """
    Determine the minimum and maximum values of a given array
    :param _band: The observed band to have values returned from
    :return: Two values, maximum and minimum of given array
    """
    _max = 100
    _min = 0
    for _coord in _band:
        if _coord[1] < _max:
            _max = _coord[1]
        if _coord[1] > _min:
            _min = _coord[1]
    return _max, _min
