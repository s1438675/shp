#!/usr/bin/env python
"""
Author: Matthew Christey-Reid
Email:  s1438675@ed.ac.uk
Date:   02/06/2020

main.py - Main entry point of the program
main(): Takes 3 arguments, FITS file, period (days) and parallax angle (mas)
"""

# Import library for passing arguments
import sys

# Import library for parsing FITS files
from astropy.io import fits

# Import libraries for determining properties
import plot
import lightcurve


def main():
    """
    Opens provided file and passes raw data, period and parallax to appropriate functions
    """
    # Attempt to parse arguments
    try:
        _args = sys.argv
        _filename = _args[1]
        _period = float(_args[2])
        _parallax = float(_args[3])
    # Throw error when incorrectly called
    except IndexError:
        print("Error loading files, please check your arguments and try again...")
        exit()

    # Open the file using the astropy FITS library
    with fits.open(_filename) as _file:
        # Read contents into array
        _data = _file[1].data
        # Split array into seperate observation bands
        _zbandraw, _ybandraw, _jbandraw, _hbandraw, _kbandraw = lightcurve.splitbands(_data)

    # Plot the bands on the same graph
    plot.plotallbands(_zbandraw, _ybandraw, _jbandraw, _hbandraw, _kbandraw, _period)
    # Plot the black body curve and determine the properties of the binary system
    plot.plotblackbody(_zbandraw, _ybandraw, _jbandraw, _hbandraw, _kbandraw, _parallax)


# Main program entry point
main()
