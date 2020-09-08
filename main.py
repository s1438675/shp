#!/usr/bin/env python
"""
Author: Matthew Christey-Reid
Email:  s1438675@ed.ac.uk
Date:   02/06/2020

main.py - Main entry point of the program
"""

# Import library for passing arguments
import sys

# Import library for parsing FITS files
from astropy.io import fits

# Import libraries for determining properties
import plot
import lightcurve


class Logger(object):
    """
    Provides functionality to write console print statements to an output file
    """
    # Create logger and output file
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("output.txt", "a")

    # When called, write terminal line to file
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    # Dummy function for compatibility
    def flush(self):
        pass


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
        _perr = float(_args[4])
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
    plot.plotblackbody(_zbandraw, _ybandraw, _jbandraw, _hbandraw, _kbandraw, _parallax, _perr)


# Initialise logger
sys.stdout = Logger()

# Main program entry point
main()
