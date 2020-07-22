# Import libraries for data processing
import math
import numpy as np
import matplotlib.pyplot as plt
from math import log10, floor
from scipy.optimize import curve_fit

# Import lightcurve functions
import lightcurve

# Import Lomb-Scargle periodogram functions
from astropy.timeseries import LombScargle

# Arrays to store common values for each observed band
VegaToAB = [0.528, 0.634, 0.938, 1.379, 1.9]  # Constant for conversion between Vega and AB magnitude systems
Wavelength = [0.9, 1.02, 1.25, 1.6, 2.2]      # Effective wavelength of observation band


def round_sig(_val):
    """
    Rounds the provided value to 2 significant figures
    :param _val: Value to be rounded
    :return: Float, original value rounded to 2 significant figures
    """
    return round(_val, 3 - int(floor(log10(abs(_val)))) - 1)


def plotallbands(_zband, _yband, _jband, _hband, _kband, _period):
    """
    Plots all observed bands to the same graph
    :param _zband: Observed z-band
    :param _yband: Observed y-band
    :param _jband: Observed j-band
    :param _hband: Observed h-band
    :param _kband: Observed k-band
    :param _period: Period of variability
    """
    # Set pyplot style to be consisten within the program
    plt.style.use('seaborn-whitegrid')
    # Frequency = 1 / Period
    _freq = 1 / _period

    # Create single dataset from all bands
    _bands = [_zband, _yband, _jband, _hband, _kband]
    # Iterate through each band and plot to screen
    i = 0
    while i < 5:
        # Array to set colours for each band
        _colours = ['-b', '-g', '-r', '-c', '-m']
        # Array to set strings for graph legend
        _legend = ['Z-band', 'Y-band', 'J-band', 'H-band', 'K-band']
        # Determine the line of best fit for each band
        _xfit, _lobf = calclobf(_bands[i], _period)
        # Plot the data in the array to screen, lightly coloured and z rank behind the line of best fit
        plt.plot(_xfit, _lobf, _colours[i], lw=1, zorder=2, label=_legend[i])
        i += 1

    # Set x-axis limit to a single period
    plt.xlim(0, 1)
    # Set graph and axis titles
    plt.xlabel("Phase")
    plt.ylabel("Magnitude")
    plt.title("Folded light curve")
    # Show the legend
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
    plt.style.use('seaborn-whitegrid')
    plt.errorbar((_band[:, 0] * _freq) % 1, _band[:, 1], _band[:, 2], fmt='.', color='gray',
                 ecolor='lightgray', capsize=0, zorder=0)
    # Plot the graph of the line of best fit
    plt.plot(_xfit, _lobf, '-k', lw=2, zorder=2)
    # Set x-axis limits to 1 period
    plt.xlim(0, 1)
    # Set graph and axis titles
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


def plotblackbody(_zband, _yband, _jband, _hband, _kband, _parallax):
    """
    Determines the black body curve and determines mass, radius and luminosities in solar units
    :param _zband: Observed z-band
    :param _yband: Observed y-band
    :param _jband: Observed j-band
    :param _hband: Observed h-band
    :param _kband: Observed k-band
    :param _parallax: Parallax angle (mas)
    """
    # Set pyplot style to be consistent within the program
    plt.style.use('seaborn-whitegrid')
    # Import raw data to plot Hertzsprung-Russell diagram
    _hrdata = inithr('hr.dat')
    # Determine distance in parsecs
    _distance = 1 / np.tan(_parallax * 10**-3)
    # Create single data array with all bands
    _bands = [_zband, _yband, _jband, _hband, _kband]
    _lambda = [0.9, 1.02, 1.22, 1.63, 2.2]
    # Set up empty arrays for each star
    _largestar = np.zeros((1, 2))
    _smallstar = np.zeros((1, 2))

    # Determine the spectral flux density from the large star
    i = 0
    while i < 5:
        # Determine the maximum and minimum values of the observed band
        _max, _min = lightcurve.maxminvals(_bands[i])
        # The large star uses the maximum flux value (smallest magnitude)
        _largestar = np.append(_largestar, np.array([_lambda[i], (magtoflux(_min, i))], ndmin=2), axis=0)
        i += 1
    # Delete first empty row of the array
    _largestar = np.delete(_largestar, 0, axis=0)

    # Determine the spectral flux density from the small star
    i = 0
    while i < 5:
        # Determine the maximum and minimum values of the observed band
        _max, _min = lightcurve.maxminvals(_bands[i])
        # Smaller star flux value is combined value minus the large star
        _smallstar = np.append(_smallstar, np.array([_lambda[i], (magtoflux(_max, i) -
                                                                  magtoflux(_min, i))], ndmin=2), axis=0)
        i += 1
    # Delete the first empty row of the array
    _smallstar = np.delete(_smallstar, 0, axis=0)

    # Determine the luminosity and effective temperature of each star
    _luma, _wiena = getwientemp(_largestar, _distance, 1)
    _lumb, _wienb = getwientemp(_smallstar, _distance, 2)

    # Calculate luminosities in solar units
    _solluma = _luma / (3.828*10**26)
    _sollumb = _lumb / (3.828*10**26)

    # Calculate masses using the mass/luminosity relation in solar mass units
    # N.B. only works as an approximation for main sequence stars, giants and dwarfs are not sutiable for this
    # approximation
    _solmassa = np.power(_solluma, 1/3.5)
    _solmassb = np.power(_sollumb, 1/3.5)

    # Calculate stellar radius in solar radii using the relationship between luminosity, surface area and temperature
    _solrada = np.sqrt(_solluma / np.power(_wiena / 5778, 4))
    _solradb = np.sqrt(_sollumb / np.power(_wienb / 5778, 4))

    # Output determined values to the screen
    print('Values for the large star:')
    print('Effective temperature: ' + str(round_sig(_wiena)))
    print('Solar luminosities: ' + str(round_sig(_solluma)))
    print('Solar radii: ' + str(round_sig(_solrada)))
    print('Solar masses: ' + str(round_sig(_solmassa)))
    print('-----------------------------------------------------')
    print('Values for the small star:')
    print('Effective temperature: ' + str(round_sig(_wienb)))
    print('Solar luminosities: ' + str(round_sig(_sollumb)))
    print('Solar radii: ' + str(round_sig(_solradb)))
    print('Solar masses: ' + str(round_sig(_solmassb)))

    # Convert from luminosity to magnitude in solar units
    _luma = -2.5 * np.log10(_luma / (3.0128 * 10**28))
    _lumb = -2.5 * np.log10(_lumb / (3.0128 * 10**28))

    # Plot Hertzsprung-Russell diagram using provided array
    plt.scatter(_hrdata[:, 1], _hrdata[:, 0], s=0.5)
    # Plot determined values for each star
    plt.scatter(_wiena, _luma, s=16, c='red')
    plt.scatter(_wienb, _lumb, s=16, c='pink')
    # Set the x and y axis limits to sensible values
    plt.xlim(3000, 10000)
    plt.ylim(-10, 20)
    # Invert both axes as convention
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    # Display to screen
    plt.show()


def getwientemp(_inputdata, _distance, _id):
    """
    Determines the effective temperature using Wien's law
    :param _inputdata: Black body curve of object
    :param _distance: Distance to object (parsecs)
    :param _id: 1 for large star, 2 for small star
    :return: Luminosity and effective surface temperature
    """
    # Maxwell-Boltzmann distribution formula probability density function
    def curve(_x, _a, _scale):
        _a1 = np.sqrt(2 / np.pi)
        _a2 = _x**2 / (2 * _a**2)
        return _scale * _a1 * (_x**2 * np.exp(-_a2)) / _a**3

    # Set pyplot style to be consistent through the program
    plt.style.use('seaborn-whitegrid')

    # Convert the distance in parsecs to metres
    _distance = 3.0857 * 10**16 * _distance
    # Create array for x and y axis data
    _xdata = _inputdata[:, 0]
    _ydata = _inputdata[:, 1]
    _ydatalum = _ydata

    # Iterate through each band and convert from Janskys to W/m^2/um
    i = 0
    while i < 5:
        _ydata[i] = 3*10**14 * (_ydata[i] * 10**-26) / (Wavelength[i]**2)
        i += 1
    # Calculate optimal values and covariance using scipy curve_fit function
    _popt, _pcov = curve_fit(curve, _xdata, _ydata)
    # Create x axis to plot curve against
    _x = np.linspace(0, 5, 100)
    # Determine y value for each point on the x axis
    _yplot = curve(_x, *_popt)
    # Plot the curve to the screen
    plt.plot(_x, _yplot)
    # Determine the area under the graph, integral gives total energy recieved per m^2
    _area = np.trapz(_yplot, dx=5/100)
    # Total luminosity found by multiplying by the surface area of a sphere with diameter of the distance
    _lum = 4 * np.pi * _distance**2 * _area
    # Peak value of Maxwell-Boltzmann distribution
    _mu = 2 * _popt[0] * np.sqrt(2 / np.pi)

    # Plot data on the graph
    plt.plot(_xdata, _ydata, '.')
    # Set axis labels
    plt.xlabel('Wavelength (um)')
    plt.ylabel('Spectral Irradiance (W m^-2 um^-1)')
    if _id == 1:
        _str = 'Large Star'
    else:
        _str = 'Small Star'

    # Calculate effective surface temperature using Wien's law
    _wien = round_sig(2898 / _mu)
    # Round luminosity to 2 significant figures
    _lum = round_sig(_lum)
    # Set graph title
    plt.suptitle('Black Body Plot for the ' + _str)
    # Display to the screen
    plt.show()

    # Returns calculated values
    return _lum, _wien


def inithr(_filename):
    """
    Parses required data for plotting a Hertzsprung-Russell diagram
    :param _filename: File containing observed data
    :return: (n x 3) size array containing magnitude, effective temperature and parallax angle
    """
    # Open file provided
    _file = open(_filename)
    # Create empty array to hold data
    _data = np.zeros((1, 3), dtype=float)

    # Iterate through the file line by line
    for _line in _file:
        # Split each line into constituent values
        _x = _line.split()
        # Append data array with each value, converted to float, convert parallax angle to distance
        _data = np.append(_data, np.array([float(_x[1]), float(_x[2]), (1 / float(_x[3]))], ndmin=2), axis=0)

    # Iterate through data array
    for _row in _data:
        np.seterr(divide='ignore')
        # Convert magnitude to luminosity
        _row[0] = _row[0] - 5 * (np.log10(_row[2]) - 1)
        # Convert B-V colour to temperature
        _row[1] = 4600 * ((1 / (0.92 * _row[1] + 1.7)) + 1 / (0.92 * _row[1] + 0.62))

    # Delete first empty row
    _data = np.delete(_data, 0, axis=0)

    # Return parsed data
    return _data


def magtoflux(_mag, _id):
    """
    Converts magnitude to flux in Janskys
    :param _mag: Magnitude of object
    :param _id: ID of observation band
    :return: Spectral flux density in Janskys
    """
    return math.pow(10, -0.4*(_mag + VegaToAB[_id] - 8.9))
