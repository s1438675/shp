import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle


# Take array from phase 0 -> 1, returns a doubled array from 0 -> 2
def doublearrayphase(_inputarray):
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


# Plot all bands on the same graph
def plotcurve(_inputkband, _inputhband, _inputjband, _inputyband, _inputzband):
    _inputkband = doublearrayphase(_inputkband)
    _inputhband = doublearrayphase(_inputhband)
    _inputjband = doublearrayphase(_inputjband)
    _inputyband = doublearrayphase(_inputyband)
    _inputzband = doublearrayphase(_inputzband)


# Plot a single band using pyplot
def plotband(_inputband):
    # Double array to phase 0 -> 2
    _inputband = doublearrayphase(_inputband)
    # Set pyplot style to use
    plt.style.use('seaborn-whitegrid')
    # Plot data as scatter graph with y-axis error bars
    plt.errorbar(_inputband[:, 0], _inputband[:, 1], yerr=_inputband[:, 2], fmt='.k')
    # Set the x-axis label
    plt.xlabel("Phase")
    # Set the y-axis label
    plt.ylabel("Magnitude")
    # Flip y-axis as convention
    plt.gca().invert_yaxis()
    # Display to screen
    plt.show()


# Plot a Lomb-Scargle line of best fit
def plotlobf(_inputband, _period):
    # Create a model with 10 terms
    _ls = LombScargle(_inputband[:, 0], _inputband[:, 1], _inputband[:, 2], nterms=10)
    # Create n linearly spaced points between phase 0 and 1
    _xfit = np.linspace(0, 1, 1000)
    # Frequency = 1 / Period
    _freq = 1 / _period
    # Plot the data in the array to screen, lightly coloured and z rank behind the line of best fit
    plt.errorbar((_inputband[:, 0] * _freq) % 1, _inputband[:, 1], _inputband[:, 2], fmt='.', color='gray',
                 ecolor='lightgray', capsize=0, zorder=0)
    # Plot the line of best fit generated
    _lobf = _ls.model(_xfit / _freq, _freq)
    plt.plot(_xfit, _lobf, '-k', lw=2, zorder=2)
    # Get maxima and minima of the fit
    _max, _min = getcritpoints(_lobf)
    # Plot maxima and minima one the line of best fit
    for i in _min:
        plt.plot(_xfit[i], _lobf[i], '.r', zorder=3)
    for i in _max:
        plt.plot(_xfit[i], _lobf[i], '.b', zorder=3)
    # Invert y-axis as convention
    plt.gca().invert_yaxis()
    # Display to screen
    plt.show()


# Calculate critical points of the line of best fit
def getcritpoints(_inputarray):
    _max, _min = [], []
    if len(_inputarray) < 3:
        return _min, _max

    NEUTRAL, RISING, FALLING = range(3)

    def get_state(a, b):
        if a < b: return RISING
        if a > b: return FALLING
        return NEUTRAL

    ps = get_state(_inputarray[0], _inputarray[1])
    begin = 1
    for i in range(2, len(_inputarray)):
        s = get_state(_inputarray[i - 1], _inputarray[i])
        if s != NEUTRAL:
            if ps != NEUTRAL and ps != s:
                if s == FALLING:
                    _max.append((begin + i - 1) // 2)
                else:
                    _min.append((begin + i - 1) // 2)
            begin = i
            ps = s
    return _min, _max
