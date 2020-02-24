import numpy as np
import matplotlib.pyplot as plt


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


def plotlobf(_inputband):
    _p = np.polynomial.Polynomial.fit(_inputband[:, 0], _inputband[:, 1], 30)
    plt.plot(*_p.linspace())
    plt.gca().invert_yaxis()
    plt.show()
