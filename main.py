import sys

from astropy.io import fits

import plot
import lightcurve


def main():
    try:
        _args = sys.argv
        _filename = _args[1]
        _period = float(_args[2])
        _parallax = float(_args[3])
    except:
        print("Error loading files, please check your arguments and try again...")
        exit()

    with fits.open(_filename) as _file:
        _data = _file[1].data
        _zbandraw, _ybandraw, _jbandraw, _hbandraw, _kbandraw = lightcurve.splitbands(_data)

    plot.plotallbands(_zbandraw, _ybandraw, _jbandraw, _hbandraw, _kbandraw, _period)
    plot.plotblackbody(_zbandraw, _ybandraw, _jbandraw, _hbandraw, _kbandraw)

main()