import numpy as np
from rosen98_gasabs import rosen98_gasabs as r98

def gasAbs(frequency, temperature, absolute_humidity, pressure):
	"""
	Provides an interface to the pamtra implementation of Rosenkranz 1998 gas
	absorption model.
	Works only with scalar inputs.

	Parameters
    ----------

    frequency : float
        Electromagnetic frequency [GHz]

    temperature : float
        Atmospheric temperature [Kelvin]

    absolute_humidity : float
    	atmospheric water vapour partial density [kg/m3]

    pressure : float
		atmospheric pressure [Pa]

    Returns
    -------
    abs_air - real
        absorption coeafficient due to dry air (Oxygen and Nitrogen) [Neper/km]

    abs_wv - real
        absorption coeafficient due to water vapour [Neper/km]

    Raises
    ------

    Nothing for now, but it would be better to implement a check on the
    frequency range and validity of temperature, absolute humidity and pressure
	"""
	_, absAir, absWv = r98(frequency, temperature, absolute_humidity, pressure)
	return absAir, absWv