import pyPamtra

def calc_pyRadarMoments(spectrum,verbose = 0, max_v = 7.885, min_v = -7.885, no_ave = 150, npeaks = 3, noise_distance_factor = 1.25, min_spectral_snr=1.2, smooth_spectrum= True,use_wider_peak=False):
  
  """
  Calculates the moments, slopes and edges of the linear radar spectrum. Refer to the Pamtra documentation for options.

  """

  pyPamtra.pyPamtraLib.report_module.verbose = verbose
  error =   pyPamtra.pyPamtraLib.settings.settings_fill_default()

  pyPamtra.pyPamtraLib.settings.radar_max_v[:]= max_v  
  pyPamtra.pyPamtraLib.settings.radar_min_v[:]= min_v  
  pyPamtra.pyPamtraLib.settings.radar_use_hildebrand= True  
  pyPamtra.pyPamtraLib.settings.radar_no_ave[:]= no_ave  
  pyPamtra.pyPamtraLib.settings.radar_npeaks= npeaks  
  pyPamtra.pyPamtraLib.settings.radar_noise_distance_factor[:]= noise_distance_factor  
  pyPamtra.pyPamtraLib.settings.radar_min_spectral_snr[:]= min_spectral_snr  
  pyPamtra.pyPamtraLib.settings.radar_smooth_spectrum  = smooth_spectrum  
  pyPamtra.pyPamtraLib.settings.radar_use_wider_peak  =  use_wider_peak  
  
  pyPamtra.pyPamtraLib.vars_index.i_f = 1

  noiseModel = 0
  
  error,spectrum_out,moments,slope,edge,quality, noise = pyPamtra.pyPamtraLib.radar_moments.radar_calc_moments(npeaks,spectrum,noiseModel)
  
  return spectrum_out,moments,slope,edge,quality,noise, pyPamtra.pyPamtraLib, error


def calc_hildebrandSekhon(spectrum, no_ave = 1,verbose=0):

  """
  Calculate the mean and maximum of noise of teh linear radar spectrum following Hildebrand and Sekhon 1974.

  Parameters
  ----------

  spectrum : array
      linear radar spectrum
  no_ave : int, optional
      number of averages (default 1)
  verbose : int, optional
      verbosity level (default 0)

  Returns
  -------

  meanNoise : float
    mean noise level in linear units
  maxNoise : float
    maximum noise level in linear units
  error : int
    Pamtra error code. Values larger than zero indicate an error
  """
  pyPamtra.pyPamtraLib.report_module.verbose = verbose
  error, meanNoise, maxNoise = pyPamtra.pyPamtraLib.radar_hildebrand_sekhon(spectrum,no_ave)
  return meanNoise, maxNoise, error