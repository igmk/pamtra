import pyRadarMoments

def calc_pyRadarMoments(spectrum,verbose = 0, max_v = 7.885, min_v = -7.885, no_ave = 150, npeaks = 3, noise_distance_factor = 1.25, min_spectral_snr=1.2, smooth_spectrum= True):
  
  pyRadarMoments.report_module.verbose = verbose
  error = pyRadarMoments.settings.settings_fill_default()

  pyRadarMoments.settings.radar_max_v= max_v  
  pyRadarMoments.settings.radar_min_v= min_v  
  pyRadarMoments.settings.radar_use_hildebrand= True  
  pyRadarMoments.settings.radar_no_ave= no_ave  
  pyRadarMoments.settings.radar_npeaks= npeaks  
  pyRadarMoments.settings.radar_noise_distance_factor= noise_distance_factor  
  pyRadarMoments.settings.radar_min_spectral_snr= min_spectral_snr  
  pyRadarMoments.settings.radar_smooth_spectrum  = smooth_spectrum  
  
  noiseModel = 0
  
  error,spectrum_out,moments,slope,edge,quality, noise = pyRadarMoments.radar_moments.radar_calc_moments(npeaks,spectrum,noiseModel)
  
  return spectrum_out,moments,slope,edge,quality,noise, pyRadarMoments, error


