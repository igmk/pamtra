function e_sat_gg_water(T)
  !    Calculates the saturation pressure over water after Goff and Gratch (1946).
  !    It is the most accurate that you can get for a temperture range from -90degC to +80degC.
  !    Source: Smithsonian Tables 1984, after Goff and Gratch 1946
  !    http://cires.colorado.edu/~voemel/vp.html
  !    http://hurri.kean.edu/~yoh/calculations/satvap/satvap.html
  !    Input: Temperature in Kelvin.

  use kinds
  implicit none

  real(kind=dbl) :: T !T in K
  real(kind=dbl) :: e_sat_gg_water !saturation pressure over water in hPa.

  e_sat_gg_water = 1013.246 * 10**( -7.90298*(373.16/T-1) &
       + 5.02808*log10(373.16/T) &
       - 1.3816e-7*(10**(11.344*(1-T/373.16))-1) &
       + 8.1328e-3 * (10**(-3.49149*(373.16/T-1))-1) )

  return
end function e_sat_gg_water
