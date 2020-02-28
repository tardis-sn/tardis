@numba.jit("float(float,float)")
def intensity_black_body_jit(nu,T):
  """
      Calculate the intensity of a black-body according to the following formula
  .. math::
      I(\\nu, T) = \\frac{2h\\nu^3}{c^2}\frac{1}
      {e^{h\\nu \\beta_\\textrm{rad}} - 1}
  Parameters
  ----------
  nu: float
      Frequency of light
  T: float
      Temperature in kelvin
  Returns
  -------
  Intensity: float
      Returns the intensity of the black body
  """
  beta_rad = 1 / (k_B_cgs * T)
  coefficient = 2 * h_cgs / c_cgs ** 2
  intensity = (coefficient * nu**3 )/ 
                          (exp(h_cgs * nu * beta_rad) -1 )
  return intensity
