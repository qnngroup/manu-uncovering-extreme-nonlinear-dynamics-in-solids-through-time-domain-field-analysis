""" FN_functions.py
These are a set of convenience functions for modeling of FN
emission from surfaces. 

Functions that offer the derivative of the emission rate are also provided
for convenience when performing small-signal sampling analysis.  
"""

import numpy as np

def J_FN_atomic(F, phi):
  """
  J_FN_atomic(F, phi)
   
  Description
  ------------------
  Calculates the Fowler-Nordheim emission rate as a function 
  of the work function (phi) and field strength (F) at a surface.  
  Note that the emission is in arbitrary units and is not intended 
  to provide absolute current density output.  
  
  Inputs
   -----------
  F --> electric field strength (atomic units)
  phi --> work function (energy in atomic units)
  
  Outputs
  -------------
  Current density (relative current density...i.e. to within a constant factor).
  """


  F_valid_range = np.where(F < 0)
  F_valid = F[F_valid_range]
  
  J = np.zeros(F.shape)
  J[F_valid_range] = F_valid**2*np.exp(-4*np.sqrt(2)*phi**1.5/(-3*F_valid))

  return J

def J_FN_SI(F, phi):
  """
  J_FN_SI(

  Simplified Standard Fowler-Nordheim-type equation for calculating current density.
  
  The point of this function is to provide output in SI units. As such, all inputs
  are expected in SI units as specified below.

  Likewise, the current density is returned in units of A/nm^2.
 
  Inputs
  ---------
    F -- field (V/nm)
    phi -- work function (eV)

  Outputs
  ---------
    Current density (A/nm^2)
  """
    
  # Calculates physical current density
  # https://en.wikipedia.org/wiki/Field_electron_emission#Recommended_form_for_simple_Fowler%E2%80%93Nordheim-type_calculations

  J = np.zeros(F.shape)

  #We need to invert F here as electron emission only occurs for negative fields.
  #The following formulas were defined assuming F is positive to emit the particle.
  F = -1*F

  #The FN function is only defined over regions where the field is nonzero.  We now
  #find the indices of those regions and only perform our calculations in those
  #regions.  Everywhere else will be left as 0.
  F_valid_range = np.where(F > 0)
  F_valid = F[F_valid_range]

  a_const = 1.541534e-6
  b_const = 6.83089  

  f = 1.43996453529595*F_valid/phi**2
  v_f = 1 - f + 1/6*f*np.log(abs(f))

  J[F_valid_range] = a_const/phi*F_valid**2*np.exp(-v_f*b_const*phi**(3/2)/F_valid)

  return J

def dJ_dF_FN_SI(F, phi):
  """
  dJ_dF_FN_SI(

  Derivative of Simplified Standard Fowler-Nordheim-type equation with respect to field.  
  
  This function is useful for determining the true transfer function of the emitter in the context of small-signal analysis.  
  
  The point of this function is to provide output in SI units. As such, all inputs
  are expected in SI units as specified below.
 
  Inputs
  ---------
    F -- field (V/nm)
    phi -- work function (eV)

  Outputs
  ---------
    Derivative of current density with respect to field (A/V/nm)
  """
    
  # Calculates physical current density
  # https://en.wikipedia.org/wiki/Field_electron_emission#Recommended_form_for_simple_Fowler%E2%80%93Nordheim-type_calculations
  dJ_dF = np.zeros(F.shape)

  #We need to invert F here as electron emission only occurs for negative fields.
  #The following formulas were defined assuming F is positive to emit the particle.
  F = -1*F

  #The FN function is only defined over regions where the field is nonzero.  We now
  #find the indices of those regions and only perform our calculations in those
  #regions.  Everywhere else will be left as 0.
  F_valid_range = np.where(F > 0)
  F_valid = F[F_valid_range]

  a_const = 1.541534e-6
  b_const = 6.83089  

  f = 1.43996453529595*F_valid/phi**2
  v_f = 1 - f + 1/6*f*np.log(abs(f))

  alpha = a_const/phi
  beta = v_f*b_const*phi**(3/2)

  dJ_dF[F_valid_range] = (2*a_const/phi*F_valid + a_const*v_f*b_const*phi**(1/2))*np.exp(-v_f*b_const*phi**(3/2)/F_valid)
    

  return dJ_dF



