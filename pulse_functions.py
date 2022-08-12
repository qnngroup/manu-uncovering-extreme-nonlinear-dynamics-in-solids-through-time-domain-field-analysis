""" pulse_functions.py
These are a set of convenience functions for creating
pulses of various forms for ultrafast simulations.
"""

import numpy as np

def cos2pulse(t, fwhm, wc, phi_ce):
  """
  cos2pulse(t, fwhm, wc, phi_ce)  

  Produce the electric field (E) and intensity envelope (A) of an ultrafast
  pulse.  Each are normalized to peak of 1) cos^2 shaped pulse. Units are 
  set by the user (i.e. t, fwhm, wc all have to be consistent).  

  Inputs:
  ------------
    t --> time 
    fwhm --> full width at half max
    wc --> central angular frequency (rad/unit time) 
    phi_ce --> CEP (in radians)

  Outputs:
  ---------------
    pulse --> dictionary of pulse field and envelope
      pulse['E'] --> Normalized field (peak of 1)
      pulse['A'] --> Normalized intensity envelope (peak of 1)
  """

  #We can write a cos2 pulse intensity envelope
  #as cos^4(alpha*t).  Then we have that:
  alpha = 2*np.arccos(2**(-0.25))/fwhm

  # We then calculate the envelope as being the cos2 pulse
  # over the time window that we are in the positive half-cycle:
  pulse = {}
  
  pulse['A'] = (abs(t) <= np.pi/2/alpha)*((np.cos(t*alpha))**4)

  #Finally, calculate the electric field:
  pulse['E'] = np.sqrt(pulse['A'])*np.cos(wc*t + phi_ce)

  return pulse

  
  
