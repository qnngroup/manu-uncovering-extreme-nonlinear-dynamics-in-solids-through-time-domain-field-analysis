""" physical_constants_atomic.py
Defines simple physical constants and convenience conversion factors
in atomic units.  All conversions are from SI units (with suffix "con").
See comments below for more specific details.
"""

from numpy import pi

#Conversion factors from si units to atomic units...
lcon = 1/(5.291772192e-11) #atomic length/meter
efieldcon = 1/(5.14220652e11) #atomic e-field/(V/m)
tcon = 1/(2.418884326505e-17) #atomic time/second
evcon = 1.0/27.211 #conversion factor, [J_a/eV]
chargecon = 1/(1.602176565e-19) #atomic charge/Coulomb

#Setup Physical Constants:
hbar = 1.0 #normalized Planck constant [J_a*s_a]
eta = 119.91699832*pi*(1.602176565**2.0)*1e-3/(4.35974417*2.418884326505) #impedance of freespace [J_a*s_a/C_a^2]
me = 1.0 #mass of electron [kg_a]
c = 2.99792458*2.418884326505*1e2/5.291772192 #speed of light [m_a/s_a]
e = 1.0 #electron charge [C_a]
