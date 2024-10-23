import matplotlib.pyplot as plt
import numpy as np
import astropy 
from astropy import constants as const
from astropy import units as u
import scipy
from scipy.integrate import quad
import scipy.integrate

const.L_sun


M = 10 * const.M_sun.value  #Mass of Black hole
Mr = 10**15  #Accretion rate

Fstart = 10**14
Fstop = 10**19
Fsteps = 100
f = np.linspace(Fstart, Fstop, Fsteps)  #Range of frequencies


Nrings = 1000
Rg = (const.G.value * M) / ((const.c.value)**2)  #
Rin = 6 * Rg  #Innermost stable orbit
Rout = (10**5) * Rg  #Outermost orbit

rin = Rin/Rg  #Scaled innermost stable orbit
rout = Rout/Rg  #Scaled innermost stable orbit

r = np.linspace(rin, rout, Nrings + 1)  
r_midpoints = (r[:-1] + r[1:]) / 2  #Array of increasingly sized disks

Log_rin = np.log(rin)
Log_rout = np.log(rout)

Log_r = np.linspace(Log_rin, Log_rout, Nrings)
Log_rMid = (Log_r[:-1] + Log_r[1:]) / 2  

#Temperature of accretion disk in terms of scaled units
def Temp(M, Mr, r_midpoints, rin):
  T = (((3 * const.G.value * M * Mr) / (8 * np.pi * const.sigma_sb.value * r_midpoints**3 * Rg**3)) * (1 - ( rin / r_midpoints )**(1/2)))**(1/4)
  return T

#Temperature of accretion disk using the log of scaled distances
def Temp_Logs(M, Mr, Log_rMid, Log_rin):
  TL = (((3 * const.G.value * M * Mr) / (8 * np.pi * const.sigma_sb.value * Log_rMid**3 * Rg**3)) * (1 - ( Log_rin / Log_rMid )**(1/2)))**(1/4)
  return TL

#Luminosity per unit frequency per unit area of an isotropically emitting blackbody
def Flux(f, M, Mr, Log_rMid, Log_rin):
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * Temp_Logs(M, Mr, Log_rMid, Log_rin)))) - 1)
  return Fv

#Integrand to calculate luminosity
def integrand(Log_rMid, M, Mr, Log_rin, f):
 
    T = Temp_Logs(M, Mr, Log_rMid, Log_rin)
    Fv = Flux(f, M, Mr, Log_rMid, Log_rin)
    return Fv * 4 * np.pi * Log_rMid * Rg**2  

#Integrating with respect to log(distance) to get luminosity per unit frequency
def Lf(f, M, Mr, Log_rin, Log_rout):
  
    return scipy.integrate.quad(integrand, Log_rin, Log_rout, args=(M, Mr, Log_rin, f))[0]

#Integrating with respect to frequency to get the total luminosity of the accretion disk
def L(M, Mr, Log_rin, Log_rout, Fstart, Fstop):
    
    result, _ = scipy.integrate.quad(Lf, Fstart, Fstop, args=(M, Mr, Log_rin, Log_rout))
    return result

#Luminosity per unit frequency in an array
Lpf = np.array([Lf(freq, M, Mr, Log_rin, Log_rout) for freq in f])


total_L = L(M, Mr, Log_rin, Log_rout, Fstart, Fstop)

print(f"Total Luminosity: {total_L:.3e} W")

'''
plt.plot(r_midpoints, Temp(M, Mr, r_midpoints, rin))
plt.title('Temperature of accretion disk at distance from centre')
plt.show()
'''

'''
plt.plot(Log_rMid, Temp_Logs(M, Mr, Log_rMid, Log_rin))
plt.title('Log(Temperature) of accretion disk at Log(distance) from centre')
plt.show()
'''

plt.loglog(f, f*Lpf)
plt.title('Spectrum of Acctretion disk of black hole')
plt.xlabel('Log(frequency)')
plt.ylabel('Log(frequency*Luminosity)')
plt.grid(True)
plt.show()






