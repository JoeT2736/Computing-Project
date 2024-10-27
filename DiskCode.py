import matplotlib.pyplot as plt
import numpy as np
import astropy 
from astropy import constants as const
from astropy import units as u
import scipy
from scipy.integrate import quad
import scipy.integrate


M = 10 * const.M_sun.value  #Mass of Black hole
Mr = 10**15  #Accretion rate
Mr2 = 10**14

Fstart = 10**14
Fstop = 10**19
Fsteps = 1000
f = np.linspace(Fstart, Fstop, Fsteps)  #Range of frequencies


Nrings = 1000
Rg = (const.G.value * M) / ((const.c.value)**2)  #Schwarzschild radius
Rin = 6 * Rg  #Innermost stable orbit
Rin2 = 1.2 * Rg  #Innermost stable orbit max spin
Rout = (10**5) * Rg  #Outermost orbit
R = np.linspace(Rin, Rout, Nrings + 1) 
R2 = np.linspace(Rin2, Rout, Nrings + 1) 
R_midpoints = (R[:-1] + R[1:]) / 2
R_midpoints2 = (R2[:-1] + R2[1:]) / 2

rin = Rin/Rg  #Scaled innermost stable orbit
rin2 = Rin2/Rg  #Scaled innermost stable orbit
rout = Rout/Rg  #Scaled innermost stable orbit

r = np.linspace(rin, rout, Nrings + 1)  
r_midpoints = (r[:-1] + r[1:]) / 2  #Array of increasingly sized disks
r2 = np.linspace(rin2, rout, Nrings + 1)  
r_midpoints2 = (r2[:-1] + r2[1:]) / 2 

Log_rin = np.log(rin)
Log_rout = np.log(rout)

Log_r = np.linspace(Log_rin, Log_rout, Nrings)
Log_rMid = (Log_r[:-1] + Log_r[1:]) / 2  

Log_rin2 = np.log(rin2)
Log_r2 = np.linspace(Log_rin2, Log_rout, Nrings)
Log_rMid2 = (Log_r2[:-1] + Log_r2[1:]) / 2  


def Temp1(M, Mr, R_midpoints, Rin):
  T = (((3 * const.G.value * M * Mr) / (8 * np.pi * const.sigma_sb.value * R_midpoints**3 * Rg**3)) * (1 - ( Rin / R_midpoints )**(1/2)))**(1/4)
  return T


#Temperature of accretion disk in terms of scaled units
def Temp(M, Mr, r_midpoints, rin):
  T = (((3 * const.G.value * M * Mr) / (8 * np.pi * const.sigma_sb.value * r_midpoints**3 * Rg**3)) * (1 - ( rin / r_midpoints )**(1/2)))**(1/4)
  return T

#Temperature of accretion disk using the log of scaled distances
def Temp_Logs(M, Mr, Log_rMid, Log_rin):
  TL = (((3 * const.G.value * M * Mr) / (8 * np.pi * const.sigma_sb.value * Log_rMid**3 * Rg**3)) * (1 - ( Log_rin / Log_rMid )**(1/2)))**(1/4)
  return TL

#different r range
def Temp_Logs2(M, Mr, Log_rMid2, Log_rin2):
  TL = (((3 * const.G.value * M * Mr) / (8 * np.pi * const.sigma_sb.value * Log_rMid2**3 * Rg**3)) * (1 - ( Log_rin2 / Log_rMid2 )**(1/2)))**(1/4)
  return TL

#Different accretion rate
def Temp_Logs3(M, Mr2, Log_rMid, Log_rin):
  TL = (((3 * const.G.value * M * Mr2) / (8 * np.pi * const.sigma_sb.value * Log_rMid**3 * Rg**3)) * (1 - ( Log_rin / Log_rMid )**(1/2)))**(1/4)
  return TL

#Luminosity per unit frequency per unit area of an isotropically emitting blackbody
def Flux(f, M, Mr, Log_rMid, Log_rin):
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * Temp_Logs(M, Mr, Log_rMid, Log_rin)))) - 1)
  return Fv

def Flux2(f, M, Mr, Log_rMid2, Log_rin2):
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * Temp_Logs2(M, Mr, Log_rMid2, Log_rin2)))) - 1)
  return Fv

def Flux3(f, M, Mr2, Log_rMid, Log_rin):
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * Temp_Logs3(M, Mr2, Log_rMid, Log_rin)))) - 1)
  return Fv

#Integrand to calculate luminosity
def integrand(Log_rMid, M, Mr, Log_rin, f):
 
    T = Temp_Logs(M, Mr, Log_rMid, Log_rin)
    Fv = Flux(f, M, Mr, Log_rMid, Log_rin)
    return Fv * 4 * np.pi * Log_rMid * Rg**2  

def integrand2(Log_rMid2, M, Mr, Log_rin2, f):
 
    T = Temp_Logs2(M, Mr, Log_rMid2, Log_rin2)
    Fv = Flux2(f, M, Mr, Log_rMid2, Log_rin2)
    return Fv * 4 * np.pi * Log_rMid2 * Rg**2  


def integrand3(Log_rMid, M, Mr2, Log_rin, f):
 
    T = Temp_Logs3(M, Mr2, Log_rMid, Log_rin)
    Fv = Flux3(f, M, Mr2, Log_rMid, Log_rin)
    return Fv * 4 * np.pi * Log_rMid * Rg**2  


#Integrating with respect to log(distance) to get luminosity per unit frequency
def Lf(f, M, Mr, Log_rin, Log_rout):
  
    return scipy.integrate.quad(integrand, Log_rin, Log_rout, args=(M, Mr, Log_rin, f))[0]

def Lf2(f, M, Mr, Log_rin2, Log_rout):
  
    return scipy.integrate.quad(integrand2, Log_rin2, Log_rout, args=(M, Mr, Log_rin2, f))[0]

def Lf3(f, M, Mr2, Log_rin, Log_rout):
  
    return scipy.integrate.quad(integrand3, Log_rin, Log_rout, args=(M, Mr2, Log_rin, f))[0]

#Integrating with respect to frequency to get the total luminosity of the accretion disk
def L(M, Mr, Log_rin, Log_rout, Fstart, Fstop):
    
    result, _ = scipy.integrate.quad(Lf, Fstart, Fstop, args=(M, Mr, Log_rin, Log_rout))
    return result

def L2(M, Mr, Log_rin2, Log_rout, Fstart, Fstop):
    
    result, _ = scipy.integrate.quad(Lf2, Fstart, Fstop, args=(M, Mr, Log_rin2, Log_rout))
    return result

def L3(M, Mr2, Log_rin, Log_rout, Fstart, Fstop):
    
    result, _ = scipy.integrate.quad(Lf3, Fstart, Fstop, args=(M, Mr2, Log_rin, Log_rout))
    return result


#Luminosity per unit frequency in an array
Lpf = np.array([Lf(freq, M, Mr, Log_rin, Log_rout) for freq in f])
Lpf2 = np.array([Lf2(freq, M, Mr, Log_rin2, Log_rout) for freq in f])
Lpf3 = np.array([Lf3(freq, M, Mr2, Log_rin, Log_rout) for freq in f])

total_L = L(M, Mr, Log_rin, Log_rout, Fstart, Fstop)
total_L2 = L(M, Mr, Log_rin2, Log_rout, Fstart, Fstop)
total_L3 = L(M, Mr2, Log_rin, Log_rout, Fstart, Fstop)

print(f"Total Luminosity for no spin, 10e15Mr: {total_L:.3e} W")
print(f"Total Luminosity for max spin: {total_L2:.3e} W")
print(f"Total Luminosity for no spin, 10e13Mr: {total_L3:.3e} W")

'''
plt.plot(R_midpoints, Temp1(M, Mr, R_midpoints, Rin))
plt.xlim(-900000, 0.06e9)
plt.title('Temperature of accretion disk at distance from innermost stable orbit')
plt.xlabel('Distance from innermost stable orbit (m)')
plt.ylabel('Temperature (K)')
plt.show()
'''

'''
plt.plot(r_midpoints, Temp(M, Mr, r_midpoints, rin))
plt.xlim(-50, 3500)
plt.title('Temperature of accretion disk at distance from innermost stable orbit (in terms of scaled units)')
plt.xlabel('Distance from innermost stable orbit')
plt.ylabel('Temperature (K)')
plt.show()
'''

'''
plt.plot(Log_rMid, Temp_Logs(M, Mr, Log_rMid, Log_rin))
plt.title('Log(Temperature) of accretion disk at Log(distance) from centre')
plt.show()
'''

plt.loglog(f, f*Lpf, label='Mr=10e15, Rin=6Rg')
plt.loglog(f, f*Lpf2, label='Mr=10e15, Rin=1.2Rg')
plt.loglog(f, f*Lpf3, label='Mr=10e14, Rin=6Rg')
plt.legend()
plt.title('Spectrum of Acctretion disk of black hole')
plt.xlabel('Log(frequency)')
plt.ylabel('Log(frequency*Luminosity)')
plt.grid(True)
plt.show()






