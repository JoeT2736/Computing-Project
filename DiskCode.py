import matplotlib.pyplot as plt
import numpy as np
import astropy 
from astropy import constants as const
from astropy import units as u
import scipy
from scipy.integrate import quad
import scipy.integrate

const.L_sun

Nrings = 100
M = 10 * 2 * 10**30  
Mr = 10**15  
Fstart = 10**14
Fstop = 10**19
Fsteps = 100
f = np.linspace(Fstart, Fstop, Fsteps)
Log_f = np.log(f)



Rg = (const.G.value * M) / ((const.c.value)**2)
Rin = 6 * Rg  
Rout = (10**5) * Rg  
R = np.linspace(Rin, Rout, Nrings + 1)
Rmid = (R[:-1] + R[1:]) / 2
Log_RMid = np.log(Rmid)

rin = Rin/Rg
rout = Rout/Rg

r = np.linspace(rin, rout, Nrings + 1) 

r_midpoints = (r[:-1] + r[1:]) / 2

Log_rin = np.log(rin)
Log_rout = np.log(rout)

Log_r = np.linspace(Log_rin, Log_rout, Nrings + 1)
Log_rMid = (Log_r[:-1] + Log_r[1:]) / 2


def Temp():
  T = (((3 * const.G.value * M * Mr) / (8 * np.pi * const.sigma_sb.value * r_midpoints**3 * Rg**3)) * (1 - ( rin / r_midpoints )**(1/2)))**(1/4)
  return T

'''print(Temp())
print(Rmid)

plt.loglog(Rmid, Temp())
plt.title('')
plt.xlabel('')
plt.ylabel('')
plt.show()'''



def Temp_Logs():
  TL = (((3 * const.G.value * M * Mr) / (8 * np.pi * const.sigma_sb.value * Log_rMid**3 * Rg**3)) * (1 - ( Log_rin / Log_rMid )**(1/2)))**(1/4)
  return TL

#print(Temp_Logs())

'''plt.plot(Log_RMid, Temp_Logs())
plt.title('')
plt.xlabel('')
plt.ylabel('')
plt.show()'''


def BlackBodyFlux():
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * Temp_Logs()))) - 1)
  return Fv

#print(BlackBodyFlux())

#plt.plot(Log_rMid, BlackBodyFlux())
#plt.show()

def LV():
  def integrand():
    T = Temp_Logs()
    Fv = BlackBodyFlux()
    return Fv * T * 4 * np.pi * Rg**2
  
  Lv = integrand() * ( Log_rout - Log_rin)
  return Lv

print(LV())



'''plt.loglog(f , f * LV())
plt.show()'''