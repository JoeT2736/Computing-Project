import matplotlib.pyplot as plt
import numpy as np
import astropy 
from astropy import constants as const
from astropy import units as u
import scipy
from scipy.integrate import quad
import scipy.integrate
import sympy



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

Log_r = np.linspace(Log_rin, Log_rout, Nrings + 1)
#Log_rMid = (Log_r[:-1] + Log_r[1:]) / 2  

Log_rin2 = np.log(rin2)
Log_r2 = np.linspace(Log_rin2, Log_rout, Nrings)
Log_rMid2 = (Log_r2[:-1] + Log_r2[1:]) / 2 

#min/max efficiency of accretion rate to energy
eta_min = 0.1
eta_max = 0.4

#Eddington limit
def L_edd(M):
  return (4 * np.pi * const.G.value * M * const.m_p.value * const.c.value)/const.sigma_T.value

#Accretion rate limit for min efficiency
Mr_edd_min = L_edd(M)/(eta_max * const.c.value**2)

#Accretion rate limit for max efficiency
Mr_edd_max = L_edd(M)/(eta_min * const.c.value**2)



'''
def Temp1(M, Mr, R_midpoints, Rin):
  T = (((3 * const.G.value * M * Mr) / (8 * np.pi * const.sigma_sb.value * R_midpoints**3)) * (1 - ( Rin / R_midpoints )**(1/2)))**(1/4)
  return T



plt.plot(R_midpoints, Temp1(M, Mr, R_midpoints, Rin))
plt.xlim(-900000, 0.06e9)
plt.title('Temperature of accretion disk at distance from innermost stable orbit')
plt.xlabel('Distance from innermost stable orbit (m)')
plt.ylabel('Temperature (K)')
plt.show()



#Temperature of accretion disk in terms of scaled units
def Temp(M, Mr, Log_rMid, Log_rin):
  T = (((3 * const.G.value * M * Mr) / (8 * np.pi * const.sigma_sb.value * Log_rMid**3 * Rg**3)) * (1 - ( Log_rin / Log_rMid )**(1/2)))**(1/4)
  return T

def Flux(f, M, Mr, Log_rMid, Log_rin):
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * Temp(M, Mr, Log_rMid, Log_rin)))) - 1)
  return Fv

def integrand(Log_rMid, M, Mr, f, Rg):
  T = (((3 * const.G.value * M * Mr) / (8 * np.pi * const.sigma_sb.value * Log_rMid**3 * Rg**3)) * (1 - ( Log_rin / Log_rMid )**(1/2)))**(1/4)
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * T))) - 1)
  return Fv * T * 4 * np.pi * Log_rMid * Rg**2 


def Lf(f, Log_rin, Log_rout):
  return scipy.integrate.quad(integrand, Log_rin, Log_rout, args=(M, Mr, f, Rg))[0]


Lpf = np.array([Lf(freq, Log_rin, Log_rout) for freq in f1])
#print(Lpf)

def L(Fstart, Fstop):
  result = scipy.integrate.quad(Lf, Fstart, Fstop, args=(Log_rin, Log_rout))[0]
  return result

Lt = scipy.integrate.trapz(f1*Lpf, f1)


total_L = L(Fstart, Fstop)
print(f"Total Luminosity for no spin, 10e15Mr: {total_L:.3e} W")
print(Lt)

plt.loglog(f1, f1*Lpf, label='Mr=10e15, Rin=6Rg')
plt.title('Spectrum of Acctretion disk of black hole')
plt.xlabel('Log(frequency)')
plt.ylabel('Log(frequency*Luminosity)')
plt.grid(True)
plt.show()
'''





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

#At Eddington limit min efficiency
def Temp_Logs4(M, Mr_edd_min, Log_rMid, Log_rin):
  TL = (((3 * const.G.value * M * Mr_edd_min) / (8 * np.pi * const.sigma_sb.value * Log_rMid**3 * Rg**3)) * (1 - ( Log_rin / Log_rMid )**(1/2)))**(1/4)
  return TL

#At Eddington limit max efficiency
def Temp_Logs5(M, Mr_edd_max, Log_rMid, Log_rin):
  TL = (((3 * const.G.value * M * Mr_edd_max) / (8 * np.pi * const.sigma_sb.value * Log_rMid**3 * Rg**3)) * (1 - ( Log_rin / Log_rMid )**(1/2)))**(1/4)
  return TL


#Luminosity per unit frequency per unit area of an isotropically emitting blackbody
def Flux(f, M, Mr, Log_rMid, Log_rin):
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * Temp_Logs(M, Mr, Log_rMid, Log_rin)))) - 1)
  return Fv

#different r range
def Flux2(f, M, Mr, Log_rMid2, Log_rin2):
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * Temp_Logs2(M, Mr, Log_rMid2, Log_rin2)))) - 1)
  return Fv

#different accretion rate
def Flux3(f, M, Mr2, Log_rMid, Log_rin):
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * Temp_Logs3(M, Mr2, Log_rMid, Log_rin)))) - 1)
  return Fv


def Flux4(f, M, Mr_edd_min, Log_rMid, Log_rin):
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * Temp_Logs4(M, Mr_edd_min, Log_rMid, Log_rin)))) - 1)
  return Fv


def Flux5(f, M, Mr_edd_max, Log_rMid, Log_rin):
  Fv = ((2 * np.pi * const.h.value * f**3)/(const.c.value**2))/((np.exp(((const.h.value * f))/(const.k_B.value * Temp_Logs5(M, Mr_edd_max, Log_rMid, Log_rin)))) - 1)
  return Fv


#Integrand to calculate luminosity
def integrand(Log_rMid, M, Mr, Log_rin, f):
 
    T = Temp_Logs(M, Mr, Log_rMid, Log_rin)
    Fv = Flux(f, M, Mr, Log_rMid, Log_rin)
    return Fv * T * 4 * np.pi * Log_rMid * Rg**2  

#different r range
def integrand2(Log_rMid2, M, Mr, Log_rin2, f):
 
    T = Temp_Logs2(M, Mr, Log_rMid2, Log_rin2)
    Fv = Flux2(f, M, Mr, Log_rMid2, Log_rin2)
    return Fv * T * 4 * np.pi * Log_rMid2 * Rg**2  

#different accretion rate
def integrand3(Log_rMid, M, Mr2, Log_rin, f):
 
    T = Temp_Logs3(M, Mr2, Log_rMid, Log_rin)
    Fv = Flux3(f, M, Mr2, Log_rMid, Log_rin)
    return Fv * T * 4 * np.pi * Log_rMid * Rg**2  


def integrand4(Log_rMid, M, Mr_edd_min, Log_rin, f):
 
    T = Temp_Logs4(M, Mr_edd_min, Log_rMid, Log_rin)
    Fv = Flux(f, M, Mr_edd_min, Log_rMid, Log_rin)
    return Fv * T * 4 * np.pi * Log_rMid * Rg**2  


def integrand5(Log_rMid, M, Mr_edd_max, Log_rin, f):
 
    T = Temp_Logs5(M, Mr_edd_max, Log_rMid, Log_rin)
    Fv = Flux(f, M, Mr_edd_max, Log_rMid, Log_rin)
    return Fv * T * 4 * np.pi * Log_rMid * Rg**2  



#Integrating with respect to log(distance) to get luminosity per unit frequency
def Lf(f, M, Mr, Log_rin, Log_rout):
  
    return scipy.integrate.quad(integrand, Log_rin, Log_rout, args=(M, Mr, Log_rin, f))[0]

#different r range
def Lf2(f, M, Mr, Log_rin2, Log_rout):
  
    return scipy.integrate.quad(integrand2, Log_rin2, Log_rout, args=(M, Mr, Log_rin2, f))[0]

#different accretion rate
def Lf3(f, M, Mr2, Log_rin, Log_rout):
  
    return scipy.integrate.quad(integrand3, Log_rin, Log_rout, args=(M, Mr2, Log_rin, f))[0]


def Lf4(f, M, Mr_edd_min, Log_rin, Log_rout):
  
    return scipy.integrate.quad(integrand4, Log_rin, Log_rout, args=(M, Mr_edd_min, Log_rin, f))[0]


def Lf5(f, M, Mr_edd_max, Log_rin, Log_rout):
  
    return scipy.integrate.quad(integrand5, Log_rin, Log_rout, args=(M, Mr_edd_max, Log_rin, f))[0]


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


def L4(M, Mr_edd_min, Log_rin, Log_rout, Fstart, Fstop):
    
    result, _ = scipy.integrate.quad(Lf4, Fstart, Fstop, args=(M, Mr_edd_min, Log_rin, Log_rout))
    return result


def L5(M, Mr_edd_max, Log_rin, Log_rout, Fstart, Fstop):
    
    result, _ = scipy.integrate.quad(Lf5, Fstart, Fstop, args=(M, Mr_edd_max, Log_rin, Log_rout))
    return result


#Luminosity per unit frequency in an array
Lpf = np.array([Lf(freq, M, Mr, Log_rin, Log_rout) for freq in f])
Lpf2 = np.array([Lf2(freq, M, Mr, Log_rin2, Log_rout) for freq in f])
Lpf3 = np.array([Lf3(freq, M, Mr2, Log_rin, Log_rout) for freq in f])
Lpf4 = np.array([Lf4(freq, M, Mr_edd_min, Log_rin, Log_rout) for freq in f])
Lpf5 = np.array([Lf5(freq, M, Mr_edd_max, Log_rin, Log_rout) for freq in f])

total_L = L(M, Mr, Log_rin, Log_rout, Fstart, Fstop)
total_L2 = L2(M, Mr, Log_rin2, Log_rout, Fstart, Fstop)
total_L3 = L3(M, Mr2, Log_rin, Log_rout, Fstart, Fstop)
total_L4 = L(M, Mr_edd_min, Log_rin, Log_rout, Fstart, Fstop)
total_L5 = L(M, Mr_edd_max, Log_rin, Log_rout, Fstart, Fstop)

print(f"Eddington limit of Luminosity for 10 solar mass: {L_edd(M):.3e} W")
print(f"Total Luminosity for no spin, min eddington limit accretion rate: {total_L4:.3e} W")
print(f"Total Luminosity for no spin, max eddington limit accretion rate: {total_L5:.3e} W")
print(f"Total Luminosity for no spin, 10e15Mr: {total_L:.3e} W")
print(f"Total Luminosity for max spin: {total_L2:.3e} W")
print(f"Total Luminosity for no spin, 10e13Mr: {total_L3:.3e} W")





'''
plt.plot(r_midpoints, Temp(M, Mr, r_midpoints, rin))
plt.xlim(-50, 3500)
plt.title('Temperature of accretion disk at distance from innermost stable orbit (in terms of scaled units)')
plt.xlabel('Distance from innermost stable orbit')
plt.ylabel('Temperature (K)')
plt.show()



plt.plot(Log_rMid, Temp_Logs(M, Mr, Log_rMid, Log_rin))
plt.title('Log(Temperature) of accretion disk at Log(distance) from centre')
plt.show()
'''

plt.loglog(f, f*Lpf, label='Mr=10e15, Rin=6Rg')
plt.loglog(f, f*Lpf2, label='Mr=10e15, Rin=1.2Rg')
plt.loglog(f, f*Lpf3, label='Mr=10e14, Rin=6Rg')
plt.loglog(f, f*Lpf4, label='Mr=min eddington, Rin=6Rg')
plt.loglog(f, f*Lpf5, label='Mr=max eddingont, Rin=6Rg')
plt.legend()
plt.title('Spectrum of Acctretion disk of black hole')
plt.xlabel('Log(frequency)')
plt.ylabel('Log(frequency*Luminosity)')
plt.grid(True)
plt.show()




import matplotlib.pyplot as plt
import numpy as np
import astropy 
from astropy import constants as const
from astropy import units as u
from astropy.modeling.models import BlackBody
from astropy import units as u
from astropy.visualization import quantity_support
import scipy
from scipy import integrate



MassBH = 10 * const.M_sun  #Mass of Black hole
AccR = 10**15 *u.kg/u.s  #Accretion rate
Mr2 = 10**14 *u.kg/(u.s*u.m**3)

Fstart = 14 
Fstop = 19
Fsteps = 10000
freq = np.logspace(Fstart, Fstop, Fsteps) *u.Hz #Range of frequencies


Nrings = 5000
Rg = (const.G * MassBH) / ((const.c)**2)  #Schwarzschild radius
Rin = 6 * Rg   #Innermost stable orbit
Rin2 = 1.2 * Rg   #Innermost stable orbit max spin
Rout = (10**5) * Rg   #Outermost orbit
R = np.linspace(Rin, Rout, Nrings + 1)
R2 = np.linspace(Rin2, Rout, Nrings + 1)
R_midpoints = ((R[1:] + R[:-1]) / 2) 
R_midpoints2 = ((R2[:-1] + R2[1:]) / 2) 



rin = Rin/Rg  #Scaled innermost stable orbit
rin2 = Rin2/Rg  #Scaled innermost stable orbit
rout = Rout/Rg  #Scaled innermost stable orbit

r = np.linspace(rin, rout, Nrings + 1)  
r_midpoints = (r[:-1] + r[1:]) / 2  #Array of increasingly sized disks

r2 = np.linspace(rin2, rout, Nrings + 1)  
r_midpoints2 = (r2[:-1] + r2[1:]) / 2 


Log_rin = np.log(rin)
Log_rout = np.log(rout)

Log_r = np.linspace(Log_rin, Log_rout, Nrings + 1)
Log_rMid = (Log_r[:-1] + Log_r[1:]) / 2  

Log_rin2 = np.log(rin2)
Log_r2 = np.linspace(Log_rin2, Log_rout, Nrings)
Log_rMid2 = (Log_r2[:-1] + Log_r2[1:]) / 2 


#Area of each disk
def Area(Radius):
  A = np.zeros((Nrings, Fsteps))
  A[0, :] = np.pi * Radius[0]**2
  for i in range(1, len(Radius)):
    A[i, :] = np.pi * (Radius[i]**2 - Radius[i-1]**2)
  return A*u.m**2
#print(Area(R_midpoints))


def Temp2(M, Ar, Radius, RIN):
  a = (3 * const.G* M * Ar)
  b = (8 * np.pi * const.sigma_sb.to(u.kg /u.s**3 /u.K**4) * Radius**3)
  c = ( RIN / Radius )**(1/2)
  T = ((a / b) * (1 - c))**(1/4)
  return T


#Temperatue equation
def Temp(M, Ar, Radius, RIN):
  a = (3 * const.G* M * Ar)
  b = (8 * np.pi * const.sigma_sb.to(u.kg /u.s**3 /u.K**4) * Radius**3 * Rg**3)
  c = ( RIN / Radius )**(1/2)
  T = ((a / b) * (1 - c))**(1/4)
  return T
#print(Temp(MassBH, AccR, r_midpoints, rin))


'''
#Blackbody flux equation using built in function from astropy
def flux(temp):
    bb = BlackBody(temp)
    return bb(freq)



#fluxm = [[0 for _ in range(Fsteps)] for _ in range(Nrings)]

#Blackbody flux at each temperature in Temperature equation 
for t in Temp(MassBH, AccR, R_midpoints, Rin):
   F = flux(t)
   F = F.to(u.W / (u.m**2 * u.Hz * u.sr))
   #print(F)
   #plt.semilogx(freq, F)
#plt.show()


#
for t in Temp(MassBH, AccR, Log_rMid, Log_rin):
   FL = flux(t)
   FL = FL.to(u.W / (u.m**2 * u.Hz * u.sr))
   #print(F)
   #plt.semilogx(freq, F)
#plt.show()
'''

#LogL = F * Area(Log_rMid)
#LogLsum = np.sum(LL, axis = 0)



#Blackbody flux not using built in function
def Flux(T):
   temps = T.reshape(-1, 1)
   freqs = freq.reshape(1, -1)
   a = (2 * np.pi *const.h * freqs**3)/(const.c**2)
   b = (const.h * freqs)/(const.k_B * temps)
   c = np.exp(b)
   Fv = a/(c-1)
   return Fv


#print(Flux(Temp(MassBH, AccR, r_midpoints, rin)))

#Using function
#LL = F * Area(Log_rMid)
#Lsum = np.sum(LL, axis = 0)
#LogL = F * Area(Log_rMid)
#LogLsum = np.sum(LL, axis = 0)


#using equation
L2 = (Flux(Temp2(MassBH, AccR, R_midpoints, Rin))) * Area(R_midpoints) 
L2sum = np.sum(L2, axis = 0)

#LogL2 = (Flux(Temp(MassBH, AccR, Log_rMid, Log_rin))) * Area(Log_rMid) 
#LogL2sum = np.sum(L2, axis = 0)

#Using scaled radius R/Rg and equation
Lscaled = (Flux(Temp(MassBH, AccR, r_midpoints, rin))) * Area(r_midpoints) 
LscaledSum = np.sum(Lscaled, axis = 0)


#totL = scipy.integrate.trapezoid(Lsum, freq)
TotL2 = scipy.integrate.trapezoid(L2sum, freq).to(u.W)
TotLscaled = scipy.integrate.trapezoid(LscaledSum, freq).to(u.W)

#TotLogL2 = scipy.integrate.trapezoid(LogL2sum, freq).to(u.W)


#print(f'Function L sum =: {totL}')
print(f'(Equation) L sum =: {TotL2}')
#print(f'Equation L sum =: {TotLogL2}')
print(f'(Equation+scaled r) L sum =: {TotLscaled}')


#plt.loglog(freq, Lsum, label = 'Function')
plt.ylim(10e1, 10e13)
plt.loglog(freq, L2sum, label = 'Equation')
#plt.loglog(freq, LogL2sum, label = 'Equation')
plt.loglog(freq, LscaledSum, label = 'Equation, Scaled r')
plt.legend()
plt.show()




'''
for t2 in Temp(MassBH, AccR, R_midpoints, Rin):
   F2 = Flux(t2)
   print(F2, t2)
   plt.semilogx(freq, F2, label = t2)
plt.legend()   
plt.show()


with quantity_support():
    plt.figure()
    plt.semilogx(f, f)
    plt.axvline(bb.nu_max.to(u.Hz, equivalencies=u.spectral()).value, ls='--')
    plt.show()
'''
