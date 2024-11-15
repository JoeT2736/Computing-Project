'''
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


def LumMilestone(M, Mr, Log_rin2, Log_rout, Fstart, Fstop):
    
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
total_LumMilestone = LumMilestone(M, Mr, Log_rin2, Log_rout, Fstart, Fstop)
total_L3 = L3(M, Mr2, Log_rin, Log_rout, Fstart, Fstop)
total_L4 = L(M, Mr_edd_min, Log_rin, Log_rout, Fstart, Fstop)
total_L5 = L(M, Mr_edd_max, Log_rin, Log_rout, Fstart, Fstop)

print(f"Eddington limit of Luminosity for 10 solar mass: {L_edd(M):.3e} W")
print(f"Total Luminosity for no spin, min eddington limit accretion rate: {total_L4:.3e} W")
print(f"Total Luminosity for no spin, max eddington limit accretion rate: {total_L5:.3e} W")
print(f"Total Luminosity for no spin, 10e15Mr: {total_L:.3e} W")
print(f"Total Luminosity for max spin: {total_LumMilestone:.3e} W")
print(f"Total Luminosity for no spin, 10e13Mr: {total_L3:.3e} W")






plt.plot(r_midpoints, Temp(M, Mr, r_midpoints, rin))
plt.xlim(-50, 3500)
plt.title('Temperature of accretion disk at distance from innermost stable orbit (in terms of scaled units)')
plt.xlabel('Distance from innermost stable orbit')
plt.ylabel('Temperature (K)')
plt.show()



plt.plot(Log_rMid, Temp_Logs(M, Mr, Log_rMid, Log_rin))
plt.title('Log(Temperature) of accretion disk at Log(distance) from centre')
plt.show()


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
'''

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



MassBH = 10 * const.M_sun.to_value()  #Mass of Black hole
MassS = const.M_sun.to_value()
AccR = 10**15 #*u.kg/u.s  #Accretion rate
AccR2 = 10**10 #*u.kg/u.s
Mr2 = 10**14 #*u.kg/(u.s*u.m**3)

Fstart = 13
Fstop = 19
Fsteps = 1000
freq = np.logspace(Fstart, Fstop, Fsteps) #*u.Hz #Range of frequencies


Nrings = 10000
Rg = (const.G.to_value() * MassBH) / ((const.c.to_value())**2)  #Schwarzschild radius
Rin = 6 * Rg   #Innermost stable orbit
Rin2 = 1.2 * Rg   #Innermost stable orbit max spin
Rout = (10**5) * Rg   #Outermost orbit
R = np.logspace(np.log10(Rin), np.log10(Rout), Nrings + 1) 
R2 = np.logspace(np.log10(Rin2), np.log10(Rout), Nrings + 1)
R = np.logspace(np.log10(Rin), np.log10(Rout), Nrings + 1) 
R2 = np.logspace(np.log10(Rin2), np.log10(Rout), Nrings + 1)
R_midpoints = ((R[1:] + R[:-1]) / 2) 
R_midpoints2 = ((R2[:-1] + R2[1:]) / 2) #midpoints for spinning blackhole

#Logspace for radius

#Logspace for rtadius

rin = Rin/Rg  #Scaled innermost stable orbit
rin2 = Rin2/Rg  #Scaled innermost stable orbit
rout = Rout/Rg  #Scaled innermost stable orbit

r = np.logspace(np.log10(rin), np.log10(rout), Nrings + 1)  
r_midpoints = (r[:-1] + r[1:]) / 2  #Array of increasingly sized disks

r2 = np.logspace(np.log10(rin2), np.log10(rout), Nrings + 1)  
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
  A[0, :] = 2 * np.pi * Radius[0]**2
  for i in range(1, len(Radius)):
    A[i, :] = 2 * np.pi * (Radius[i]**2 - Radius[i-1]**2)
  return A #*u.m**2
#print(Area(R_midpoints))


#for none scaled units (radius)
#for none scaled units (radius)
def Temp2(M, Ar, Radius, RIN):
  a = (3 * const.G.to_value()* M * Ar)
  b = (8 * np.pi * const.sigma_sb.to_value() * Radius**3)
  c = ( RIN / Radius )**(1/2)
  T = ((a / b) * (1 - c))**(1/4)
  return T
#print(Temp2(MassBH, AccR, R_midpoints, Rin))

#Temperatue equation for scaled units
#Temperatue equation for scaled units
def Temp(M, Ar, Radius, RIN):
  a = (3 * const.G.to_value()* M * Ar)
  b = (8 * np.pi * const.sigma_sb.to_value() * Radius**3 * Rg**3)
  c = ( RIN / Radius )**(1/2)
  T = ((a / b) * (1 - c))**(1/4)
  return T
#print(Temp(MassBH, AccR, r_midpoints, rin))

#Temp without viscous forces
def Temp3(M, Ar, Radius):
  a = (const.G.to_value()* M * Ar)
  b = (8 * np.pi * const.sigma_sb.to_value() * Radius**3)
  T = (a / b)**(1/4)
  return T



'''
plt.figure(figsize=(10, 8))
plt.plot(R_midpoints, Temp3(MassBH, AccR, R_midpoints), label = 'Non-Viscous')
plt.plot(R_midpoints, Temp2(MassBH, AccR, R_midpoints, Rin), label = 'Viscous')
plt.xlim(0, 2e6)
plt.ylabel('Temperature (K)', fontsize=20)
plt.xlabel('Distance from innermost stable orbit (m)', fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.show()
'''

'''
#Blackbody flux equation using built in function from astropy
def flux(Temp):
    bb = BlackBody(Temp*u.K)
    return bb(freq)
'''


'''
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
   a = (2 * np.pi *const.h.to_value() * freqs**3)/(const.c.to_value()**2)
   b = (const.h.to_value() * freqs)/(const.k_B.to_value() * temps)
   c = np.exp(b)
   Fv = a/(c-1)
   return Fv



Te = Temp2(MassBH, AccR, R_midpoints, Rin)
areas = Area(R_midpoints)



#print(Flux(Temp(MassBH, AccR, r_midpoints, rin)))

#Using function
#LL = F * Area(Log_rMid)
#Lsum = np.sum(LL, axis = 0)
#LogL = F * Area(Log_rMid)
#LogLsum = np.sum(LL, axis = 0)




#Milestone values
LumMilestone = (Flux(Temp2(MassBH, AccR, R_midpoints, Rin))) * Area(R_midpoints)  
LumMilestonesum = np.sum(LumMilestone, axis = 0)


#milestone values, max spin black hole
L3 = (Flux(Temp2(MassBH, AccR, R_midpoints2, Rin2))) * Area(R_midpoints2) 
L3sum = np.sum(L3, axis = 0)

#Object mass = mass of sun
LumSunMass = (Flux(Temp2(MassS, AccR, R_midpoints, Rin))) * Area(R_midpoints)
LumSunMassSum = np.sum(LumSunMass, axis = 0)

#mass of sun with spin
#mass of sun with spin
LumSunMassSpin = (Flux(Temp2(MassS, AccR, R_midpoints2, Rin2))) * Area(R_midpoints2)
LumSunMassSumSpin = np.sum(LumSunMassSpin, axis = 0)

#Different accretion rate
#Different accretion rate
LumdAccr = (Flux(Temp2(MassBH, AccR2, R_midpoints, Rin))) * Area(R_midpoints)
LumdAccrSum = np.sum(LumdAccr, axis = 0)

#D Accr with spin
#D Accr with spin
LumdAccrSpin = (Flux(Temp2(MassBH, AccR2, R_midpoints2, Rin2))) * Area(R_midpoints2)
LumdAccrSumSpin = np.sum(LumdAccrSpin, axis = 0)






#Using scaled radius R/Rg and equation
Lscaled = (Flux(Temp(MassBH, AccR, r_midpoints, rin))) * Area(r_midpoints)  
LscaledSum = np.sum(Lscaled, axis = 0)



TotLumMilestone = scipy.integrate.trapezoid(LumMilestonesum, freq)*u.W  #.to(u.W)
TotLscaled = scipy.integrate.trapezoid(LscaledSum, freq)*u.W  #.to(u.W)
TotL3 = scipy.integrate.trapezoid(L3sum, freq)*u.W  #.to(u.W)
TotLumSunMass = scipy.integrate.trapezoid(LumSunMassSum, freq)*u.W  #.to(u.W)
TotLumSunMassSpin = scipy.integrate.trapezoid(LumSunMassSumSpin, freq)*u.W  #.to(u.W)
TotLumdAccr = scipy.integrate.trapezoid(LumdAccrSum, freq)*u.W  #.to(u.W)
TotLumdAccrSpin = scipy.integrate.trapezoid(LumdAccrSumSpin, freq)*u.W  #.to(u.W)




'''
print(f'(Milestone) L sum = {TotLumMilestone}')
print(f'(Equation+scaled r) L sum = {TotLscaled}')
print(f'(Milestone + spin) L sum = {TotL3}')
print(f'(Sun mass) L sum = {TotLumSunMass}')
print(f'(Sun mass + spin) L sum = {TotLumSunMassSpin}')
print(f'(D Accr) L sum = {TotLumdAccr}')
print(f'(D Accr + spin) L sum = {TotLumdAccrSpin}')
'''
'''
spectra = []
for temp in Temp2(MassBH, AccR, R_midpoints, Rin):
  area = Area(R_midpoints)
  spectra.append(Flux(temp) * area)

print(spectra)
'''


#spectra = Flux(5000)

'''
fig, ax = plt.subplots()
plt.figure(figsize=(10,6))
plt.ylim(10e20, 10e31)
#for i, temp in enumerate(Temp2(MassBH, AccR, R_midpoints, Rin)):
  #plt.loglog(freq, freq * spectra[i], label=f'{temp} K')

MileStone, = plt.loglog(freq, freq * LumMilestonesum, linestyle='-', color = 'blue')
MaxSpin, = plt.loglog(freq, freq * L3sum, linestyle='--', color = 'blue')
#SunMass, = plt.loglog(freq, freq * LumSunMassSum, linestyle='-', color = 'green')
#SunMassSpin, = plt.loglog(freq, freq * LumSunMassSumSpin, linestyle='--', color = 'green')
dAccr, = plt.loglog(freq, freq * LumdAccrSum, linestyle='-', color = 'red')
dAccrSpin, = plt.loglog(freq, freq * LumdAccrSumSpin, linestyle='--', color = 'red')

#Legend (linestyle and colour)
#Legend (linestyle and colour)
no_spin_line = plt.Line2D([0], [0], color='black', linestyle='-', label='No Spin')
max_spin_line = plt.Line2D([0], [0], color='black', linestyle='--', label='Max Spin')
mass_10msun = plt.Line2D([0], [0], color='blue', linestyle='-', label=r'AccR=10$^{15}$Kg/s')
mass_1msun = plt.Line2D([0], [0], color='red', linestyle='-', label=r'AccR=10$^{10}$Kg/s')
plt.legend(handles=[no_spin_line, max_spin_line, mass_10msun, mass_1msun])

#plt.loglog(freq, LogLumMilestonesum, label = 'Equation')
#plt.loglog(freq, LscaledSum, label = 'Equation, Scaled r')
plt.ylabel('Luminosity (W)', fontsize=16)
plt.xlabel('Frequency (Hz)', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()
'''

#plot each blackbodies on same plot
#make spectrums to location in ring, draw blackhole?



fig, ax = plt.subplots()
plt.figure(figsize=(10,6))
plt.ylim(20, 32)
#for i, temp in enumerate(Temp2(MassBH, AccR, R_midpoints, Rin)):
  #plt.loglog(freq, freq * spectra[i], label=f'{temp} K')

MileStone, = plt.plot(np.log10(freq), np.log10(freq * LumMilestonesum), linestyle='-', color = 'blue')
MaxSpin, = plt.plot(np.log10(freq), np.log10(freq * L3sum), linestyle='--', color = 'blue')
#SunMass, = plt.loglog(freq, freq * LumSunMassSum, linestyle='-', color = 'green')
#SunMassSpin, = plt.loglog(freq, freq * LumSunMassSumSpin, linestyle='--', color = 'green')
dAccr, = plt.plot(np.log10(freq), np.log10(freq * LumdAccrSum), linestyle='-', color = 'red')
dAccrSpin, = plt.plot(np.log10(freq), np.log10(freq * LumdAccrSumSpin), linestyle='--', color = 'red')

#Legend (linestyle and colour)
#Legend (linestyle and colour)
no_spin_line = plt.Line2D([0], [0], color='black', linestyle='-', label='No Spin')
max_spin_line = plt.Line2D([0], [0], color='black', linestyle='--', label='Max Spin')
mass_10msun = plt.Line2D([0], [0], color='blue', linestyle='-', label=r'AccR=10$^{15}$Kg/s')
mass_1msun = plt.Line2D([0], [0], color='red', linestyle='-', label=r'AccR=10$^{10}$Kg/s')
plt.legend(handles=[no_spin_line, max_spin_line, mass_10msun, mass_1msun], fontsize=12)

#plt.loglog(freq, LogLumMilestonesum, label = 'Equation')
#plt.loglog(freq, LscaledSum, label = 'Equation, Scaled r')
plt.ylabel(r'$\log_{10}(fL_f)$ W', fontsize=16)
plt.xlabel(r'$\log_{10}(f)$ Hz', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
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


'''
plt.figure(figsize=(10,6))
plt.ylim(20, 32)
index = [100, 2000, 3500, 5000, 6500, 8000, 9500]
colors = ['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet']
for j, color in zip(index, colors):
  #print(Flux(Te[i]))
  L = Flux(Te[j]) * areas[j]
  Lsum = np.sum(L, axis = 0)
  temp_label = f'{Te[j]:.2g}'
  temp_label = temp_label.replace('e+0',r'\times 10^{')+'}'
  plt.plot(np.log10(freq), np.log10(freq * Lsum), linestyle='--', color = color, label = f'${temp_label}$ K')
plt.plot(np.log10(freq), np.log10(freq * LumMilestonesum), linestyle='-', color = 'black')
plt.legend(fontsize=12)
plt.ylabel(r'$\log_{10}(fL_f)$ W', fontsize=16)
plt.xlabel(r'$\log_{10}(f)$ Hz', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()
'''
