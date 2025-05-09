#how far away can you see balck hole as function of accretion rate?????
#set black hole at a distance then look at received flux ^^^
#flux vs distance for observable BHs
#extra non thermal process in black hole not in this modeled >>>> Much higher energies than possible for thermal
#what AGNs can be seen????
#filters???? V-band??? B-band???
#size of black holes observables in a galaxy?
#Ledd vs mass (add in observed black holes -> do they fit this prediction???)
#Ledd -> M_dot_critical
#is Ledd actual limit??
#variation in accretion rates??
#this code is for static disc -> is it good for moving disc??? accretion rate as function of time???
#radiation pressure >>>> hot flow
#winds
#optically thick to thin disc (becomes thin close to BH)
#Rakshit etal    SDSS

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
from astropy.cosmology import WMAP9 as cosmo
from astropy.io.votable import parse
import pandas as pd
#from astroquery.simbad import Simbad
import scipy.integrate
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator
from scipy.stats import gaussian_kde
from matplotlib import rc
import sys

#def votable_to_pandas(votable_file):
 #   votable = parse(votable_file)
  ##  table = votable.get_first_table().to_table(use_names_over_ids=True)
    #return table.to_pandas()


#result_table = Simbad.query_object("M [1-9]", wildcard=True) 




'''
# Extracted data
data = """
0.1611	3.266E-17	0.00000	0.00000
0.1663	1.150E-18	0.00000	0.00000
1.865	1.544E-17	0.00000	0.00000
1.925	5.382E-17	0.00000	0.00000
1.987	1.264E-17	0.00000	0.00000
2.051	1.088E-17	0.00000	0.00000
2.051	2.002E-19	0.00000	0.00000
2.257	5.658E-17	0.00000	0.00000
3.202	2.675E-17	0.00000	0.00000
3.306	1.982E-17	0.00000	0.00000
6.245	4.633E-17	0.00000	0.00000
7.322	4.633E-17	0.00000	0.00000
7.803	8.870E-17	0.00000	0.00000
13.83	5.658E-17	0.00000	0.00000
15.21	9.324E-17	0.00000	0.00000
15.70	7.635E-17	0.00000	0.00000
176.1	1.015E-23	0.00000	0.00000
1631	2.246E-12	0.00000	0.00000
3501	1.830E-11	0.00000	0.00000
3851	1.110E-11	0.00000	0.00000
5464	6.379E-11	0.00000	0.00000
5641	2.897E-13	0.00000	0.00000
5641	1.296E-12	0.00000	0.00000
5641	5.491E-11	0.00000	0.00000
5823	4.497E-11	0.00000	2961000
9089	4.497E-11	0.00000	0.00000
9383	4.774E-13	0.00000	0.00000
9686	6.411E-12	0.00000	0.00000
9686	9.094E-12	0.00000	0.00000
14180	3.015E-11	0.00000	0.00000
14180	3.503E-11	0.00000	0.00000
16110	1.740E-11	0.00000	0.00000
16630	3.186E-12	0.00000	0.00000
17170	4.751E-12	0.00000	0.00000
44580	1.432E-12	0.00000	0.00000
44580	1.173E-12	0.00000	0.00000
47510	6.411E-12	0.00000	0.00000
50630	7.484E-13	0.00000	0.00000
52260	8.230E-12	0.00000	0.00000
52260	1.749E-12	0.00000	0.00000
52260	1.607E-15	0.00000	0.00000
63250	3.046E-13	0.00000	0.00000
98730	3.366E-13	0.00000	0.00000
98730	2.742E-12	0.00000	0.00000
98730	1.749E-12	0.00000	0.00000
101900	1.656E-11	0.00000	0.00000
101900	3.701E-12	0.00000	0.00000
101900	7.484E-13	0.00000	0.00000
105200	1.227E-11	0.00000	0.00000
149200	2.234E-11	0.00000	0.00000
149200	3.382E-14	0.00000	0.00000
154100	2.729E-11	0.00000	0.00000
218600	2.125E-11	0.00000	0.00000
233000	3.015E-11	0.00000	0.00000
291100	1.740E-11	0.00000	0.00000
291100	2.596E-11	0.00000	0.00000
291100	1.749E-12	0.00000	0.00000
291100	2.883E-12	0.00000	0.00000
300500	3.871E-11	0.00000	0.00000
375400	9.512E-11	0.00000	0.00000
400100	4.497E-11	0.00000	0.00000
413000	9.999E-11	0.00000	0.00000
413000	1.227E-11	0.00000	0.00000
413000	5.250E-12	0.00000	0.00000
644700	1.820E-10	0.00000	0.00000
644700	2.457E-10	0.00000	0.00000
665500	2.223E-10	0.00000	0.00000
858400	3.664E-10	0.00000	0.00000
886100	4.944E-10	0.00000	0.00000
886100	4.049E-10	0.00000	0.00000
944300	4.751E-12	0.00000	0.00000
944300	2.854E-10	0.00000	0.00000
1217000	4.049E-10	0.00000	0.00000
1217000	5.464E-10	0.00000	0.00000
1257000	3.485E-10	0.00000	0.00000
1339000	4.994E-12	0.00000	0.00000
1570000	3.701E-12	0.00000	0.00000
1674000	5.744E-10	0.00000	0.00000
1674000	8.271E-13	0.00000	0.00000
1783000	1.283E-10	0.00000	0.00000
1783000	3.015E-11	0.00000	0.00000
1841000	8.271E-13	0.00000	0.00000
1962000	2.742E-12	0.00000	0.00000
2025000	4.703E-10	0.00000	0.00000
2091000	4.408E-17	0.00000	0.00000
2091000	2.469E-11	0.00000	0.00000
2159000	2.742E-12	0.00000	0.00000
2159000	1.296E-12	0.00000	0.00000
2228000	7.868E-13	0.00000	0.00000
2531000	1.923E-11	0.00000	0.00000
2612000	1.664E-12	0.00000	0.00000
2697000	7.120E-13	0.00000	0.00000
2697000	3.851E-10	0.00000	0.00000
2874000	9.608E-13	0.00000	0.00000
2874000	6.706E-11	0.00000	0.00000
2967000	1.349E-10	0.00000	0.00000
2967000	3.871E-11	0.00000	0.00000
2967000	2.246E-12	0.00000	0.00000
2967000	1.749E-12	0.00000	0.00000
3063000	3.520E-12	0.00000	0.00000
3162000	8.608E-11	0.00000	0.00000
3478000	2.854E-10	0.00000	0.00000
3591000	1.167E-11	0.00000	0.00000
3707000	4.069E-11	0.00000	0.00000
3707000	8.694E-13	0.00000	0.00000
3827000	1.943E-13	0.00000	0.00000
3827000	2.223E-10	0.00000	0.00000
3950000	4.299E-12	0.00000	0.00000
3950000	1.061E-12	0.00000	0.00000
3950000	2.635E-14	0.00000	0.00000
3950000	1.432E-12	0.00000	0.00000
4078000	1.179E-13	0.00000	0.00000
4078000	5.574E-14	0.00000	0.00000
4078000	4.320E-13	0.00000	0.00000
4781000	1.647E-10	0.00000	0.00000
4781000	8.694E-13	0.00000	0.00000
4781000	1.296E-12	0.00000	0.00000
4935000	2.481E-12	0.00000	0.00000
5095000	9.560E-12	0.00000	0.00000
5260000	4.110E-13	0.00000	0.00000
5260000	5.830E-13	0.00000	0.00000
6571000	4.320E-13	0.00000	0.00000
6571000	2.897E-13	0.00000	0.00000
6784000	4.110E-13	0.00000	0.00000
6784000	1.362E-12	0.00000	0.00000
6784000	2.136E-12	0.00000	0.00000
6784000	3.046E-13	0.00000	0.00000
7704000	2.042E-13	0.00000	0.00000
7704000	1.758E-13	0.00000	0.00000
7953000	4.969E-11	0.00000	0.00000
7953000	3.170E-11	0.00000	0.00000
7953000	5.859E-14	0.00000	0.00000
8210000	4.277E-11	0.00000	0.00000
8475000	4.798E-14	0.00000	0.00000
8749000	3.538E-13	0.00000	0.00000
8749000	1.943E-13	0.00000	0.00000
9032000	3.719E-13	0.00000	0.00000
9625000	1.672E-13	0.00000	0.00000
9625000	1.179E-13	0.00000	0.00000
9625000	8.738E-14	0.00000	0.00000
1.202E7	3.170E-11	0.00000	0.00000




"""

# Split the data into lines
lines = data.strip().split('\n')

# Initialize lists
L4 = []
M4 = []
T4 = []

# Extract values into lists
for line in lines:
    values = line.split()
    L4.append(float(values[0]))
    M4.append(float(values[1]))
    T4.append(values[2])

# Print the lists to verify
print("L4:", L4)
print("M4:", M4)
#print("T4:", T4)
'''




MassSMBH1 =  10**4 * const.M_sun.to_value()
MassSMBH2 =  10**8 * const.M_sun.to_value()
MassSMBH3 =  10**11 * const.M_sun.to_value()
MassBH =  10 * const.M_sun.to_value()  #Mass of Black hole
Mass2 = 100 * const.M_sun.to_value()
AccR = 10**15 #*u.kg/u.s  #Accretion rate
AccR2 = 10**16 #*u.kg/u.s
Mr2 = 10**14 #*u.kg/(u.s*u.m**3)

#print((AccR/const.M_sun.to_value())*365.25*24*60*60)

Fstart = 0.01
Fstop = 30
Fsteps = 1000
freq = np.logspace(Fstart, Fstop, Fsteps) #* u.Hz #Range of frequencies


Nrings = 10000
#Rg = (const.G.to_value() * MassBH) / ((const.c.to_value())**2)  #Schwarzschild radius

def Rg(M):
  return (const.G.to_value() * M) / ((const.c.to_value())**2)

'''
Rin = 6 * Rg()   #Innermost stable orbit
Rin2 = 1.2 * Rg()   #Innermost stable orbit max spin
Rout = (10**5) * Rg()   #Outermost orbit
R = np.logspace(np.log10(Rin), np.log10(Rout), Nrings + 1) 
R2 = np.logspace(np.log10(Rin2), np.log10(Rout), Nrings + 1)
R = np.logspace(np.log10(Rin), np.log10(Rout), Nrings + 1) 
R2 = np.logspace(np.log10(Rin2), np.log10(Rout), Nrings + 1)
R_midpoints = 10**((np.log10(R[1:]) + np.log10(R[:-1])) / 2) 
R_midpoints2 = 10**((np.log10(R2[:-1]) + np.log10(R2[1:])) / 2) #midpoints for spinning blackhole
'''


def Rin(M):
  return 6 * Rg(M)   #Innermost stable orbit

def Rin2(M):
  return 1.2 * Rg(M)   #Innermost stable orbit max spin

def Rout(M):
  return (10**5) * Rg(M)   #Outermost orbit



def R(M):
  return np.logspace(np.log10(Rin(M)), np.log10(Rout(M)), Nrings + 1)

def R_midpoints(M):
  return 10**((np.log10(R(M)[1:]) + np.log10(R(M)[:-1])) / 2)


'''
print(R_midpoints(MassBH))
print(R_midpoints(Mass2))
'''


#scaled radius 
def r(M):
  return R_midpoints(M)/Rg(M)

def rin(M):
  return Rin(M)/Rg(M)






#Logspace for radius

#Logspace for rtadius
'''
rin = Rin/Rg()  #Scaled innermost stable orbit
rin2 = Rin2/Rg()  #Scaled innermost stable orbit
rout = Rout/Rg() #Scaled innermost stable orbit

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
'''

#Area of each disk
def Area(Radius):
  A = np.zeros((Nrings, Fsteps))
  A[0, :] = 2 * np.pi * Radius[0]**2
  for i in range(1, len(Radius)):
    A[i, :] = 2 * np.pi * (Radius[i]**2 - Radius[i-1]**2)
  return A #*u.m**2



#for none scaled units (radius)
def Temp2(M, Ar, Radius, Rin):
    a = (3 * const.G.to_value() * M * Ar)
    b = (8 * np.pi * const.sigma_sb.to_value() * Radius**3)
    c = (Rin / Radius)**0.5

    T = ((a / b) * (1 - c))**0.25
    return T



'''
#Temperatue equation for scaled units
def Temp(M, Ar, Radius, RIN):
  a = (3 * const.G.to_value()* M * Ar)
  b = (8 * np.pi * const.sigma_sb.to_value() * Radius**3 * Rg(M)**3)
  c = ( RIN / Radius )**(1/2)
  T = ((a / b) * (1 - c[:, np.newaxis]))**(1/4)
  return T
#print(Temp(MassBH, AccR, r(MassBH), rin(MassBH)))
'''

#Schwarzschild radius
def R_s(M):
  return (2 * const.G.to_value() * M) / const.c.to_value()**2

'''
#Temp from Frank et al. 1992       i think this is for max spinning black hole
def Temp3(M, Ar, Radius, Rs):
  a = (3 * const.G.to_value()* M * Ar)
  b = (8 * np.pi * const.sigma_sb.to_value() * Radius**3)
  c = ( Rs / Radius )**(1/2)
  T = ((a / b) * (1 - c))**(1/4)
  return T
'''



#T2 = Temp2(MassBH, AccR, R_midpoints(MassBH), Rin(MassBH))
#T3 = Temp3(MassBH, AccR, R_midpoints)

#print(T2.shape)





Akkr = 1e-8*const.M_sun.to_value()/(365.25*24*60*60)
Akkr2 = 1e-7*const.M_sun.to_value()/(365.25*24*60*60)

#GOOD TEMPERATURE PLOT MIGHT NEED THIS

#fig, ax1 = plt.subplots()
plt.figure(figsize=(8, 6))

#fig.patch.set_facecolor('#D5D5D5') 
#ax1.set_facecolor('#D5D5D5') 

#print((Rin(MassBH) + R_midpoints(MassBH))/1000)

plt.plot((Rin(MassBH) + R_midpoints(MassBH))/1000, Temp2(MassBH, Akkr, R_midpoints(MassBH), Rin(MassBH))/1e6, color = 'black', 
         label = 'Mass=10$M_{\odot}$, AccR=10$^{-8}$$M_{\odot}$Yr$^{-1}$', linewidth = 2)
plt.plot((Rin(MassBH) + R_midpoints(MassBH))/1000, Temp2(MassBH, Akkr2, R_midpoints(MassBH), Rin(MassBH))/1e6, color = 'blue', linestyle = ':',
         label = 'Mass=10$M_{\odot}$, AccR=10$^{-8}$$M_{\odot}$Yr$^{-1}$', linewidth = 2)
plt.plot((Rin(Mass2) + R_midpoints(Mass2))/1000, Temp2(Mass2, Akkr, R_midpoints(Mass2), Rin(Mass2))/1e6, color = 'red', linestyle = '--',
         label = 'Mass=100$M_{\odot}$, AccR=10$^{-7}$$M_{\odot}$Yr$^{-1}$', linewidth = 2)
plt.xlim(-0.003e4, 0.007e6)
plt.ylabel('Temperature (MK)', fontsize=18)
plt.xlabel('Distance from centre of Black hole (km)', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=16, length=7, width=2, top=True, right=True) 
plt.tick_params(axis='both', which='minor', labelsize=10, length=4, width=1.5, top=True, right=True) 
plt.legend(fontsize=16, frameon=True, loc='best')
plt.minorticks_on()
plt.tight_layout()
plt.savefig('Temperature', dpi=300, bbox_inches='tight')







#Blackbody flux not using built in function
def Flux(T):
   temps = T.reshape(-1, 1)
   tem = np.maximum(temps, 1e-10)
   freqs = freq.reshape(1, -1)
   a = (2 * np.pi *const.h.to_value() * freqs**3)/(const.c.to_value()**2)
   b = (const.h.to_value() * freqs)/(const.k_B.to_value() * tem) 
   #b = b.astype(np.float32)
   #b = np.clip(b, -100, 100)
   c = np.exp(b)
   Fv = a/(c-1)
   return Fv

'''
def Flux(T):
    temps = np.maximum(T, 1e-10)  # Prevent division by zero
    a = (2 * np.pi * const.h.to_value() * freq**3) / (const.c.to_value()**2)

    chunk_size = 100  # Process in batches to avoid memory errors
    Fv = np.zeros((T.shape[0], len(freq)))  # Preallocate output array
    
    for i in range(0, T.shape[0], chunk_size):
        temp_chunk = temps[i:i+chunk_size]  # Get a batch
        temp_chunk = temp_chunk.reshape(-1, 1)  # Ensure correct shape for broadcasting

        b = (const.h.to_value() * freq) / (const.k_B.to_value() * temp_chunk)
        #b = np.clip(b, -100, 100)  # Prevent overflow
        c = np.exp(b)
        
        Fv[i:i+chunk_size, :] = a / (c - 1)

    return Fv
'''


  
#print(Flux(Temp2(MassBH, AccR, R_midpoints(MassBH), Rin(MassBH))))

'''
#so i can change the frequency bands
def Flux3(T, freq):
   temps = T.reshape(-1, 1)
   freqs = freq.reshape(1, -1)
   a = (2 * np.pi *const.h.to_value() * freqs**3)/(const.c.to_value()**2)
   b = (const.h.to_value() * freqs)/(const.k_B.to_value() * temps)
   c = np.exp(b)
   Fv = a/(c-1)
   return Fv
'''


#Luminosity eddington limit
def L_edd(Mass):
  return (4 * np.pi * const.G.value * Mass * const.m_p.value * const.c.value)/const.sigma_T.value

eta = 0.1  #standard efficiency with which mass is assumed to produce energy near the event horizon

#max accretion rate
def AccR_Edd(Mass):
  return L_edd(Mass)/(eta * const.c.value**2)

#print(AccR_Edd(70e6*const.M_sun.to_value()))


#print(Flux(Temp(MassBH, AccR, r_midpoints, rin)))

#Using function
#LL = F * Area(Log_rMid)
#Lsum = np.sum(LL, axis = 0)
#LogL = F * Area(Log_rMid)
#LogLsum = np.sum(LL, axis = 0)

'''
def Lsum3(M, Ar, Radius, RIN, freq):
  F = Flux3(Temp2(M, Ar, Radius, RIN), freq)
  L = F * Area(Radius)
  Lsum = np.sum(L, axis = 0)
  return Lsum
'''



def Lsum(M, Ar, Radius, RIN):
  F = Flux(Temp2(M, Ar, Radius, RIN))
  L = F * Area(Radius)
  Lsum = np.sum(L, axis = 0)
  return Lsum

'''
#L from temp3
def Lsum2(M, Ar, Radius, Rs):
  F = Flux(Temp3(M, Ar, Radius, Rs))
  L = F * Area(Radius)
  Lsum = np.sum(L, axis = 0)
  return Lsum


def Lsum4(M, Ar, Radius, RIN):
  F = Flux(Temp(M, Ar, Radius, RIN))
  L = F * Area(Radius) * 4 * np.pi * Rg(M)**2  #idk how the scaled units work 
  Lsum = np.sum(L, axis = 0)
  return Lsum
'''

#print(Lsum4(MassBH, AccR, r(MassBH), rin(MassBH)))



#Flux to luminosity using redshift
#def Luminosity(F, z):
#  return F * 4 * np.pi * cosmo.luminosity_distance(z).to(u.m).value**2


'''
#lum in erg s-1 Hz-1
def Lv(Rs, i, M, Accr, freq):
  a = (3 * const.G.to_value() * M * Accr)
  b = 8 * np.pi * const.sigma_sb.to_value() * Rs**3
  c = (a / b)**(1/4)
  d = c**(8/3)
  e = Rs**2
  return 2.4e-18 * e * d * np.cos(i) * freq**(1/3)


#Shukara Sunyaev 1976
def AcRate(L, M):
  return (L/(0.25*(2 * Rg(M)/Rin(M))))**(1/2)

#print(AcRate(Luminosity(5.19e-11, -0.00015677512474313146), 70e6*const.M_sun.to_value()))

#Bolometric corrections
def Lum_BolCorr(Lband, K):
  return Lband * K
'''


#flux received 
#def Flux2(L, z):
#  return L / (4 * np.pi * cosmo.luminosity_distance(z).to(u.m).value**2)

'''
F1M81 = Flux2(Lsum(70e6*const.M_sun.to_value(), 0.00001*AccR_Edd(70e6*const.M_sun.to_value()), 
R_midpoints(70e6*const.M_sun.to_value()), Rin(70e6*const.M_sun.to_value())), -0.00015677512474313146)


#Not very accurate, not sure why
#accretion rate knowing mass and luminoisty
#Mass in solar mass    Luminosity in 10^44 erg/s      Wavelegth in Angstrom      Accretion in solar mass per year
def M_dot(Mass, Luminosity, Wavelength):
  return np.exp(1.5*np.log(Wavelength*Luminosity) + 2*np.log(Wavelength) + 0.213 - np.log(Mass))

#print(M_dot(7e7, 7.89e-1*Luminosity(0.0001*6e-8, -0.00015677512474313146)/10e44, 100000)*const.M_sun.to_value()/(365.25*24*60*60))

#bondi accretion rate
def Bondi(Mass, Density, SoundSpeed):
  return (4 * np.pi * const.G.to_value()**2 * Mass**2 * Density) / SoundSpeed**3


F2M81 = Flux2(Lsum(70e6*const.M_sun.to_value(), M_dot(7e7, 1*Luminosity(0.0001*6e-8, -0.00015677512474313146)/10e44, 10000)*const.M_sun.to_value()/(365.25*24*60*60), 
R_midpoints(70e6*const.M_sun.to_value()), Rin(70e6*const.M_sun.to_value())), -0.00015677512474313146)


F3M81 = Flux2(Lsum2(70e6*const.M_sun.to_value(), M_dot(7e7, 1*Luminosity(0.0001*6e-8, -0.00015677512474313146)/10e44, 10000)*const.M_sun.to_value()/(365.25*24*60*60),
R_midpoints(70e6*const.M_sun.to_value()), R_s(70e6*const.M_sun.to_value())), -0.00015677512474313146)


F4M81 = Flux2(Lsum2(70e6*const.M_sun.to_value(),0.00001*AccR_Edd(70e6*const.M_sun.to_value()),
R_midpoints(70e6*const.M_sun.to_value()), R_s(70e6*const.M_sun.to_value())), -0.00015677512474313146)


F5M81 = Flux2(Lsum(70e6*const.M_sun.to_value(), AcRate(Luminosity(5.19e-11, -0.00015677512474313146), 70e6*const.M_sun.to_value()),
R_midpoints(70e6*const.M_sun.to_value()), Rin(70e6*const.M_sun.to_value())), -0.00015677512474313146)

F6M81 = Flux2(Lsum2(70e6*const.M_sun.to_value(), AcRate(Luminosity(5.19e-11, -0.00015677512474313146), 70e6*const.M_sun.to_value()),
R_midpoints(70e6*const.M_sun.to_value()), R_s(70e6*const.M_sun.to_value())), -0.00015677512474313146)



VBand = const.c.to_value()/0.55e-6
#np.logspace(np.log10(const.c.to_value()/0.44e-6), 
#np.log10(const.c.to_value()/0.44e-6), 1000)

F1M81_VBand = np.sum(Flux2(Lsum3(70e6*const.M_sun.to_value(),  AcRate(Luminosity(5.19e-11, -0.00015677512474313146), 70e6*const.M_sun.to_value()),
R_midpoints(70e6*const.M_sun.to_value()), Rin(70e6*const.M_sun.to_value()), VBand), -0.00015677512474313146), axis = 0)

#conversion to magnitude
def mag(f):
  return -2.5*np.log10(f) - 48.6

#print(mag(F1M81_VBand))

em_spectrum = {
    'Radio': (1e6 / 1e9, 3e9 / 1e9),
    'Microwave': (3e9 / 1e9, 3e11 / 1e9),
    'Infrared': (3e11 / 1e9, 4e14 / 1e9),
    'Optical': (4e14 / 1e9, 7.5e14 / 1e9),
    'Ultraviolet': (7.5e14 / 1e9, 3e16 / 1e9),
    'X-ray': (3e16 / 1e9, 3e19 / 1e9),
    'Gamma-ray': (3e19 / 1e9, 3e24 / 1e9)
}


#refeeence flux 
F0U = 1.81e-20 * 1000
F0B = 4.26e-20 * 1000
F0V = 3.64e-20 * 1000
F0R = 3.08e-20 * 1000
F0I = 2.55e-20 * 1000
F0J = 1.60e-20 * 1000
F0H = 1.08e-20 * 1000
F0K = 0.67e-20 * 1000

F0s = [F0U, F0B, F0V, F0R, F0I, F0J, F0H, F0K]
'''
'''
mags = {13.551, 12.18, 11.48, 12.284, 11.613, 11.081, 10.448, 8.5, 7.788, 7.381}

#Magnitude to flux
def fluxx(mag, F0):
  return 10**(mag/(-2.5))*F0

print(fluxx(13.551, F0U))

#Total flux received from all magnitude bands
def Tfluxx(mag, F0):
  return np.sum(fluxx(mag, F0))




NGC4151 = Flux2(Lsum(20e6*const.M_sun.to_value(), AcRate(Luminosity(Tfluxx(mags, F0s), 0.00315215), 20e6*const.M_sun.to_value()),
R_midpoints(20e6*const.M_sun.to_value()), Rin(20e6*const.M_sun.to_value())), 0.00315215)


#Plot of a Seyfert 1 and Seyfert 2
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

ax1.set_xlim(-0.3, 7)
ax1.set_ylim(-21, -10)

ax1.plot(np.log10((freq/10e9)), np.log10(freq * NGC4151), color='blue', linestyle='-', label='M81')

colors = {
    'Radio': 'lightblue',
    'Microwave': 'lightgreen',
    'Infrared': 'lightcoral',
    'Optical': 'lightyellow',
    'Ultraviolet': 'lightpink',
    'X-ray': 'lightgray',
    'Gamma-ray': 'lightcyan'
}

for band, (start, end) in em_spectrum.items():
    ax1.axvspan(np.log10(start), np.log10(end), color=colors[band], alpha=0.3)


ax1.set_ylabel(r'$\log_{10}(vF_v)$ W m$^{-2}$', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(v)$ GHz', fontsize=16)

ax1.legend(fontsize=12, loc='upper left')
plt.show()
'''


'''
Lsun = 3.828e33  #Luminosity of sun in erg s-1

#luminosity given b-band magnitude
def Luminosity2(mag, z):
  return 4.26e-20*(10**(-mag/2.5)) * 4 * np.pi * cosmo.luminosity_distance(z).to(u.m).value**2

#print(Luminosity2(7.89, -0.00015677512474313146))

#L in terms of B-band luminosity(erg s-1)/10e9*L_sun
#M in terms of mass/10e8*M_sun
def Accretion(M, L):
  return ((L/(13.8 * M**(2/3)))**(3/2))/1000

#print(Accretion(0.7, Luminosity2(7.89, -0.00015677512474313146)))


#Plot of radio-loud 3C 273 (think this is a blazar) 
F3C273 = Flux2(Lsum(886e6*const.M_sun.to_value(), AcRate(Luminosity(2.54e-13, 0.15756751), 886e6*const.M_sun.to_value()),
R_midpoints(886e6*const.M_sun.to_value()), Rin(886e6*const.M_sun.to_value())), 0.15756751)


#Plot of Seyfert 1 NGC 6814
FNGC6814 = Flux2(Lsum(3.9e6*const.M_sun.to_value(), AcRate(Luminosity(1.23e-12, 0.00579177), 3.9e6*const.M_sun.to_value()),
R_midpoints(3.9e6*const.M_sun.to_value()), Rin(3.9e6*const.M_sun.to_value())), 0.00579177)
'''









M81x = np.array([0.1611, 0.1663, 1.865, 1.925, 1.987, 2.051, 2.051, 2.257, 3.202, 3.306, 6.245, 7.322, 7.803, 13.83, 15.21, 15.7, 176.1, 1631.0, 3501.0, 3851.0, 5464.0, 5641.0, 5641.0, 5641.0, 5823.0, 9089.0, 9383.0, 9686.0, 9686.0, 14180.0, 14180.0, 16110.0, 16630.0, 17170.0, 44580.0, 44580.0, 47510.0, 50630.0, 52260.0, 52260.0, 52260.0, 63250.0, 98730.0, 98730.0, 98730.0, 101900.0, 101900.0, 101900.0, 105200.0, 149200.0, 149200.0, 154100.0, 218600.0, 233000.0, 291100.0, 291100.0, 291100.0, 291100.0, 300500.0, 375400.0, 400100.0, 413000.0, 413000.0, 413000.0, 644700.0, 644700.0, 665500.0, 858400.0, 886100.0, 886100.0, 944300.0, 944300.0, 1217000.0, 1217000.0, 1257000.0, 1339000.0, 1570000.0, 1674000.0, 1674000.0, 1783000.0, 1783000.0, 1841000.0, 1962000.0, 2025000.0, 2091000.0, 2091000.0, 2159000.0, 2159000.0, 2228000.0, 2531000.0, 2612000.0, 2697000.0, 2697000.0, 2874000.0, 2874000.0, 2967000.0, 2967000.0, 2967000.0, 2967000.0, 3063000.0, 3162000.0, 3478000.0, 3591000.0, 3707000.0, 3707000.0, 3827000.0, 3827000.0, 3950000.0, 3950000.0, 3950000.0, 3950000.0, 4078000.0, 4078000.0, 4078000.0, 4781000.0, 4781000.0, 4781000.0, 4935000.0, 5095000.0, 5260000.0, 5260000.0, 6571000.0, 6571000.0, 6784000.0, 6784000.0, 6784000.0, 6784000.0, 7704000.0, 7704000.0, 7953000.0, 7953000.0, 7953000.0, 8210000.0, 8475000.0, 8749000.0, 8749000.0, 9032000.0, 9625000.0, 9625000.0, 9625000.0, 12020000.0])
M81y = np.array([3.266e-17, 1.15e-18, 1.544e-17, 5.382e-17, 1.264e-17, 1.088e-17, 2.002e-19, 5.658e-17, 2.675e-17, 1.982e-17, 4.633e-17, 4.633e-17, 8.87e-17, 5.658e-17, 9.324e-17, 7.635e-17, 1.015e-23, 2.246e-12, 1.83e-11, 1.11e-11, 6.379e-11, 2.897e-13, 1.296e-12, 5.491e-11, 4.497e-11, 4.497e-11, 4.774e-13, 6.411e-12, 9.094e-12, 3.015e-11, 3.503e-11, 1.74e-11, 3.186e-12, 4.751e-12, 1.432e-12, 1.173e-12, 6.411e-12, 7.484e-13, 8.23e-12, 1.749e-12, 1.607e-15, 3.046e-13, 3.366e-13, 2.742e-12, 1.749e-12, 1.656e-11, 3.701e-12, 7.484e-13, 1.227e-11, 2.234e-11, 3.382e-14, 2.729e-11, 2.125e-11, 3.015e-11, 1.74e-11, 2.596e-11, 1.749e-12, 2.883e-12, 3.871e-11, 9.512e-11, 4.497e-11, 9.999e-11, 1.227e-11, 5.25e-12, 1.82e-10, 2.457e-10, 2.223e-10, 3.664e-10, 4.944e-10, 4.049e-10, 4.751e-12, 2.854e-10, 4.049e-10, 5.464e-10, 3.485e-10, 4.994e-12, 3.701e-12, 5.744e-10, 8.271e-13, 1.283e-10, 3.015e-11, 8.271e-13, 2.742e-12, 4.703e-10, 4.408e-17, 2.469e-11, 2.742e-12, 1.296e-12, 7.868e-13, 1.923e-11, 1.664e-12, 7.12e-13, 3.851e-10, 9.608e-13, 6.706e-11, 1.349e-10, 3.871e-11, 2.246e-12, 1.749e-12, 3.52e-12, 8.608e-11, 2.854e-10, 1.167e-11, 4.069e-11, 8.694e-13, 1.943e-13, 2.223e-10, 4.299e-12, 1.061e-12, 2.635e-14, 1.432e-12, 1.179e-13, 5.574e-14, 4.32e-13, 1.647e-10, 8.694e-13, 1.296e-12, 2.481e-12, 9.56e-12, 4.11e-13, 5.83e-13, 4.32e-13, 2.897e-13, 4.11e-13, 1.362e-12, 2.136e-12, 3.046e-13, 2.042e-13, 1.758e-13, 4.969e-11, 3.17e-11, 5.859e-14, 4.277e-11, 4.798e-14, 3.538e-13, 1.943e-13, 3.719e-13, 1.672e-13, 1.179e-13, 8.738e-14, 3.17e-11])





















#using scipy optimize to fit the accretion disc model to the data
'''
def model(M, z, AccR):
  return Flux2(Lsum(M, AccR, R_midpoints(M), Rin(M)), z)

popt = scipy.optimize.curve_fit(model, M81x*10e9, M81y)


print(popt)
'''


#print(np.log10(M81x), np.log10(M81y))
#print(len(M81x), len(M81y))


















































########## Matrix attempt for luminosity increase with mass and with accretion rate ##########


'''
# Plot the spectrum for each mass
plt.figure(figsize=(8, 6))
plt.xlim(9, 18)
plt.ylim(24, 31)


for M in MassList:
    LBol = Lsum(M, Ar, R_midpoints(M), Rin(M))  # Assume Accretion Rate = 0.1
    plt.plot(np.log10(freq), np.log10(freq * LBol), label=f"M = {M:.1e} M☉")  # Plot freq × Luminosity(freq)
    #print(scipy.integrate.trapezoid(L, freq))
    TotalL = scipy.integrate.trapezoid(LBol, freq)

plt.xlabel("Frequency (Hz)")
plt.ylabel(r"$\nu L_\nu$ (W)")
plt.title("Spectra of Accretion Disks for Different Black Hole Masses")
plt.legend(loc="best", fontsize=8)
#plt.show()


def FluxV(T):
  temps = T.reshape(-1, 1)
  freqs = VisibleF.reshape(1, -1)
  a = (2 * np.pi *const.h.to_value() * freqs**3)/(const.c.to_value()**2)
  b = (const.h.to_value() * freqs)/(const.k_B.to_value() * temps)
  c = np.exp(b)
  Fv = a/(c-1)
  return Fv

def FluxX(T):
  temps = T.reshape(-1, 1)
  freqs = XrayF.reshape(1, -1)
  a = (2 * np.pi *const.h.to_value() * freqs**3)/(const.c.to_value()**2)
  b = (const.h.to_value() * freqs)/(const.k_B.to_value() * temps)
  c = np.exp(b)
  Fv = a/(c-1)
  return Fv

for M in MassList: 
  LVis = np.sum(FluxV(Temp2(M, Ar, R_midpoints(M), Rin(M))) * Area(R_midpoints(M)), axis = 0)
  TotLVis = scipy.integrate.trapezoid(LVis, VisibleF)

for M in MassList:
  LXray = np.sum(FluxX(Temp2(M, Ar, R_midpoints(M), Rin(M))) * Area(R_midpoints(M)), axis = 0)
  TotLXray = scipy.integrate.trapezoid(LXray, XrayF)

'''




VisibleF = (np.logspace(np.log10(420), np.log10(564), 1000)/10e9)*const.c.to_value()
XrayF = (np.logspace(np.log10(1), np.log10(100), 1000)/10e12)*const.c.to_value()
MassList = np.logspace(3, 8, 100) * const.M_sun.to_value() 
Mass = 10*const.M_sun.to_value()
Mass2 = 100*const.M_sun.to_value()
Mass3 = 5*const.M_sun.to_value()
Ar = 1e-8*const.M_sun.to_value()/(365.25*24*60*60)
Ar2 = 1e-7*const.M_sun.to_value()/(365.25*24*60*60)
Ar3 = 1e-9*const.M_sun.to_value()/(365.25*24*60*60)
AccList = np.logspace(np.log10(1e4), np.log10(1e18), 100)
Halpha = const.c.to_value()*656.3e-9
Hbeta = const.c.to_value()*486.1e-9
OIII = const.c.to_value()*500.7e-9
HeII = const.c.to_value()*468.6e-9



############################ 2D cases ############################
'''
################## Changing mass ################## 

MmBol = np.empty(len(MassList))
MmBol2 = np.empty(len(MassList))
MmBol3 = np.empty(len(MassList))
MmBol_Edd = np.empty(len(MassList))
MmBol_solar = np.empty(len(MassList))
MmVis= np.zeros(len(MassList)) 
MmVis2= np.zeros(len(MassList))
MmVis3= np.zeros(len(MassList))
MmXray = np.zeros(len(MassList))
MmXray2 = np.zeros(len(MassList))
MmXray3 = np.zeros(len(MassList))
MmVX = np.zeros(len(MassList))
MmVX2 = np.zeros(len(MassList))
MmVX3 = np.zeros(len(MassList))

for i, M in enumerate(MassList):
    LBol = Lsum(M, Ar, R_midpoints(M), Rin(M))  
    LBol2 = Lsum(M, Ar2, R_midpoints(M), Rin(M))
    LBol3 = Lsum(M, Ar3, R_midpoints(M), Rin(M))
    MmBol[i] = scipy.integrate.trapezoid(LBol, freq)
    MmBol2[i] = scipy.integrate.trapezoid(LBol2, freq)
    MmBol3[i] = scipy.integrate.trapezoid(LBol3, freq)
    MmBol_Edd[i] = MmBol[i]/L_edd(M)
    MmBol_solar[i] = MmBol[i]/const.L_sun.to_value()

    maskVis = (freq >= 4.3e14) & (freq <= 7.5e14)
    maskXray = (freq >= 1e16) & (freq <= 1e19)

    LVis = scipy.integrate.trapezoid(LBol[maskVis], freq[maskVis])
    LVis2 = scipy.integrate.trapezoid(LBol2[maskVis], freq[maskVis])
    LVis3 = scipy.integrate.trapezoid(LBol3[maskVis], freq[maskVis])
    LXray = scipy.integrate.trapezoid(LBol[maskXray], freq[maskXray])
    LXray2 = scipy.integrate.trapezoid(LBol2[maskXray], freq[maskXray])
    LXray3 = scipy.integrate.trapezoid(LBol3[maskXray], freq[maskXray])
    LHalpha = scipy.integrate.trapezoid(LBol[(freq >= Halpha - 10e-9) & (freq <= Halpha + 10e-9)], freq[(freq >= Halpha - 10e-9) & (freq <= Halpha + 10e-9)])
    LHbeta = scipy.integrate.trapezoid(LBol[(freq >= Hbeta - 10e-9) & (freq <= Hbeta + 10e-9)], freq[(freq >= Hbeta - 10e-9) & (freq <= Hbeta + 10e-9)])
    LOIII = scipy.integrate.trapezoid(LBol[(freq >= OIII - 10e-9) & (freq <= OIII + 10e-9)], freq[(freq >= OIII - 10e-9) & (freq <= OIII + 10e-9)])
    LHeII = scipy.integrate.trapezoid(LBol[(freq >= HeII - 10e-9) & (freq <= HeII + 10e-9)], freq[(freq >= HeII - 10e-9) & (freq <= HeII + 10e-9)])

    MmVis[i] = LVis 
    MmXray[i] = LXray
    MmVX[i] = LVis/LXray if LXray > 0 else np.nan
    MmVis2[i] = LVis2
    MmXray2[i] = LXray2
    MmVX2[i] = LVis2/LXray2 if LXray2 > 0 else np.nan
    MmVis3[i] = LVis3
    MmXray3[i] = LXray3
    MmVX3[i] = LVis3/LXray3 if LXray3 > 0 else np.nan


mLBolmatrix = np.vstack((MassList, MmBol)).T
mLVismatrix = np.vstack((MassList, MmVis)).T
mLXraymatrix = np.vstack((MassList, MmXray)).T
mLVXmatrix = np.vstack((MassList, MmVX)).T

###################################################


################## Changing accretion rate ##################

AmBol = np.empty(len(AccList))
AmBol2 = np.empty(len(AccList))
AmBol3 = np.empty(len(AccList))
AmBol_solar = np.empty(len(AccList))
AmBol_Edd = np.empty(len(AccList))
AmVis= np.zeros(len(AccList)) 
AmVis2= np.zeros(len(AccList))
AmVis3= np.zeros(len(AccList))
AmXray = np.zeros(len(AccList))
AmXray2 = np.zeros(len(AccList))
AmXray3 = np.zeros(len(AccList))
AmVX = np.zeros(len(AccList))
AmVX2 = np.zeros(len(AccList))
AmVX3 = np.zeros(len(AccList))


AmHalpha = np.zeros(len(AccList))
AmHbeta = np.zeros(len(AccList))
AmOIII = np.zeros(len(AccList))
AmHeII = np.zeros(len(AccList))
AmHB_OIII = np.zeros(len(AccList))


for j, Acc in enumerate(AccList):
  LBol = Lsum(Mass, Acc, R_midpoints(Mass), Rin(Mass))
  LBol2 = Lsum(Mass2, Acc, R_midpoints(Mass), Rin(Mass))
  LBol3 = Lsum(Mass3, Acc, R_midpoints(Mass), Rin(Mass))
  AmBol[j] = scipy.integrate.trapezoid(LBol, freq)
  AmBol2[j] = scipy.integrate.trapezoid(LBol2, freq)
  AmBol3[j] = scipy.integrate.trapezoid(LBol3, freq)
  AmBol_Edd[j] = AmBol[j]/L_edd(Mass)
  AmBol_solar[j] = AmBol[j]/const.L_sun.to_value()

  maskVis = (freq >= 4.3e14) & (freq <= 7.5e14)
  maskXray = (freq >= 1e16) & (freq <= 1e19)

  LVis = scipy.integrate.trapezoid(LBol[maskVis], freq[maskVis])
  LVis2 = scipy.integrate.trapezoid(LBol2[maskVis], freq[maskVis])
  LVis3 = scipy.integrate.trapezoid(LBol3[maskVis], freq[maskVis])
  LXray = scipy.integrate.trapezoid(LBol[maskXray], freq[maskXray])
  LXray2 = scipy.integrate.trapezoid(LBol2[maskXray], freq[maskXray])
  LXray3 = scipy.integrate.trapezoid(LBol3[maskXray], freq[maskXray])
  LHalpha = scipy.integrate.trapezoid(LBol[(freq >= Halpha - 10e-9) & (freq <= Halpha + 10e-9)], freq[(freq >= Halpha - 10e-9) & (freq <= Halpha + 10e-9)])
  LHbeta = scipy.integrate.trapezoid(LBol[(freq >= Hbeta - 10e-9) & (freq <= Hbeta + 10e-9)], freq[(freq >= Hbeta - 10e-9) & (freq <= Hbeta + 10e-9)])
  LOIII = scipy.integrate.trapezoid(LBol[(freq >= OIII - 10e-9) & (freq <= OIII + 10e-9)], freq[(freq >= OIII - 10e-9) & (freq <= OIII + 10e-9)])
  LHeII = scipy.integrate.trapezoid(LBol[(freq >= HeII - 10e-9) & (freq <= HeII + 10e-9)], freq[(freq >= HeII - 10e-9) & (freq <= HeII + 10e-9)])

  AmVis[j] = LVis 
  AmXray[j] = LXray
  AmVX[j] = LVis/LXray if LXray > 0 else np.nan
  AmVis2[j] = LVis2
  AmXray2[j] = LXray2
  AmVX2[j] = LVis2/LXray2 if LXray2 > 0 else np.nan
  AmVis3[j] = LVis3
  AmXray3[j] = LXray3
  AmVX3[j] = LVis3/LXray3 if LXray3 > 0 else np.nan
  AmHalpha[j] = LHalpha
  AmHbeta[j] = LHbeta
  AmOIII[j] = LOIII
  AmHeII[j] = LHeII
  AmHB_OIII[j] = LHbeta/LOIII if LOIII > 0 else np.nan


aLBolematrix = np.vstack((AccList, AmBol)).T
aLVismatrix = np.vstack((AccList, AmVis)).T
aLXraymatrix = np.vstack((AccList, AmXray)).T
aLVXmatrix = np.vstack((AccList, AmVX)).T

###########################################################


plt.figure(figsize=(8, 6))
#plt.plot(np.log10(AccList/AccR_Edd(Mass)), np.log10(AmVis), linestyle='-', label='Visible')
#plt.plot(np.log10(AccList/AccR_Edd(Mass)), np.log10(AmXray), linestyle='-', label='X-ray')
#plt.plot(np.log10(AccList/const.M_sun.to_value()), np.log10(AmBol_solar), linestyle='-', label='Total Luminosity')
plt.plot(np.log10(AccList/AccR_Edd(Mass)), np.log10(AmVX), linestyle='-', color = 'black', label = r'10$M_{\odot}$', linewidth = 2)
plt.plot(np.log10(AccList/AccR_Edd(Mass2)), np.log10(AmVX2), linestyle='--', color = 'red', label = r'100$M_{\odot}$', linewidth = 2)
plt.plot(np.log10(AccList/AccR_Edd(Mass3)), np.log10(AmVX3), linestyle=':', color = 'blue', label = r'5$M_{\odot}$', linewidth = 2)
plt.xlabel(r"Log($\dot{M}$) [$\dot{M}$$_{Edd}$]", fontsize = 18)
plt.ylabel(r"Log($L_{V}$/$L_{X}$)", fontsize = 18)
plt.legend(fontsize=18, frameon=True, loc='best')
plt.tick_params(axis='both', which='major', labelsize=16, length=7, width=2, top=True, right=True) 
plt.tick_params(axis='both', which='minor', labelsize=10, length=4, width=1.5, top=True, right=True) 
plt.minorticks_on()
plt.tight_layout()
plt.savefig('Luminosity_Ratio_vs_Acccretion_Rate.png', dpi=300, bbox_inches='tight')

#            Mass=10$M_{\odot}$, AccR=10$^{-7}$$\dot{M}$$_{Edd}$') #5x10$^{-15}$$M_{\odot}$Yr$^{-1}$


######### Limit the lines so that they all start and stop at the same x-axis value #########


plt.figure(figsize=(8, 6))
#plt.plot(np.log10(MassList/const.M_sun.to_value()), np.log10(MmVis), linestyle='-', label='Visible')
#plt.plot(np.log10(MassList/const.M_sun.to_value()), np.log10(MmXray), linestyle='-', label='X-ray')
#plt.plot(np.log10(MassList/const.M_sun.to_value()), np.log10(MmBol_solar), linestyle='-', label='Total Luminosity')
plt.plot(np.log10(MassList/const.M_sun.to_value()), np.log10(MmVX), linestyle='-', color = 'black', label = r'10$^{-8}$$M_{\odot}$Yr$^{-1}$', linewidth = 2)
plt.plot(np.log10(MassList/const.M_sun.to_value()), np.log10(MmVX2), linestyle='--', color = 'red', label = r'10$^{-7}$$M_{\odot}$Yr$^{-1}$', linewidth = 2)
plt.plot(np.log10(MassList/const.M_sun.to_value()), np.log10(MmVX3), linestyle=':', color = 'blue', label = r'10$^{-9}$$M_{\odot}$Yr$^{-1}$', linewidth = 2)
plt.xlabel(r"Log(Mass) [$M_{\odot}$]", fontsize = 18)
plt.ylabel(r"Log($L_{V}$/$L_{X}$)", fontsize = 18)
plt.legend(fontsize=18, frameon=True, loc='best')
plt.tick_params(axis='both', which='major', labelsize=16, length=7, width=2, top=True, right=True) 
plt.tick_params(axis='both', which='minor', labelsize=10, length=4, width=1.5, top=True, right=True) 
plt.minorticks_on()
plt.tight_layout()
plt.savefig('Luminosity_Ratio_vs_Mass.png', dpi=300, bbox_inches='tight')
#plt.show()
'''





'''
for acc, lum in aLBolematrix:
    print(f"{acc:.2e}   |   {lum:.2e}")

for acc, lum in aLVismatrix:
    print(f"{acc:.2e}   |   {lum:.2e}")

for acc, lum in aLXraymatrix:
    print(f"{acc:.2e}   |   {lum:.2e}")


print("Mass (M☉)   |   Total Luminosity (erg/s)")
print("----------------------------------------")

for mass, lum in mLVismatrix:
    print(f"{mass/const.M_sun.to_value():.2e}   |   {lum:.2e}")

for mass, lum in mLXraymatrix:
    print(f"{mass/const.M_sun.to_value():.2e}   |   {lum:.2e}")

for mass, lum in mLBolmatrix:
    print(f"{mass/const.M_sun.to_value():.2e}   |   {lum:.2e}")
'''

########################################################################################












'''

#################################### 3D cases ##########################################
mBolematrix = np.zeros((len(MassList), len(AccList)))
mVismatrix = np.zeros((len(MassList), len(AccList)))
mXraymatrix = np.zeros((len(MassList), len(AccList)))
mVisXray = np.zeros((len(MassList), len(AccList)))

for i, M in enumerate(MassList):
  for j, Acc in enumerate(AccList):
    L_tot = Lsum(M, Acc, R_midpoints(M), Rin(M))

    maskVis = (freq >= 4.3e14) & (freq <= 7.5e14)
    maskXray = (freq >= 1e16) & (freq <= 1e19)

    LVis = scipy.integrate.trapezoid(L_tot[maskVis], freq[maskVis])
    LXray = scipy.integrate.trapezoid(L_tot[maskXray], freq[maskXray])

    #mBolematrix[i, j] = scipy.integrate.trapezoid(L_tot, freq)
    #mVismatrix[i, j] = LVis
    #mXraymatrix[i, j] = LXray
    mVisXray[i, j] = LVis / LXray if LXray > 0 else np.nan



MassGrid, AccGrid = np.meshgrid(MassList, AccList, indexing='ij')

AccGridEdd = AccGrid/AccR_Edd(MassGrid)
AccGridSolar = (AccGrid/const.M_sun.to_value())*24*60*60*365.25

#print(AccR_Edd(MassGrid))
#for M in MassGrid:
#  L_Bol = Lsum(M, AccGrid, R_midpoints(M), Rin(M))

#L_Bol = Lsum(MassGrid, Ar, R_midpoints(MassGrid), Rin(MassGrid))


#print(L_Bol)


#Bol = np.vstack((MassGrid.ravel(), AccGrid.ravel(), mBolematrix.ravel())).T
#Vis = np.vstack((MassGrid.ravel(), AccGrid.ravel(), mVismatrix.ravel())).T
#Xray = np.vstack((MassGrid.ravel(), AccGrid.ravel(), mXraymatrix.ravel())).T
VisXray = np.vstack((MassGrid.ravel(), AccGridEdd.ravel(), mVisXray.ravel())).T
#LBols = np.vstack((MassGrid.ravel(), AccGrid.ravel(), L_Bol.ravel())).T


#h = np.vstack((MassGrid.ravel(), AccGrid.ravel(), mVisXray.ravel())).T

#for mass, acc, visxray in Bol:
#  print(f"{mass/const.M_sun.to_value():.2e}   |   {acc:.2e}   |   {visxray:.2e}")#

#for mass, acc, visxray in Vis:
#  print(f"{mass:.2e}   |   {acc:.2e}   |   {visxray:.2e}")

#for mass, acc, visxray in Xray:
#  print(f"{mass:.2e}   |   {acc:.2e}   |   {visxray:.2e}")

#for mass, acc, visxray in VisXray:
#  print(f"{mass:.2e}   |   {acc:.2e}   |   {visxray:.2e}")


small_value = 1e-30
mVisXray_safe = np.maximum(mVisXray, small_value)  # Avoid log10(0)
mVisXray_safe = np.minimum(mVisXray_safe, 10e150)


plt.figure(figsize=(8, 6))

cmap = plt.get_cmap("inferno")
plt.pcolormesh(np.log10(MassGrid/const.M_sun.to_value()), np.log10(AccGridSolar), np.log10(mVisXray), shading='nearest', cmap=cmap) #shading options: flat, gouraud, nearest, auto
cbar = plt.colorbar()
plt.ylabel(r'Log(M) ($M_{\odot}$)', fontsize = 18)
cbar.set_label(r"Log($L_{V}$/$L_{X}$)", fontsize = 18)
plt.xlabel(r"Log($\dot{M}$) ($M_{\odot}$Yr$^{-1}$)", fontsize = 18)
plt.tick_params(axis='both', which='major', labelsize=16, length=7, width=2, top=True, right=True)
plt.tick_params(axis='both', which='minor', labelsize=10, length=4, width=1.5, top=True, right=True)
plt.minorticks_on()
plt.tight_layout()
plt.savefig('Accretion_Rate.png', dpi=300, bbox_inches='tight')
'''




'''
rc('font',size=28)
rc('font',family='serif')
rc('axes',labelsize=32)

# 3D Plot 
fig = plt.figure(figsize=(13.5,8))
ax = fig.add_subplot(121, projection='3d')
plt.subplots_adjust(left=0, right=0.90, top=1, bottom=0, wspace=0.22)
#ax.plot_surface(np.log10(MassGrid/const.M_sun.to_value()), np.log10(AccGrid), np.log10(mBolematrix))
#ax.plot_surface(np.log10(MassGrid), np.log10(AccGrid), np.log10(mVismatrix), cmap='viridis')
#ax.plot_surface(np.log10(MassGrid), np.log10(AccGrid), np.log10(mXraymatrix), cmap='viridis')
sur = ax.plot_surface(np.log10(MassGrid/const.M_sun.to_value()), np.log10(AccGridSolar), np.log10(mVisXray_safe), cmap='jet')
cont = (np.log10(mVisXray_safe).min()-np.log10(mVisXray_safe).max())*0.4
cset = ax.contour(np.log10(MassGrid/const.M_sun.to_value()),np.log10(AccGridSolar),np.log10(mVisXray_safe).reshape(np.log10(MassGrid/const.M_sun.to_value()).shape),zdir='z',offset=cont)
#ax.plot3D(np.log10(AccList), np.log10(AmVis), np.log10(AmXray))#, cmap='viridis')
ax1 = fig.add_subplot(122, projection='3d')
sur2 = ax1.plot_surface(np.log10(MassGrid/const.M_sun.to_value()), np.log10(AccGridEdd), np.log10(mVisXray_safe), cmap='viridis')
cont2 = (np.log10(mVisXray_safe).min()-np.log10(mVisXray_safe).max())*0.4
cset2 = ax1.contour(np.log10(MassGrid/const.M_sun.to_value()),np.log10(AccGridEdd),np.log10(mVisXray_safe).reshape(np.log10(MassGrid/const.M_sun.to_value()).shape),zdir='z',offset=cont2)



ax.xaxis.set_rotate_label(False)
ax.yaxis.set_rotate_label(False)
ax.zaxis.set_rotate_label(False)
ax.view_init(elev=10, azim=135)
ax.grid(False)
ax.xaxis.pane.set_edgecolor('black')
ax.yaxis.pane.set_edgecolor('black')
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
[t.set_va('center') for t in ax.get_yticklabels()]
[t.set_ha('left') for t in ax.get_yticklabels()]
[t.set_va('center') for t in ax.get_xticklabels()]
[t.set_ha('right') for t in ax.get_xticklabels()]
[t.set_va('center') for t in ax.get_zticklabels()]
[t.set_ha('left') for t in ax.get_zticklabels()]
ax.xaxis._axinfo['tick']['inward_factor'] = 0
ax.xaxis._axinfo['tick']['outward_factor'] = 0.4
ax.yaxis._axinfo['tick']['inward_factor'] = 0
ax.yaxis._axinfo['tick']['outward_factor'] = 0.4
ax.zaxis._axinfo['tick']['inward_factor'] = 0
ax.zaxis._axinfo['tick']['outward_factor'] = 0.4
ax.zaxis._axinfo['tick']['outward_factor'] = 0.4

ax.set_ylim(-19, -4)
ax.set_xlabel(r'Log(Mass) ($M_{\odot}$)', fontsize = 16)
ax.set_zlabel(r"Log($L_{V}$/$L_{X}$)", fontsize = 16)
ax.set_ylabel(r"Log(Accretion Rate) ($M_{\odot}$Yr$^{-1}$)", fontsize = 16)
#ax.tick_params(axis='x', which='major', labelsize=14, length=6, width=1.5, top=True, right=True, left=True, bottom=True)
#ax.tick_params(axis='x', which='minor', labelsize=10, length=3, width=1, top=True, right=True, left=True, bottom=True)
#ax.tick_params(axis='y', which='major', labelsize=14, length=6, width=1.5, top=True, right=True, left=True, bottom=True)
#ax.tick_params(axis='y', which='minor', labelsize=10, length=3, width=1, top=True, right=True, left=True, bottom=True)
#ax.tick_params(axis='z', which='major', labelsize=14, length=6, width=1.5, top=True, right=True, left=True, bottom=True)
#ax.tick_params(axis='z', which='minor', labelsize=10, length=3, width=1, top=True, right=True, left=True, bottom=True)
#ax.minorticks_on()

ax1.set_xlabel(r'Log(Mass) ($M_{\odot}$)', fontsize = 16)
ax1.set_zlabel(r"Log($L_{V}$/$L_{X}$)", fontsize = 16)
ax1.set_ylabel(r"Log(Accretion Rate) ($M_{\odot}$Yr$^{-1}$)", fontsize = 16)

ax1.xaxis.set_rotate_label(False)
ax1.yaxis.set_rotate_label(False)
ax1.zaxis.set_rotate_label(False)
ax1.view_init(elev=10, azim=135)
ax1.grid(False)
ax1.xaxis.pane.set_edgecolor('black')
ax1.yaxis.pane.set_edgecolor('black')
ax1.xaxis.pane.fill = False
ax1.yaxis.pane.fill = False
ax1.zaxis.pane.fill = False
[t.set_va('center') for t in ax1.get_yticklabels()]
[t.set_ha('left') for t in ax1.get_yticklabels()]
[t.set_va('center') for t in ax1.get_xticklabels()]
[t.set_ha('right') for t in ax1.get_xticklabels()]
[t.set_va('center') for t in ax1.get_zticklabels()]
[t.set_ha('left') for t in ax1.get_zticklabels()]
ax1.xaxis._axinfo['tick']['inward_factor'] = 0
ax1.xaxis._axinfo['tick']['outward_factor'] = 0.4
ax1.yaxis._axinfo['tick']['inward_factor'] = 0
ax1.yaxis._axinfo['tick']['outward_factor'] = 0.4
ax1.zaxis._axinfo['tick']['inward_factor'] = 0
ax1.zaxis._axinfo['tick']['outward_factor'] = 0.4
ax1.zaxis._axinfo['tick']['outward_factor'] = 0.4
'''

#plt.tight_layout()

#fig = plt.figure(figsize=(12, 9))
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(np.log10(MassGrid/const.M_sun.to_value()), np.log10(L_tot), np.log10(mVisXray), cmap='viridis')

#ax.set_zlim(0, 200)

#ax.set_xlabel('log10(Mass) [M☉]')
#ax.set_zlabel('log10(Luminosity [visible/x-ray])')
#ax.set_ylabel('log10(Accretion Rate) [M☉/yr]')
#plt.show()

'''
plt.figure(figsize=(8, 6))
#ax = fig.add_subplot(111, projection='3d')
plt.plot(np.log10(AmXray), np.log10(AmVis), linestyle='-', label='Visible/X-ray')
#plt.plot(AccList, AmVX)
plt.xlabel('Log(L X-ray)')
plt.ylabel('Log(L Optical)')
'''




########################################################################################






'''
#Milestone values
LumMilestone = (Flux(Temp2(MassBH, AccR, R_midpoints(MassBH), Rin(MassBH)))) * Area(R_midpoints(MassBH))  
LumMilestonesum = np.sum(LumMilestone, axis = 0)

LumSMBH1 = (Flux(Temp2(MassSMBH1, AccR, R_midpoints(MassSMBH1), Rin(MassSMBH1)))) * Area(R_midpoints(MassSMBH1))  
LumSMBH1sum = np.sum(LumSMBH1, axis = 0)

LumSMBH2 = (Flux(Temp2(MassSMBH2, AccR, R_midpoints(MassSMBH2), Rin(MassSMBH2)))) * Area(R_midpoints(MassSMBH2))  
LumSMBH2sum = np.sum(LumSMBH2, axis = 0)

LumSMBH3 = (Flux(Temp2(MassSMBH3, AccR, R_midpoints(MassSMBH3), Rin(MassSMBH3)))) * Area(R_midpoints(MassSMBH3))  
LumSMBH3sum = np.sum(LumSMBH3, axis = 0)



TotLumMilestone = scipy.integrate.trapezoid(LumMilestonesum, freq)*u.W
TotLumSMBH1 = scipy.integrate.trapezoid(LumSMBH1sum, freq)*u.W
TotLumSMBH2 = scipy.integrate.trapezoid(LumSMBH2sum, freq)*u.W
TotLumSMBH3 = scipy.integrate.trapezoid(LumSMBH3sum, freq)*u.W



MassList2 = np.logspace(3, 11, 100) * const.M_sun.to_value()
AccList2 = np.logspace(np.log10(10e4), np.log10(10e20), 100)

Mass = 10*const.M_sun.to_value()
Ar = 10**17


LvissM = np.zeros(len(MassList2))
LvissA = np.zeros(len(AccList2))
'''

'''

fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

for i, M in enumerate(MassList2):
  BolLums = Lsum(M, Ar, R_midpoints(M), Rin(M))

  maskVis = (freq >= 4.3e14) & (freq <= 7.5e14)

  LVis = scipy.integrate.trapezoid(BolLums[maskVis], freq[maskVis])
  #print(LVis)
  LvissM[i] = LVis

ax1.plot(np.log10(MassList2/const.M_sun.to_value()), np.log10(LvissM), linestyle='-', label=f"M = {i:.1e} M☉")  



fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

for j, A in enumerate(AccList2):
  BolLums = Lsum(Mass, A, R_midpoints(Mass), Rin(Mass))

  maskVis = (freq >= 4.3e14) & (freq <= 7.5e14)

  LVis = scipy.integrate.trapezoid(BolLums[maskVis], freq[maskVis])
  LvissA[j] = LVis

ax1.plot(np.log10(AccList2/AccR_Edd(Mass)), np.log10(LvissA), linestyle='-', label=f"Acc = {j:.1e} M☉/s")

'''







###########################  SED plot USE THIS ###########################

'''
LumDAcc1sum = Lsum(MassBH, 0.0001*AccR_Edd(MassBH), R_midpoints(MassBH), Rin(MassBH))
LumDAcc2sum = Lsum(MassBH, 0.1*AccR_Edd(MassBH), R_midpoints(MassBH), Rin(MassBH))
LumMilestonesum = Lsum(MassBH, 0.0000001*AccR_Edd(MassBH), R_midpoints(MassBH), Rin(MassBH))
LumSMBH1sum = Lsum(MassSMBH1, 0.0000001*AccR_Edd(MassBH), R_midpoints(MassSMBH1), Rin(MassSMBH1))
LumSMBH2sum = Lsum(MassSMBH2, 0.0000001*AccR_Edd(MassBH), R_midpoints(MassSMBH2), Rin(MassSMBH2))









em_spectrum = {
    'Optical': (4e14 / 1e9, 7.5e14 / 1e9),
    'X-ray': (3e16 / 1e9, 3e19 / 1e9),
}





fig, ax1 = plt.subplots()
fig.set_size_inches(8,6)

#fig.patch.set_facecolor('#D5D5D5') 
#ax1.set_facecolor('#D5D5D5') 

ax1.set_xlim(-2, 9)
ax1.set_ylim(9, 32)

d, = ax1.plot(np.log10(freq/10e9), np.log10(freq * LumMilestonesum), linestyle='-', color='black', label=r'Mass=10$M_{\odot}$, AccR=10$^{-7}$$\dot{M}$$_{Edd}$', linewidth=2) #5x10$^{-15}$$M_{\odot}$Yr$^{-1}$')
#e, = ax1.plot(np.log10(freq), np.log10(freq * LumEddsum), linestyle='-', color='purple', label='Eddington Limit (10)')
a, = ax1.plot(np.log10(freq/10e9), np.log10(freq * LumSMBH1sum), linestyle='--', color='red', label=r'Increaing Mass +10$^4$$M_{\odot}$', linewidth=2)
b, = ax1.plot(np.log10(freq/10e9), np.log10(freq * LumSMBH2sum), linestyle='--', color='red', linewidth=2)#, label='10^8')
h, = ax1.plot(np.log10(freq/10e9), np.log10(freq * LumDAcc1sum), linestyle=':', color='blue', label=r'Increasing Accretion Rate +10$^{-3}$$\dot{M}$$_{Edd}$', linewidth=2)
i, = ax1.plot(np.log10(freq/10e9), np.log10(freq * LumDAcc2sum), linestyle=':', color='blue', linewidth=2)#, label='') 



#c, = ax1.plot(np.log10(freq), np.log10(freq * LumSMBH3sum), linestyle='-', color='red', label='10^11')
#f, = ax1.plot(np.log10(freq), np.log10(freq * LumSunMass2um), linestyle='-', color='brown', label='10^4')
#g, = ax1.plot(np.log10(freq), np.log10(Lsum4(MassBH, AccR, r(MassBH), rin(MassBH))), linestyle='-', color='orange', label='scaled')


#dAccr, = ax1.plot(np.log10(freq), np.log10(freq * LumSunMass2um), linestyle='-', color='green')
#dAccrSpin, = ax1.plot(np.log10(freq), np.log10(freq * LumSunMass2umSpin), linestyle=':', color='green')



original = plt.Line2D([0], [0], color='black', linestyle='-', label=r'Mass=10$M_{\odot}$, AccR=10$^{15}$$M_{\odot}$', linewidth=2)
mass = plt.Line2D([0], [0], color='red', linestyle='--', label=r'Increasing Mass +10$^4$$M_{\odot}$', linewidth=2)
accretion = plt.Line2D([0], [0], color='blue', linestyle=':', label=r'Increasing Accretion Rate +10$^{15}$$M_{\odot}$', linewidth=2)
mass_10msun = plt.Line2D([0], [0], color='blue', linestyle='-', label=r'Mass=10$M_{\odot}$', linewidth=2)
mass_1msun = plt.Line2D([0], [0], color='green', linestyle='-', label=r'Mass=10$^{4}$$M_{\odot}$', linewidth=2),

ax1.legend(handles=[d, a, h], fontsize=16, loc='upper right')

colors = {
    'Optical': 'gray',
    'X-ray': 'gray',
}

for band, (start, end) in em_spectrum.items():
    ax1.axvspan(np.log10(start), np.log10(end), color=colors[band], alpha=0.2)



ax1.set_ylabel(r"$\log_{10}(\nu L_\nu)$ (W)", fontsize=20)
ax1.set_xlabel(r'$\log_{10}(\nu)$ (GHz)', fontsize=20)
ax1.minorticks_on()
ax1.tick_params(axis='x', which='major', bottom=True, top=False, labelsize=16, length=7, width=2)
ax1.tick_params(axis='x', which='minor', bottom=True, top=False, labelsize=10, length=4, width=1.5)
ax1.tick_params(axis='y', which='major', right=True, top=False, labelsize=16, length=7, width=2)
ax1.tick_params(axis='y', which='minor', right=True, top=False, labelsize=10, length=4, width=1.5)
ax1.legend(fontsize=14, loc='upper left')

# Create a secondary x-axis
ax2 = ax1.twiny()

freq_ticks = ax1.get_xticks()

# Convert frequency ticks to wavelength ticks
wavelength_ticks = np.log10(3*(const.c.value / (10**freq_ticks)))

# Create custom labels
custom_labels = [f'{tick:.0f}' if tick != 0 else '-e' for tick in wavelength_ticks]

# Set the secondary x-axis limits and labels
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(freq_ticks)
ax2.set_xticklabels(custom_labels)
ax2.set_xlabel(r'$\log_{10}(\lambda)$ m', fontsize=20)
ax2.minorticks_on()
ax2.tick_params(axis='x', which='major', top=True, labelsize=16, length=7, width=2)
ax2.tick_params(axis='x', which='minor', right=True, top=True, labelsize=16, length=4, width=1.5)


peak_indices = {
    'MileStone': np.argmax(LumSMBH1sum),
    'MaxSpin': np.argmax(LumSMBH2sum),
    
}

peak_frequencies = {key: freq[idx] for key, idx in peak_indices.items()}
peak_y_values = {key: np.log10(freq[peak_indices[key]] * value[peak_indices[key]]) for key, value in zip(peak_indices.keys(),
[LumSMBH1sum, LumSMBH2sum])}


#print("Peak frequencies (Hz) and corresponding Luminosities (W):")
#for key in peak_frequencies:
#    print(f"{key}: Frequency = {peak_frequencies[key]:.2e} Hz, Luminosity = {10**peak_y_values[key]:.2e}")

plt.savefig('SED_plot.png', dpi=300, bbox_inches='tight')
'''
plt.tight_layout()  

plt.show()





















'''
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

ax1.set_xlim(-0.3, 9)
ax1.set_ylim(-21, -10)

#ax1.plot(np.log10((freq/10e9)), np.log10(freq * F3C273), color='black', linestyle='-', label='Radio-loud 3C 273')



#Plot of Seyfert 1 NGC 6814
FNGC6814 = Flux2(Lsum(3.9e6*const.M_sun.to_value(), AcRate(Luminosity(1.23e-12, 0.00579177), 3.9e6*const.M_sun.to_value()),
R_midpoints(3.9e6*const.M_sun.to_value()), Rin(3.9e6*const.M_sun.to_value())), 0.00579177)



#ax1.plot(np.log10((freq/10e9)), np.log10(freq * FNGC6814), color='green', linestyle='-', label='Seyfert 1 NGC 6814')


#Plot of Seyfert 2 M81
ax1.plot(np.log10((freq/10e9)), np.log10(freq * F5M81), color='purple', linestyle='-', label='Seyfer 2 M81')
ax1.scatter(np.log10(M81x), np.log10(M81y), color='black', marker='o', label='Data points')

#print(np.sum(F3C273), np.sum(FNGC6814), np.sum(F5M81))

colors = {
    'Radio': 'lightblue',
    'Microwave': 'lightgreen',
    'Infrared': 'lightcoral',
    'Optical': 'lightyellow',
    'Ultraviolet': 'lightpink',
    'X-ray': 'lightgray',
    'Gamma-ray': 'lightcyan'
}

for band, (start, end) in em_spectrum.items():
    ax1.axvspan(np.log10(start), np.log10(end), color=colors[band], alpha=0.3)


num_points = 1000
random_x = np.random.uniform(-0.3, 7, num_points)
random_y = np.random.uniform(-21, -10, num_points)

# Plot the random points (invisible)
#ax1.errorbar(random_x, random_y, color='red', marker='o', yerr = 0.5)


ax1.set_ylabel(r'$\log_{10}(vF_v)$ W m$^{-2}$', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(v)$ GHz', fontsize=16)

ax1.legend(fontsize=12, loc='upper left')
plt.show()
'''


'''
# Calculate residuals
observed_y = np.log10(freq * F3C273)
predicted_y = np.polyval(np.polyfit(np.log10((freq/10e9)), observed_y, 1), np.log10((freq/10e9)))
residuals = observed_y - predicted_y

# Plot residuals
fig, ax2 = plt.subplots()
fig.set_size_inches(10, 6)

ax2.scatter(np.log10((freq/10e9)), residuals, color='blue', marker='o', label='Residuals')
ax2.axhline(0, color='red', linestyle='--')

ax2.set_ylabel('Residuals', fontsize=16)
ax2.set_xlabel(r'$\log_{10}(v)$ GHz', fontsize=16)

ax2.legend(fontsize=12, loc='upper left')
plt.show()

print(scipy.stats.chisquare(freq * F3C273, random_y))
'''

#Plot of M81 perhaps (Seyfert 2 AGN)
'''
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

ax1.set_xlim(-0.3, 7)
ax1.set_ylim(-21, -10)

ax1.plot(np.log10((freq/10e9)), np.log10(freq * F1M81), color='blue', linestyle='-', label='M81')
ax1.plot(np.log10((freq/10e9)), np.log10(freq * F2M81), color='black', linestyle='-', label='Using eq accretion')
ax1.plot(np.log10((freq/10e9)), np.log10(freq * F3M81), color='black', linestyle='--', label='Using Frank')
ax1.plot(np.log10((freq/10e9)), np.log10(freq * F4M81), color='blue', linestyle='--', label='Using Frank')
ax1.plot(np.log10((freq/10e9)), np.log10(freq * F5M81), color='purple', linestyle='-', label='L from S,S 1976')
ax1.plot(np.log10((freq/10e9)), np.log10(freq * F6M81), color='purple', linestyle='--', label='L from S,S 1976 and Using Frank')


colors = {
    'Radio': 'lightblue',
    'Microwave': 'lightgreen',
    'Infrared': 'lightcoral',
    'Optical': 'lightyellow',
    'Ultraviolet': 'lightpink',
    'X-ray': 'lightgray',
    'Gamma-ray': 'lightcyan'
}

for band, (start, end) in em_spectrum.items():
    ax1.axvspan(np.log10(start), np.log10(end), color=colors[band], alpha=0.3)


ax1.set_ylabel(r'$\log_{10}(vF_v)$ W m$^{-2}$', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(v)$ GHz', fontsize=16)

ax1.legend(fontsize=12, loc='upper left')
plt.show()
'''

#peak frequencies and luminosities
'''
peak_indices = {
    'M81': np.argmax(freq * F1M81),
    'Using eq accretion': np.argmax(freq * F2M81),
    'Using Frank': np.argmax(freq * F3M81),
    'Using Frank 2': np.argmax(freq * F4M81),
    'L from S,S 1976': np.argmax(freq * F5M81),
    'L from S,S 1976 and Using Frank': np.argmax(freq * F6M81)
}

peak_y_values = {key: np.log10(freq[idx] * value[idx]) for key, idx, value in zip(peak_indices.keys(), peak_indices.values(), [F1M81, F2M81, F3M81, F4M81, F5M81, F6M81])}

print("Peak frequencies (Hz) and corresponding Luminosities (W):")
for key in peak_y_values:
    print(f"{key}: Frequency = {freq[peak_indices[key]]:.2e} Hz, Luminosity = {10**peak_y_values[key]:.2e}")

'''


#plot of peak frequencies and mass of agn
'''
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

#ax1.set_xlim(6, 10)
#ax1.set_ylim(0, 12)

masses = np.logspace(6, 10, 100)
peak_freqs = np.zeros_like(masses)
peak_freqs2 = np.zeros_like(masses)
peak_freqs3 = np.zeros_like(masses)
peak_freqs4 = np.zeros_like(masses)
for i, mass in enumerate(masses):
    F = Flux2(Lsum2(mass*const.M_sun.to_value(), AcRate(Luminosity(5.19e-11, 0.005), mass*const.M_sun.to_value()),
    R_midpoints(mass*const.M_sun.to_value()), Rin(mass*const.M_sun.to_value())), 0.005)
    peak_freqs[i] = freq[np.argmax(freq * F)]
    F2 = Flux2(Lsum2(mass*const.M_sun.to_value(), AccR_Edd(mass*const.M_sun.to_value()),
    R_midpoints(mass*const.M_sun.to_value()), Rin(mass*const.M_sun.to_value())), 0.005)
    peak_freqs2[i] = freq[np.argmax(freq * F2)]
    F3 = Flux2(Lsum2(mass*const.M_sun.to_value(), AcRate(Luminosity(5.19e-11, 0.05), mass*const.M_sun.to_value()),
    R_midpoints(mass*const.M_sun.to_value()), Rin(mass*const.M_sun.to_value())), 0.05)
    peak_freqs3[i] = freq[np.argmax(freq * F3)]
    F4 = Flux2(Lsum2(mass*const.M_sun.to_value(), AcRate(Luminosity(5.19e-9, 0.005), mass*const.M_sun.to_value()),
    R_midpoints(mass*const.M_sun.to_value()), Rin(mass*const.M_sun.to_value())), 0.005)
    peak_freqs4[i] = freq[np.argmax(freq * F4)]

ax1.plot(np.log10(masses), np.log10(peak_freqs), color='purple', linestyle='-', label='S, S Accretion rate')
ax1.plot(np.log10(masses), np.log10(peak_freqs2), color='black', linestyle='-', label='Eddington limit')
ax1.plot(np.log10(masses), np.log10(peak_freqs3), color='blue', linestyle='-', label='Larger redshift')
ax1.plot(np.log10(masses), np.log10(peak_freqs4), color='red', linestyle='-', label='Larger luminosity')

ax1.set_ylabel(r'$\log_{10}(\nu_{\text{peak}})$ Hz', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(M_{\text{AGN}})$ M$_\odot$', fontsize=16)

plt.show()
'''



#Luminosity vs mass
'''
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

#ax1.set_xlim(6, 10)
#ax1.set_ylim(0, 12)

masses = np.logspace(6, 10, 100)
Ls = np.zeros_like(masses)
Ls2 = np.zeros_like(masses)
Ls3 = np.zeros_like(masses)
Ls4 = np.zeros_like(masses)
for i, mass in enumerate(masses):
    L = Lsum2(mass*const.M_sun.to_value(), AcRate(Luminosity(5.19e-11, 0.005), mass*const.M_sun.to_value()),
    R_midpoints(mass*const.M_sun.to_value()), Rin(mass*const.M_sun.to_value()))
    Ls[i] = Luminosity(L, 0.005)
    L2 = Lsum2(mass*const.M_sun.to_value(), AccR_Edd(mass*const.M_sun.to_value()),
    R_midpoints(mass*const.M_sun.to_value()), Rin(mass*const.M_sun.to_value()))
    Ls2[i] = Luminosity(L2, 0.005)
    L3 = Lsum2(mass*const.M_sun.to_value(), AcRate(Luminosity(5.19e-11, 0.05), mass*const.M_sun.to_value()),
    R_midpoints(mass*const.M_sun.to_value()), Rin(mass*const.M_sun.to_value()))
    Ls3[i] = Luminosity(L3, 0.05)
    L4 = Lsum2(mass*const.M_sun.to_value(), AcRate(Luminosity(5.19e-9, 0.005), mass*const.M_sun.to_value()),
    R_midpoints(mass*const.M_sun.to_value()), Rin(mass*const.M_sun.to_value()))
    Ls4[i] = Luminosity(L4, 0.005)
    

ax1.plot(np.log10(masses), np.log10(Ls), color='purple', linestyle='-', label='S, S Accretion rate')
ax1.plot(np.log10(masses), np.log10(Ls2), color='black', linestyle='-', label='Eddington limit')
ax1.plot(np.log10(masses), np.log10(Ls3), color='blue', linestyle='-', label='Larger redshift')
ax1.plot(np.log10(masses), np.log10(Ls4), color='red', linestyle='-', label='Larger luminosity')

ax1.set_ylabel(r'$\log_{10}(\nu_{\text{peak}})$ Hz', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(M_{\text{AGN}})$ M$_\odot$', fontsize=16)

plt.show()
'''






#def rho():


#def c_s():


#def Bondi_accr(M, rho, c_s):





#L in terms of B-band luminosity(erg s-1)/10e9*L_sun
#M in terms of mass/10e8*M_sun
#def Accretion(M, L):
#  return (L/(13.8 * M**(2/3)))**(3/2)



#Calculations using ergs
'''
def Terg(M, Ar, Radius, Rin):
  a = (3 * const.G.to_value()* M * Ar)**(1/4)
  b = (8 * np.pi * const.sigma_sb.to_value() * Rin**3)**(1/4)
  c = ( Radius / Rin )**(-3/4)
  T = ((a / b) * (1 - c))
  return T

def Tstar(M, Ar, Rin):
  a = (3 * const.G.to_value()* M * Ar)
  b = (8 * np.pi * const.sigma_sb.to_value() * Rin**3)
  return (a / b)**(1/4)

def Ferg(M, incl, T, frequency):
  return 2.4e-18 * Rg(M)**2 * (T**(8/3)) * np.cos(incl) * frequency**(1/3)

#frequency between 3900 and 4900 Angstrom
def Bband(M, Ar, incl):
  return 0.5312e44 * M**(2/3) * Ar**(2/3) * np.cos(incl)

R14 = Radius/1e14 #cm
Ms = M/(10e8*const.M_sun.to_value())
L9 = Bband(M, Ar, incl)/(10e9*const.L_sun.to_value())
Mdot26 = Ar/(10e26) #g/s

def L9(M, Ar, incl):
  return 13.8 * M**(2/3) * Ar**(2/3) * np.cos(incl)
'''




'''
#Plot of Eddington limited accretion disc
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

#ax1.set_xlim(-21, -10)
#ax1.set_ylim(0, 50)
ax1.set_ylim(-25, -5)

LumEdd = Lsum(70e6*const.M_sun.to_value(), AccR_Edd(70e6*const.M_sun.to_value()),
R_midpoints(70e6*const.M_sun.to_value()), Rin(70e6*const.M_sun.to_value()))

FluxLumEdd = Flux2(LumEdd, -0.00015677512474313146)

#print(np.sum(FluxLumEdd, axis = 0))

#LumAcRate = (Flux(Temp2(MassBH, AcRate(MassBH), R_midpoints(MassBH), Rin(MassBH)))) * Area(R_midpoints(MassBH))

#ax1.plot(np.log10((freq/10e9)), np.log10(freq * LumEdd), color='blue', linestyle='-')

ax1.plot(np.log10(freq), np.log10(freq * FluxLumEdd), color='red', linestyle='-', label='Eddington limit M81')

ax1.plot(np.log10(freq), np.log10(freq * F5M81), color='purple', linestyle='-', label='Seyfer 2 M81')
ax1.plot(np.log10(freq), np.log10(freq * FNGC6814), color='green', linestyle='-', label='Seyfert 1 NGC 6814')
ax1.plot(np.log10(freq), np.log10(freq * F3C273), color='black', linestyle='-', label='Radio-loud 3C 273')

ax1.set_ylabel(r'$\log_{10}(vF_v)$ W m$^{-2}$', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(Lambda)$ Micrometers', fontsize=16)
ax1.legend(fontsize=12, loc='upper left')

plt.show()
'''
#(const.c.to_value()/freq)*10e6



'''
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

Masses = 1e3*const.M_sun.to_value()
Masses1 = []
Mend = 1e9*const.M_sun.to_value()
AccretionR = []

#increase the mass by eddington accretion rate
def Mup(dt, M):
  return AccR_Edd(M)*dt

while Masses < Mend:
  Masses += Mup(60*60*24*365.25*1e5, Masses)
  AccretionR.append(AccR_Edd(Masses))
  Masses1.append(Masses)
  #print(Masses[0]/(const.M_sun.to_value()))#, (AccR_Edd(Masses)*60*60*24*365.25)/(const.M_sun.to_value()))

Masses1 = np.array(Masses1)
AccretionR = np.array(AccretionR)
#print(Masses1[0]/(const.M_sun.to_value()))

ax1.plot(np.log10(Masses1/(const.M_sun.to_value())), np.log10((AccretionR*(60*60*24*365.25))/(const.M_sun.to_value())), color='blue', linestyle='-')

ax1.set_ylabel(r'$\log_{10}(\dot{M})$ ($M_\odot$ yr$^{-1}$)', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(M_{\text{AGN}})$ M$_\odot$', fontsize=16)
#plt.show()




#Time taken to reach supermassive black hole (10^9) 
def TimeTaken(M, Mf):
  t = 0
  dt = 60*60*24*365.25*1e5
  for i in range(len(M)):
    while M[i] < Mf:
      M[i] += AccR_Edd(M[i])*dt
      t += dt
  return t



Massss = np.linspace(0.001 * const.M_sun.to_value(), 1e5 * const.M_sun.to_value(), 1000)
time_taken = [TimeTaken([m], 1e9 * const.M_sun.to_value()) for m in Massss]

fig, ax2 = plt.subplots()
fig.set_size_inches(10, 6)
ax2.plot(np.log10(Massss/const.M_sun.to_value()), np.log10(np.array(time_taken)/(60*60*24*365.25)), color='blue', linestyle='-')

ax2.set_xlabel(r'$\log_{10}(M_{\text{AGN}})$ M$_\odot$', fontsize=16)
ax2.set_ylabel(r'$\log_{10}(\text{Time Taken})$ (years)', fontsize=16)

plt.show()
'''


















#Milestone values
LumMilestone = (Flux(Temp2(MassBH, AccR, R_midpoints(MassBH), Rin(MassBH)))) * Area(R_midpoints(MassBH))  
LumMilestonesum = np.sum(LumMilestone, axis = 0)

'''
#milestone values, max spin black hole
L3 = (Flux(Temp2(MassBH, AccR, R_midpoints2, Rin2))) * Area(R_midpoints2) 
L3sum = np.sum(L3, axis = 0)



#Object mass = mass of sun
LumSunMass = (Flux(Temp2(Mass2, AccR, R_midpoints(Mass2), Rin(Mass2)))) * Area(R_midpoints(Mass2))
LumSunMass2um = np.sum(LumSunMass, axis = 0)



#mass of sun with spin
LumSunMass2pin = (Flux(Temp2(Mass2, AccR, R_midpoints2, Rin2))) * Area(R_midpoints2)
LumSunMass2umSpin = np.sum(LumSunMass2pin, axis = 0)


#Different accretion rate
LumdAccr = (Flux(Temp2(MassBH, AccR2, R_midpoints, Rin))) * Area(R_midpoints)
LumdAccrSum = np.sum(LumdAccr, axis = 0)



#D Accr with spin
LumdAccrSpin = (Flux(Temp2(MassBH, AccR2, R_midpoints2, Rin2))) * Area(R_midpoints2)
LumdAccrSumSpin = np.sum(LumdAccrSpin, axis = 0)
'''

'''
#Luminosity eddington limit
def L_edd(Mass):
  return (4 * np.pi * const.G.value * Mass * const.m_p.value * const.c.value)/const.sigma_T.value

eta = 0.1  #standard efficiency with which mass is assumed to produce energy near the event horizon

#max accretion rate
def AccR_Edd(Mass):
  return L_edd(Mass)/(eta * const.c.value**2)


#print(f'Eddington limit accretion = {AccR_Edd}')


LumEdd = (Flux(Temp2(MassBH, AccR_Edd(MassBH), R_midpoints(MassBH), Rin(MassBH)))) * Area(R_midpoints(MassBH))  
LumEddsum = np.sum(LumEdd, axis = 0)
'''



'''
#Using scaled radius R/Rg and equation
Lscaled = (Flux(Temp(MassBH, AccR, r_midpoints, rin))) * Area(r_midpoints)  
LscaledSum = np.sum(Lscaled, axis = 0)
'''

'''
L_non_viscous = Flux((Temp3(MassBH, AccR, R_midpoints))) * Area(R_midpoints) 
Lsum_non_viscous = np.sum(L_non_viscous, axis=0)
'''


LumSMBH1 = (Flux(Temp2(MassSMBH1, AccR, R_midpoints(MassSMBH1), Rin(MassSMBH1)))) * Area(R_midpoints(MassSMBH1))  
LumSMBH1sum = np.sum(LumSMBH1, axis = 0)

LumSMBH2 = (Flux(Temp2(MassSMBH2, AccR, R_midpoints(MassSMBH2), Rin(MassSMBH2)))) * Area(R_midpoints(MassSMBH2))  
LumSMBH2sum = np.sum(LumSMBH2, axis = 0)

LumSMBH3 = (Flux(Temp2(MassSMBH3, AccR, R_midpoints(MassSMBH3), Rin(MassSMBH3)))) * Area(R_midpoints(MassSMBH3))  
LumSMBH3sum = np.sum(LumSMBH3, axis = 0)



TotLumMilestone = scipy.integrate.trapezoid(LumMilestonesum, freq)*u.W  #.to(u.W)
#TotLscaled = scipy.integrate.trapezoid(LscaledSum, freq)*u.W  #.to(u.W)
#TotL3 = scipy.integrate.trapezoid(L3sum, freq)*u.W  #.to(u.W)
#TotLumSunMass = scipy.integrate.trapezoid(LumSunMass2um, freq)*u.W  #.to(u.W)
#TotLumSunMass2pin = scipy.integrate.trapezoid(LumSunMass2umSpin, freq)*u.W  #.to(u.W)
#TotLumdAccr = scipy.integrate.trapezoid(LumdAccrSum, freq)*u.W  #.to(u.W)
#TotLumdAccrSpin = scipy.integrate.trapezoid(LumdAccrSumSpin, freq)*u.W  #.to(u.W)
#TotLumMilestone2 = scipy.integrate.trapezoid(Lsum_non_viscous, freq)*u.W  #.to(u.W)
TotLumSMBH1 = scipy.integrate.trapezoid(LumSMBH1sum, freq)*u.W
TotLumSMBH2 = scipy.integrate.trapezoid(LumSMBH2sum, freq)*u.W
TotLumSMBH3 = scipy.integrate.trapezoid(LumSMBH3sum, freq)*u.W
#TotLumEdd = scipy.integrate.trapezoid(LumEddsum, freq)*u.W





#print(f'(2) L sum = {TotLumSunMass}')
#print(f'(Milestone) L sum = {TotLumMilestone}')
#print(f'(1) L sum = {TotLumSMBH1}')
#print(f'(2) L sum = {TotLumSMBH2}')
#print(f'(3) L sum = {TotLumSMBH3}')
#print(f'(Eddington) L sum = {TotLumEdd}')



'''
MassAGN = np.linspace(10**6, 10**10, 1000) * const.M_sun.to_value()
LumAGN = (Flux(Temp2(MassAGN, AccR_Edd, R_midpoints, Rin))) * Area(R_midpoints) 
LumAGNsum = np.sum(LumAGN, axis = 0)
TotLumAGN = scipy.integrate.trapezoid(LumAGNsum, freq)


#flux received
def f(Lum, dist):
  return Lum/(4*np.pi * dist**2)

#flux to magnitude with no filter with respect to a zero point source
def mag(f):
  return m_zero_point - 2.5 * np.log10(f)
'''





'''
print(f'(Milestone) L sum = {TotLumMilestone}')
print(f'(Milestone) L sum = {TotLumSunMass}')
print(f'(1) L sum = {TotLumMilestone2}')
print(TotLumMilestone2/TotLumMilestone)
print(TotLumSunMass/TotLumMilestone)
print(TotL3/TotLumMilestone)
print(TotLumSunMass2pin/TotLumSunMass)
'''



#IMPORTAMNT PLOT
'''
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

#fig.patch.set_facecolor('#D5D5D5') 
#ax1.set_facecolor('#D5D5D5') 

ax1.set_xlim(7, 20)
ax1.set_ylim(15, 35)

d, = ax1.plot(np.log10(freq), np.log10(freq * LumMilestonesum), linestyle='-', color='black', label='10')
#e, = ax1.plot(np.log10(freq), np.log10(freq * LumEddsum), linestyle='-', color='purple', label='Eddington Limit (10)')
a, = ax1.plot(np.log10(freq), np.log10(freq * LumSMBH1sum), linestyle='-', color='green', label='10^7')
b, = ax1.plot(np.log10(freq), np.log10(freq * LumSMBH2sum), linestyle='-', color='blue', label='10^9')
c, = ax1.plot(np.log10(freq), np.log10(freq * LumSMBH3sum), linestyle='-', color='red', label='10^11')
#f, = ax1.plot(np.log10(freq), np.log10(freq * LumSunMass2um), linestyle='-', color='brown', label='10^4')
#g, = ax1.plot(np.log10(freq), np.log10(Lsum4(MassBH, AccR, r(MassBH), rin(MassBH))), linestyle='-', color='orange', label='scaled')


#dAccr, = ax1.plot(np.log10(freq), np.log10(freq * LumSunMass2um), linestyle='-', color='green')
#dAccrSpin, = ax1.plot(np.log10(freq), np.log10(freq * LumSunMass2umSpin), linestyle=':', color='green')


no_spin_line = plt.Line2D([0], [0], color='black', linestyle='-', label='No Spin')
max_spin_line = plt.Line2D([0], [0], color='black', linestyle=':', label='Max Spin')
max_spin_line2 = plt.Line2D([0], [0], color='red', linestyle='--', label='SMBH')
mass_10msun = plt.Line2D([0], [0], color='blue', linestyle='-', label=r'Mass=10$M_{\odot}$')
mass_1msun = plt.Line2D([0], [0], color='green', linestyle='-', label=r'Mass=10$^{4}$$M_{\odot}$')
ax1.legend(handles=[no_spin_line, max_spin_line, mass_10msun, mass_1msun, max_spin_line2], fontsize=12)



ax1.set_ylabel(r'$\log_{10}(vL_v)$ W', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(v)$ Hz', fontsize=16)
ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)
ax1.legend(fontsize=12)

# Create a secondary x-axis
ax2 = ax1.twiny()

freq_ticks = ax1.get_xticks()

# Convert frequency ticks to wavelength ticks
wavelength_ticks = np.log10(const.c.value / (10**freq_ticks))

# Create custom labels
custom_labels = [f'{tick:.0f}' if tick != 0 else '-e' for tick in wavelength_ticks]

# Set the secondary x-axis limits and labels
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(freq_ticks)
ax2.set_xticklabels(custom_labels)
ax2.set_xlabel(r'$\log_{10}(\lambda/3)$ m', fontsize=16)
ax2.tick_params(axis='x', labelsize=14)

peak_indices = {
    'MileStone': np.argmax(LumSMBH1sum),
    'MaxSpin': np.argmax(LumSMBH2sum),
    'NonVisous': np.argmax(LumSMBH3sum),
    
}

peak_frequencies = {key: freq[idx] for key, idx in peak_indices.items()}
peak_y_values = {key: np.log10(freq[peak_indices[key]] * value[peak_indices[key]]) for key, value in zip(peak_indices.keys(),
[LumSMBH1sum, LumSMBH2sum, LumSMBH3sum])}


print("Peak frequencies (Hz) and corresponding Luminosities (W):")
for key in peak_frequencies:
    print(f"{key}: Frequency = {peak_frequencies[key]:.2e} Hz, Luminosity = {10**peak_y_values[key]:.2e}")


plt.show()
'''




#Log(Lum/Lum_edd) vs log(AccR) as seen in 'slim accretion discs' by Abramowicz and buddies
#doesnt look same, linear even past eddington limit, wrote by chat so could be wrong, could be using different equations to the paper (prob)
#should become flat past eddington limit

'''
accretion_rate_ratios = np.logspace(np.log10(0.001), np.log10(10), 10)  # From 0.001 to 10 times the Eddington accretion rate
accretion_rates = accretion_rate_ratios * AccR_Edd

#print(accretion_rates)

luminosities = np.zeros(len(accretion_rate_ratios))
luminosities_NV = np.zeros(len(accretion_rate_ratios))


# Calculate the luminosity for each accretion rate
for i, acc_rate in enumerate(accretion_rates):
    Lum = Flux(Temp2(MassBH, acc_rate, R_midpoints, Rin)) * Area(R_midpoints)
    luminosities[i] = scipy.integrate.trapezoid(np.sum(Lum, axis=0), freq)

for i, acc_rate in enumerate(accretion_rates):
    Lum2 = Flux(Temp3(MassBH, acc_rate, R_midpoints)) * Area(R_midpoints)
    luminosities_NV[i] = scipy.integrate.trapezoid(np.sum(Lum2, axis=0), freq)
  
L_Le = luminosities*u.W / TotLumEdd
L_Le_NV = luminosities_NV*u.W / TotLumEdd


# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(np.log10(accretion_rate_ratios), np.log10(L_Le), linestyle='-', color='b', label='Viscous')
plt.plot(np.log10(accretion_rate_ratios), np.log10(L_Le_NV), linestyle='-', color='r', label='Non-Viscous')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Log(Accretion Rate / Eddington Accretion Rate)')
plt.ylabel('Log(Luminosity / Eddington Luminosity)')
plt.legend()
plt.show()
'''



#Spectrums of black holes with different accretion rates
'''
accretion_rates = [0.1, 1, 10]  # Different accretion rates (in units of Eddington accretion rate)

LLowEdd = (Flux(Temp2(MassSMBH2, 0.1*AccR_Edd(MassSMBH2), R_midpoints(MassSMBH2), Rin(MassSMBH2)))) * Area(R_midpoints(MassSMBH2))  
LLowEddSum = np.sum(LLowEdd, axis = 0)
LAboveEdd = (Flux(Temp2(MassSMBH2, 3*AccR_Edd(MassSMBH2), R_midpoints(MassSMBH2), Rin(MassSMBH2)))) * Area(R_midpoints(MassSMBH2))  
LAEddSum = np.sum(LAboveEdd, axis = 0)
L1Edd = (Flux(Temp2(MassSMBH2, AccR_Edd(MassSMBH2), R_midpoints(MassSMBH2), Rin(MassSMBH2)))) * Area(R_midpoints(MassSMBH2))  
L1EddSum = np.sum(L1Edd, axis = 0)
LFarLessEdd = (Flux(Temp2(MassBH, 0.0001*AccR_Edd(MassSMBH2), R_midpoints(MassSMBH2), Rin(MassSMBH2)))) * Area(R_midpoints(MassSMBH2))  
LFLEddSum = np.sum(LFarLessEdd, axis = 0)
#LumEdd2 = (Flux(Temp2(30*MassBH, AccR_Edd, R_midpoints(MassSMBH2), Rin))) * Area(R_midpoints(MassSMBH2))  
#LumEddsum2 = np.sum(LumEdd2, axis = 0)
#LumEdd3 = (Flux(Temp2(60*MassBH, AccR_Edd, R_midpoints(MassSMBH2), Rin))) * Area(R_midpoints(MassSMBH2))  
#LumEddsum3 = np.sum(LumEdd3, axis = 0)

LumM81 = (Flux(Temp2(70e6 * const.M_sun.to_value(), 0.001*AccR_Edd(70e6 * const.M_sun.to_value()), R_midpoints(70e6 * const.M_sun.to_value()), 
Rin(70e6 * const.M_sun.to_value())))) * Area(R_midpoints(70e6 * const.M_sun.to_value()))  
LumM81Sum = np.sum(LumM81, axis = 0)

fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

#ax1.set_xlim(9, 20)
#ax1.set_ylim(15, 45)

#lowEdd, = ax1.plot(np.log10(freq), np.log10(freq * LLowEddSum), color='blue', linestyle='-', label='Accretion rate = 0.1 Edd')
#Edd1, = ax1.plot(np.log10(freq), np.log10(freq * LumEddsum), color='red', linestyle='-', label='Accretion rate = 1 Edd')
#Edd11, = ax1.plot(np.log10(freq), np.log10(freq * L1EddSum), color='black', linestyle='-', label='Accretion rate = 1 Edd')
#Edd2, = ax1.plot(np.log10(freq), np.log10(freq * LumEddsum2), color='red', linestyle=':')
#Edd3, = ax1.plot(np.log10(freq), np.log10(freq * LumEddsum3), color='red', linestyle='--')
#AboveEdd, = ax1.plot(np.log10(freq), np.log10(freq * LAEddSum), color='purple', linestyle='-', label='Accretion rate = 3 Edd')
#FarLessThanEdd, = ax1.plot(np.log10(freq), np.log10(freq * LFLEddSum), color='green', linestyle='-', label='Accretion rate = 0.001 Edd')

ax1.plot(np.log10(freq), np.log10(freq * LumM81Sum), color='blue', linestyle='-', label='M81')

mass1_line = plt.Line2D([0], [0], color='black', linestyle='-', label='Mass = MassBH')
mass30_line = plt.Line2D([0], [0], color='black', linestyle=':', label='Mass = 30*MassBH')
mass60_line = plt.Line2D([0], [0], color='black', linestyle='--', label='Mass = 60*MassBH')

ax1.legend(fontsize=12)


ax1.set_ylabel(r'$\log_{10}(vL_v)$ W', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(v)$ Hz', fontsize=16)
ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)

# Create a secondary x-axis
ax2 = ax1.twiny()

freq_ticks = ax1.get_xticks()

# Convert frequency ticks to wavelength ticks
wavelength_ticks = np.log10(const.c.value / (10**freq_ticks))

# Create custom labels
custom_labels = [f'{tick:.0f}' if tick != 0 else '-e' for tick in wavelength_ticks]

# Set the secondary x-axis limits and labels
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(freq_ticks)
ax2.set_xticklabels(custom_labels)
ax2.set_xlabel(r'$\log_{10}(\lambda/3)$ m', fontsize=16)
ax2.tick_params(axis='x', labelsize=14)

plt.show()
'''







'''
print(f'(Milestone) L sum = {TotLumMilestone}')
print(f'(Equation+scaled r) L sum = {TotLscaled}')
print(f'(Milestone + spin) L sum = {TotL3}')
print(f'(Sun mass) L sum = {TotLumSunMass}')
print(f'(Sun mass + spin) L sum = {TotLumSunMass2pin}')
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







#plot each blackbodies on same plot
#make spectrums to location in ring, draw blackhole?



#Plot used in poster, luminosity*freq vs frequency
'''
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

fig.patch.set_facecolor('#D5D5D5') 
ax1.set_facecolor('#D5D5D5') 

ax1.set_xlim(11, 20)
ax1.set_ylim(15, 38)

MileStone, = ax1.plot(np.log10(freq), np.log10(freq * LumMilestonesum), linestyle='-', color='blue')
MaxSpin, = ax1.plot(np.log10(freq), np.log10(freq * L3sum), linestyle=':', color='blue')
NonVisous, = ax1.plot(np.log10(freq), np.log10(freq * Lsum_non_viscous), linestyle='--', color='red', label='Non-Viscous')
dAccr, = ax1.plot(np.log10(freq), np.log10(freq * LumSunMass2um), linestyle='-', color='green')
dAccrSpin, = ax1.plot(np.log10(freq), np.log10(freq * LumSunMass2umSpin), linestyle=':', color='green')

no_spin_line = plt.Line2D([0], [0], color='black', linestyle='-', label='No Spin')
max_spin_line = plt.Line2D([0], [0], color='black', linestyle=':', label='Max Spin')
max_spin_line2 = plt.Line2D([0], [0], color='red', linestyle='--', label='Non-Viscous')
mass_10msun = plt.Line2D([0], [0], color='blue', linestyle='-', label=r'Mass=10$M_{\odot}$')
mass_1msun = plt.Line2D([0], [0], color='green', linestyle='-', label=r'Mass=10$^{4}$$M_{\odot}$')
ax1.legend(handles=[no_spin_line, max_spin_line, mass_10msun, mass_1msun, max_spin_line2], fontsize=12)

ax1.set_ylabel(r'$\log_{10}(vL_v)$ W', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(v)$ Hz', fontsize=16)
ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)

# Create a secondary x-axis
ax2 = ax1.twiny()

freq_ticks = ax1.get_xticks()

# Convert frequency ticks to wavelength ticks
wavelength_ticks = np.log10(const.c.value / (10**freq_ticks))

# Create custom labels
custom_labels = [f'{tick:.0f}' if tick != 0 else '-e' for tick in wavelength_ticks]

# Set the secondary x-axis limits and labels
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(freq_ticks)
ax2.set_xticklabels(custom_labels)
ax2.set_xlabel(r'$\log_{10}(\lambda/3)$ m', fontsize=16)
ax2.tick_params(axis='x', labelsize=14)

peak_indices = {
    'MileStone': np.argmax(LumMilestonesum),
    'MaxSpin': np.argmax(L3sum),
    'NonVisous': np.argmax(Lsum_non_viscous),
    'dAccr': np.argmax(LumSunMass2um),
    'dAccrSpin': np.argmax(LumSunMass2umSpin)
}

peak_frequencies = {key: freq[idx] for key, idx in peak_indices.items()}
peak_y_values = {key: np.log10(freq[peak_indices[key]] * value[peak_indices[key]]) for key, value in zip(peak_indices.keys(), [LumMilestonesum, L3sum, Lsum_non_viscous, LumSunMass2um, LumSunMass2umSpin])}

print("Peak frequencies (Hz) and corresponding Luminosities (W):")
for key in peak_frequencies:
    print(f"{key}: Frequency = {peak_frequencies[key]:.2e} Hz, Y-value = {10**peak_y_values[key]:.2e}")


plt.show()
'''


#DO HOW LONG IT TAKES TO GET A CONSTANT ANSWER FOR LUMINOSITY
#LUMINOSITY ANSWER/LUMINOSITY FOR SAID NUMBER OF RINGS MAYBES?


#same as last one, but doesnt fully work
'''
plt.figure(figsize=(10, 6))
plt.xlim(11, 19)
plt.ylim(15, 32)
plt.plot(np.log10(freq), np.log10(freq * LumMilestonesum),  linestyle='-', color='blue', label='Viscous')
plt.plot(np.log10(freq), np.log10(freq * Lsum_non_viscous), linestyle='-', color='red', label='Non-Viscous')

# Add labels and legend
plt.xlabel(r'$\log_{10}(v)$ Hz', fontsize=16)
plt.ylabel(r'$\log_{10}(vL_v)$ W', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.show()
'''

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



#AND USE THIS ONE
#shows how individual temp spectrums make up the larger spectrum for the whole accreion disc
'''
Te = Temp2(MassBH, AccR, R_midpoints, Rin)
areas = Area(R_midpoints)



fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

fig.patch.set_facecolor('#D5D5D5')  # Adjust this hex color for accuracy
ax1.set_facecolor('#D5D5D5') 

ax1.set_xlim(11, 19)
ax1.set_ylim(15, 32)
index2 = [100, 2000, 3500, 5000, 6500, 8000, 9500]
colors = ['violet', 'indigo', 'blue', 'green', 'yellow', 'orange', 'red']
for j, color in zip(index2, colors):
    L = Flux(Te[j]) * areas[j]
    Lsum = np.sum(L, axis=0)
    temp_label = f'{Te[j]:.2g}'
    temp_label = temp_label.replace('e+0', r'\times 10^{').replace('e-0', r'\times 10^{-') + '}'
    ax1.plot(np.log10(freq), np.log10(freq * Lsum), linestyle='--', color=color, label=f'${temp_label}$ K')
ax1.plot(np.log10(freq), np.log10(freq * LumMilestonesum), linestyle='-', color='black')
ax1.legend(fontsize=12)
ax1.set_ylabel(r'$\log_{10}(vL_v)$ W', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(v)$ Hz', fontsize=16)
ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)

ax2 = ax1.twiny()

freq_ticks = ax1.get_xticks()

# Convert frequency ticks to wavelength ticks
wavelength_ticks = np.log10(const.c.value / (10**freq_ticks))

# Create custom labels
custom_labels = [f'{tick:.0f}' if tick != 0 else '-e' for tick in wavelength_ticks]

# Set the secondary x-axis limits and labels
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(freq_ticks)
ax2.set_xticklabels(custom_labels)
ax2.set_xlabel(r'$\log_{10}(\lambda/3)$ m', fontsize=16)
ax2.tick_params(axis='x', labelsize=14)


plt.show()
'''



# Plot the luminosity spectrum for each accretion rate but doesnt work, Max_AccR not defined
'''
accretion_rates1 = [0.1 * Max_AccR, 0.5 * Max_AccR, Max_AccR, 2 * Max_AccR]
accretion_rates2 = np.logspace(np.log10(0.01 * Max_AccR), np.log10(2000 * Max_AccR), 5)

# Define a list of colors
colors = ['red', 'orange', 'green', 'blue']

# Plot the luminosity spectrum for each accretion rate
plt.figure(figsize=(10, 6))
plt.xlim(13, 19)
plt.ylim(20, 35)
for Ar, color in zip(accretion_rates1, colors):
    Te = Temp2(MassBH, Ar, R_midpoints, Rin)
    areas = Area(R_midpoints)
    L = Flux(Te) * areas
    Lsum = np.sum(L, axis=0)
    acc_label = f'{Ar:.2e}'
    acc_label = acc_label.replace('e+0', r'\times 10^{').replace('e', r'\times 10^{').replace('+', '') + '}'
    plt.plot(np.log10(freq), np.log10(freq * Lsum), linestyle='--', color=color, label=f'Accretion Rate: ${acc_label}$ kg/s')

# Assuming LumMilestonesum is defined elsewhere in your code
# Plot the black line for comparison
plt.plot(np.log10(freq), np.log10(freq * LumMilestonesum), linestyle='-', color='black', label='Milestone')

# Add labels and legend
plt.xlabel(r'$\log_{10}(f)$ Hz', fontsize=16)
plt.ylabel(r'$\log_{10}(fL_f)$ W', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.show()
'''

#spectrums for a range of masses
'''
masses = np.logspace(np.log10(const.M_sun.to_value()), np.log10(100 * const.M_sun.to_value()), 5)

# Define a list of colors
colors = ['red', 'orange', 'yellow', 'green', 'blue']

# Plot the luminosity spectrum for each Eddington-limited accretion rate
plt.figure(figsize=(10, 6))
plt.xlim(11, 19)
plt.ylim(15, 37)
for M, color in zip(masses, colors):
    Max_AccR = L_edd(M) / (0.4 * const.c.value**2)  # Eddington-limited accretion rate
    Te = Temp2(M, Max_AccR, R_midpoints, Rin)
    areas = Area(R_midpoints)
    L = Flux(Te) * areas
    Lsum = np.sum(L, axis=0)
    mass_label = f'{M:.2e}'
    mass_label = mass_label.replace('e+0', r'\times 10^{').replace('e', r'\times 10^{').replace('+', '') + '}'
    plt.plot(np.log10(freq), np.log10(freq * Lsum), linestyle='--', color=color, label=f'Mass:${mass_label}$ kg')

# Assuming LumMilestonesum is defined elsewhere in your code
# Plot the black line for comparison
plt.plot(np.log10(freq), np.log10(freq * LumMilestonesum), linestyle='-', color='black', label='Milestone')

# Add labels and legend
plt.xlabel(r'$\log_{10}(f)$ Hz', fontsize=16)
plt.ylabel(r'$\log_{10}(fL_f)$ W', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.show()
'''


#kinda the same as last one
'''
masses = np.logspace(np.log10(10 * const.M_sun.to_value()), np.log10(1000900 * const.M_sun.to_value()), 5)

# Define a list of colors
colors = ['red', 'orange', 'yellow', 'green', 'blue']

# Plot the luminosity spectrum for each black hole mass
plt.figure(figsize=(10, 6))
plt.xlim(11, 25)
plt.ylim(15, 40)
for M, color in zip(masses, colors):
    Te = Temp2(M, AccR, R_midpoints, Rin)
    areas = Area(R_midpoints)
    L = Flux(Te) * areas
    Lsum = np.sum(L, axis=0)
    mass_label = f'{M:.1e}'
    mass_label = mass_label.replace('e+0', r'\times 10^{').replace('e+0', r'\times 10^{').replace('e', r'\times 10^{').replace('+', '') + '}'
    plt.plot(np.log10(freq), np.log10(freq * Lsum), linestyle='--', color=color, label=f'Mass: ${mass_label}$ kg')

# Assuming LumMilestonesum is defined elsewhere in your code
# Plot the black line for comparison
plt.plot(np.log10(freq), np.log10(freq * LumMilestonesum), linestyle='-', color='black', label='Milestone')

# Add labels and legend
plt.xlabel(r'$\log_{10}(f)$ Hz', fontsize=16)
plt.ylabel(r'$\log_{10}(fL_f)$ W', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.show()
'''


#same plot as poster one, but includes super massive black hole
'''
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

fig.patch.set_facecolor('#D5D5D5') 
ax1.set_facecolor('#D5D5D5') 

ax1.set_xlim(11, 22)
ax1.set_ylim(15, 40)

MileStone, = ax1.plot(np.log10(freq), np.log10(freq * LumMilestonesum), linestyle='-', color='blue')
MaxSpin, = ax1.plot(np.log10(freq), np.log10(freq * L3sum), linestyle=':', color='blue')
SMBH, = ax1.plot(np.log10(freq), np.log10(freq * LumSMBH1sum), linestyle='--', color='red', label='SMBH')
#dAccr, = ax1.plot(np.log10(freq), np.log10(freq * LumSunMass2um), linestyle='-', color='green')
#dAccrSpin, = ax1.plot(np.log10(freq), np.log10(freq * LumSunMass2umSpin), linestyle=':', color='green')

no_spin_line = plt.Line2D([0], [0], color='black', linestyle='-', label='No Spin')
max_spin_line = plt.Line2D([0], [0], color='black', linestyle=':', label='Max Spin')
max_spin_line2 = plt.Line2D([0], [0], color='red', linestyle='--', label='SMBH')
mass_10msun = plt.Line2D([0], [0], color='blue', linestyle='-', label=r'Mass=10$M_{\odot}$')
mass_1msun = plt.Line2D([0], [0], color='green', linestyle='-', label=r'Mass=10$^{4}$$M_{\odot}$')
ax1.legend(handles=[no_spin_line, max_spin_line, mass_10msun, mass_1msun, max_spin_line2], fontsize=12)

ax1.set_ylabel(r'$\log_{10}(vL_v)$ W', fontsize=16)
ax1.set_xlabel(r'$\log_{10}(v)$ Hz', fontsize=16)
ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)

# Create a secondary x-axis
ax2 = ax1.twiny()

freq_ticks = ax1.get_xticks()

# Convert frequency ticks to wavelength ticks
wavelength_ticks = np.log10(const.c.value / (10**freq_ticks))

# Create custom labels
custom_labels = [f'{tick:.0f}' if tick != 0 else '-e' for tick in wavelength_ticks]

# Set the secondary x-axis limits and labels
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(freq_ticks)
ax2.set_xticklabels(custom_labels)
ax2.set_xlabel(r'$\log_{10}(\lambda/3)$ m', fontsize=16)
ax2.tick_params(axis='x', labelsize=14)

peak_indices = {
    'MileStone': np.argmax(LumMilestonesum),
    'MaxSpin': np.argmax(L3sum),
    'NonVisous': np.argmax(Lsum_non_viscous),
    'dAccr': np.argmax(LumSunMass2um),
    'dAccrSpin': np.argmax(LumSunMass2umSpin)
}

peak_frequencies = {key: freq[idx] for key, idx in peak_indices.items()}
peak_y_values = {key: np.log10(freq[peak_indices[key]] * value[peak_indices[key]]) for key, value in zip(peak_indices.keys(), [LumMilestonesum, L3sum, Lsum_non_viscous, LumSunMass2um, LumSunMass2umSpin])}

print("Peak frequencies (Hz) and corresponding Luminosities (W):")
for key in peak_frequencies:
    print(f"{key}: Frequency = {peak_frequencies[key]:.2e} Hz, Luminosity = {10**peak_y_values[key]:.2e}")


plt.show()
'''
'''
#Magnitude to flux equation (flux in watts per square kilometer)
def flux(mag):
    return 10**((-mag)/2.5) * flux_ref_B
    #.to(u.J / (u.km**2 * u.s * u.Hz)) #change this depending on filter used



#Luminosity of AGN from the apparent magnitude
def Lum(mag, dist):
    return flux(mag) * (4*np.pi * dist**2)


#flux_ref = #changes depending on filter used
#values given in erg/(s*cm^2*Hz)
flux_ref_U = 1.81e-20
flux_ref_B = 4.26e-20 
#* u.erg / (u.s * u.cm**2 * u.Hz)
flux_ref_V = 3.64e-20
flux_ref_R = 3.08e-20
flux_ref_I = 1.554e-20


D = np.array([
    113, 672, 740, 209, 115, 145, 142, 96.7, 461, 287,
    155, 1170, 368, 758, 283, 740, 133, 410, 395, 138,
    606, 522, 1510, 251, 151, 283
])

mag = np.array([
    16.16, 16.21, 17.81, 15.13, 15.59, 16.64, 14.50, 15.84, 16.70, 16.75,
    18.08, 17.64, 17.16, 15.49, 17.17, 16.66, 16.13, 16.70, 15.36, 17.62,
    16.14, 16.36, 17.54, 16.74, 13.96, 18.73
])


#List of Black hole mass and luminoisty
L_bol1 = np.log10((10**np.array([
    45.34, 44.88, 44.91, 45.23, 44.78, 44.57, 44.71, 44.69, 45.03, 44.63,
    44.99, 43.86, 44.29, 44.41, 43.56, 43.73, 44.09, 44.83, 45.28, 45.39,
    45.93, 45.93, 45.36, 46.16, 45.81, 45.01, 45.83, 45.50, 45.58, 45.19,
    45.66]) * 10**(-7))/const.L_sun.to_value())

log_L2 = np.log10(10**np.array([
    9.53, 11.09, 10.50, 10.49, 9.78, 9.57, 10.40, 9.52, 
    10.56, 10.13, 9.05, 10.94, 10.18, 11.45, 9.94, 10.99, 9.69,
    10.46, 10.99, 9.13, 11.02, 10.83, 11.22, 10.01, 10.71, 9.32]))

L_bol = np.concatenate((L_bol1, log_L2), axis = None)

M_BH1 = np.array([
    7.42, 8.55, 8.27, 7.91, 6.77, 7.86, 6.82, 6.69, 7.86, 7.20,
    7.60, 7.64, 7.36, 6.94, 6.13, 7.13, 6.91, 8.03, 6.84, 7.58,
    8.41, 8.24, 7.38, 8.24, 7.49, 8.56, 7.90, 8.48, 7.57, 7.92,
    8.62]) 
    
M_BH2 = np.array([7.15, 8.59, 8.57, 8.41, 7.68, 7.74, 8.18, 7.72, 
    8.84, 7.97, 7.40, 8.44, 8.16, 8.95, 7.86, 8.64,
    7.54, 8.65, 9.11, 7.69, 8.45, 8.77, 8.89, 8.46, 8.16, 7.58])

M_BH = np.concatenate((M_BH1, M_BH2), axis = None)






#Spatially Resolved Kinematics
L1 = np.log10((10**np.array([44.98, 43.45]) / 3.826E+33))
M1 = np.array([7.23, 7.62]) 
T1 = np.array(["SY2", "SY2"])
Z1 = np.array([0.004, 0.001])

#Reverberaiton Mapping
L2 = np.log10((10**np.array([45.34, 44.88, 44.91, 45.23, 44.78, 44.57, 44.71, 44.69, 45.03, 44.63, 
        44.99, 43.86, 44.29, 44.41, 43.56, 43.73, 44.09, 44.83, 45.28, 45.39, 
        45.93, 45.36, 46.16, 45.81, 45.01, 45.83, 45.50, 45.58, 45.19, 45.66,
        45.52, 46.56, 45.47, 47.35, 46.33]) / 3.826E+33))

M2 = np.array([7.42, 8.55, 8.27, 7.91, 6.77, 7.86, 6.82, 6.69, 7.86, 7.20, 
       7.60, 7.64, 7.36, 6.94, 6.13, 7.13, 6.91, 8.03, 6.84, 7.58, 
       8.41, 7.38, 8.24, 7.49, 8.56, 7.90, 8.48, 7.57, 7.92, 8.62,
       7.88, 8.31, 7.74, 7.22, 8.23])

T2 = np.array(["SY1", "SY1", "SY1", "SY1", "SY1", "SY1", "SY1", "SY1", "SY1", "SY1",
        "SY1", "SY1", "SY1", "SY1", "SY1", "SY1", "SY1", "SY1", "SY1", "RQQ",
        "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", 
        "RQQ", "RQQ", "RQQ", "RLQ", "RLQ"])

Z2 = np.array([0.033, 0.056, 0.032, 0.047, 0.016, 0.022, 0.035, 0.026, 0.034, 0.026, 
     0.032, 0.004, 0.009, 0.010, 0.002, 0.003, 0.009, 0.017, 0.016, 0.142, 
     0.155, 0.100, 0.064, 0.239, 0.085, 0.064, 0.155, 0.087, 0.089, 0.086, 
     0.129, 0.114, 0.292, 0.061, 0.158, 0.371])


#Optical Luminoisty
L3 = np.log10((10**np.array([45.84, 44.40, 43.92, 45.47, 45.62, 45.05, 45.45, 45.51, 46.22, 45.51, 
         46.03, 46.02, 45.39, 45.63, 45.13, 45.93, 45.13, 44.98, 45.31, 46.54, 
         45.78, 41.02, 45.57, 45.83, 46.69, 46.44, 46.64, 45.22, 44.74, 46.84, 
         44.92, 44.94, 47.72, 46.01, 46.32, 46.47, 47.40, 47.00, 46.15, 46.12,
         45.32, 45.12, 46.36, 45.36, 46.34, 46.43, 45.69, 46.30, 47.16, 45.58,
         46.41, 45.97, 46.94, 46.54, 46.23, 45.99, 45.26, 44.63, 46.26, 46.35,
         46.59, 46.54, 46.21, 46.71, 46.63, 44.87, 46.2, 44.55, 45.8, 46.49, 
         46.33, 44.97, 44.25, 46.26, 44.08, 45.98, 45.56, 46.41, 45.81, 45.83, 
         46.63, 46.48, 45.61, 46.1, 45.52, 45.86, 45.81, 47.11, 46.48, 46.19, 
         43.94, 46.16, 46.93, 44.54, 46.38, 45.86, 46.0, 44.94, 46.99, 45.47, 
         46.68, 46.89, 45.78, 47.21, 44.01, 45.63, 46.07, 45.85, 46.78, 46.23, 
         46.21, 46.68, 45.54, 46.31, 45.32, 46.23, 46.84, 45.83, 45.75, 46.76, 
         46.17, 46.23, 46.65, 45.67, 46.62, 47.17, 46.11, 45.47, 47.27, 46.96,
         46.55, 46.22, 45.56, 47.07, 45.92, 45.94, 45.01]) / 3.826E+33))

M3 = np.array([8.10, 6.54, 7.28, 8.90, 7.70, 6.67, 7.86, 8.03, 8.94, 7.79, 
        9.08, 8.21, 8.29, 8.00, 7.29, 8.06, 8.10, 7.91, 6.37, 8.71, 
        8.34, 6.73, 8.38, 9.52, 8.73, 8.74, 9.13, 8.57, 7.56, 8.97, 
        6.54, 7.29, 8.52, 8.60, 8.98, 9.07, 9.47, 9.03, 8.79, 8.53, 
        8.13, 7.42, 8.88, 7.55, 8.53, 9.58, 9.02, 8.68, 9.41, 8.74, 
        8.67, 8.00, 9.40, 7.96, 8.52, 7.90, 7.72, 8.14, 9.28, 8.46,
        9.00, 8.07, 9.10, 8.79, 8.89, 8.36, 8.75, 7.8, 6.83, 9.31, 
        8.61, 7.5, 6.9, 9.82, 6.72, 8.73, 9.02, 8.41, 9.0, 8.41, 9.28, 
        9.04, 8.42, 8.43, 8.83, 8.3, 8.15, 9.44, 9.73, 8.07, 6.46, 
        8.82, 8.98, 7.99, 8.65, 8.93, 8.72, 7.25, 9.57, 7.28, 9.18, 
        9.42, 7.76, 9.62, 6.63, 8.04, 8.07, 8.22, 9.85, 9.14, 8.89, 
        8.91, 6.48, 8.63, 7.48, 9.62, 9.13, 8.73, 8.19, 9.61, 8.94, 
        8.74, 7.68, 7.59, 8.87, 9.24, 7.14, 7.59, 9.17, 9.16, 9.3, 
        8.93, 7.31, 9.31, 8.72, 8.78, 8.39])

T3 = np.array(["SY1", "SY1", "SY1", "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", 
         "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", "RQQ", "RLQ", "RLQ", "RLQ", 
         "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", 
         "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ",
         "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", 
         "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", 
         "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", "RLQ", 'RLQ', 'RLQ', 'RLQ', 'RLQ', 
         'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 
         'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 
         'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 
         'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 
         'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 
         'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 
         'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ', 'RLQ'])

Z3 = np.array([0.036, 0.013, 0.005, 0.171, 0.164, 0.109, 0.155, 0.176, 0.190, 0.185,
              0.086, 0.177, 0.165, 0.184, 0.164, 0.267, 0.146, 0.406, 0.333, 0.717,
              0.395, 0.343, 0.637, 0.425, 0.859, 0.367, 0.831, 0.258, 0.226, 0.669,
              0.186, 0.510, 2.224, 0.888, 0.852, 0.571, 0.574, 0.915, 0.454, 0.781, 0.444, 
              0.405, 0.954, 0.194, 0.759, 0.545, 0.452, 0.324, 0.654, 0.455, 0.846, 0.191, 
              0.631, 0.871, 0.684, 0.668, 0.427, 0.052, 0.698, 0.348, 0.712, 0.901, 0.24, 
              0.612, 0.794, 0.197, 0.312, 0.525, 0.595, 0.311, 0.355, 0.157, 0.497, 0.734, 
             0.266, 0.554, 0.656, 0.334, 0.258, 0.381, 0.789, 0.24, 0.751, 0.633, 0.321, 
             0.536, 0.19, 0.286, 0.332, 0.72, 0.313, 0.803, 0.368, 0.314, 0.905, 0.219, 
             0.361, 0.266, 0.412, 0.097, 1.401, 0.988, 0.75, 0.594, 0.751, 0.879, 0.449, 
             0.206, 0.293, 0.714, 0.691, 0.657, 0.46, 0.302, 0.303, 0.626, 0.24, 0.104, 
             1.012, 0.524, 0.932, 0.501, 0.2, 0.213, 0.698, 0.672, 0.298, 0.901, 0.655, 
             0.237, 0.859, 0.926, 0.741, 0.671, 0.735, 0.673, 0.576, 0.173, 0.21])


#Stellar velocity dispersions
L4 = np.log10((10**np.array([44.45, 43.67, 43.54, 43.54, 44.61, 44.27, 42.52, 44.33, 
               43.84, 45.04, 44.02, 44.37, 43.38, 44.69, 44.1, 44.05, 43.92, 43.08, 
               44.27, 43.47, 43.64, 43.38, 43.79, 45.39, 43.03, 43.81, 44.12, 
               43.04, 44.05, 43.6, 44.3, 44.19, 44.66, 43.86, 43.93, 43.6, 44.2, 
               44.54, 44.59, 43.37, 44.27, 45.15, 44.44, 44.52, 44.11, 44.75, 44.39, 
               44.53, 44.55, 44.27, 45.24, 44.84, 44.53, 44.54, 44.13, 44.39, 44.48
               ]) / 3.826E+33))

M4 = np.array([6.92, 8.21, 6.09, 8.95, 7.47, 7.02, 7.65, 7.51, 8.19, 8.51, 7.18, 
               7.88, 7.24, 7.88, 8.3, 7.3, 7.72, 6.06, 6.77, 7.53, 6.83, 7.4, 
               6.95, 8.04, 6.51, 6.79, 6.39, 7.25, 6.94, 7.6, 7.99, 7.38, 8.08, 
               6.88, 7.28, 6.59, 7.16, 8.65, 7.87, 7.6, 7.21, 7.56, 7.28, 6.92, 
               7.56, 7.62, 8.09, 7.64, 7.01, 6.83, 7.54, 8.0, 7.74, 8.23, 7.15, 
               7.69, 7.7])

T4 = np.array(['SY1', 'SY1', 'SY1', 'SY1', 'SY1', 'SY1', 'SY2', 'SY2', 'SY2', 
               'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 
               'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 
               'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 
               'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 
               'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 'SY2', 
               'SY2', 'SY2', 'SY2']) 

Z4 = np.array([0.005, 0.002, 0.004, 0.003, 0.029, 0.016, 0.002, 0.014, 0.005, 
               0.018, 0.009, 0.013, 0.003, 0.015, 0.008, 0.006, 0.008, 0.004, 0.028, 0.009, 
               0.003, 0.004, 0.002, 0.023, 0.004, 0.008, 0.009, 0.008, 0.007, 0.028, 0.006,
               0.023, 0.03, 0.013, 0.017, 0.006, 0.016, 0.014, 0.037, 0.01, 0.015, 0.029, 0.017, 
               0.023, 0.014, 0.024, 0.042, 0.017, 0.012, 0.015, 0.018, 0.023, 0.011, 0.025, 0.016, 
               0.016, 0.029])
'''

'''
#Luminosity vs mass of black hole
Masses = np.linspace(10**6, 10**10, 100) * const.M_sun.to_value()
luminosities = []
LEdd = []
for i in Masses:
    LumAGN = (Flux(Temp2(i, AccR, R_midpoints(i), Rin(i)))) * Area(R_midpoints(i))
    LumAGNsum = np.sum(LumAGN, axis=0)
    TotLumAGN = scipy.integrate.trapezoid(LumAGNsum, freq)
    luminosities.append(TotLumAGN)
    L = L_edd(i)
    LEdd.append(L)

fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)
#ax1.set_ylim(28, 45)
#ax1.plot(np.log10(luminosities), np.log10((Masses)/const.M_sun.to_value()), linestyle='-', color='black', label='')
ax1.plot(np.log10(LEdd/const.L_sun.to_value()), np.log10((Masses)/const.M_sun.to_value()), linestyle='-', color='black', label='Eddington limit')
#ax1.errorbar(L_bol1, M_BH1, fmt='o', color='red', label='Observed Data')
#ax1.errorbar(log_L2, M_BH2, fmt='o', color='red', label='Observed Data')

colours = {
    "SY1": "blue",
    "SY2": "purple",
    "RQQ": "green",
    "RLQ": "red"
}

# Plot Spatially Resolved Kinematics
for t, l, m, z in zip(T1, L1, M1, Z1):
    ax1.scatter(l, m, color=colours[t], marker="*")

# Plot Reverberation Mapping
for t, l, m, z in zip(T2, L2, M2, Z2):
    ax1.scatter(l, m, color=colours[t], marker="s")

# Plot Optical Luminosity
for t, l, m, z in zip(T3, L3, M3, Z3):
    ax1.scatter(l, m, color=colours[t], marker="^")

# Plot Stellar Velocity Dispersion
for t, l, m, Z in zip(T4, L4, M4, Z4):
    ax1.scatter(l, m, color=colours[t], marker="D")
'''

'''
# Line of best fit
coefficients = np.polyfit(L_bol, M_BH, 1)
polynomial = np.poly1d(coefficients)
x_fit = np.linspace(min(M_BH), max(M_BH), 100)
y_fit = polynomial(x_fit)
ax1.plot(x_fit, y_fit, linestyle='--', color='green', label='Best Fit Line')
'''
#ax1.set_xscale('log')
#ax1.set_yscale('log')
'''
legend_handles = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='Seyfert 1'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=10, label='Seyfert 2'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label='RQQ'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='RLQ')
]

# Add custom legend handles to the plot
ax1.legend(handles=legend_handles, fontsize=12)

ax1.set_ylabel(r'$\log\left(\frac{Mass}{M_\odot}\right)$', fontsize=16)
ax1.set_xlabel(r'$\log\left(\frac{Luminosity}{L_\odot}\right)$', fontsize=16)
plt.show()



#Redshift vs mass of black hole
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)
#ax1.plot(Z1, M1, linestyle='-', color='black', label='')
#ax1.plot(Z2, M2, linestyle='-', color='black', label='')
#ax1.plot(Z3, M3, linestyle='-', color='black', label='')
#ax1.plot(Z4, M4, linestyle='-', color='black', label='')

colours = {
    "SY1": "blue",
    "SY2": "orange",
    "RQQ": "purple", #radio quiet quasar
    "RLQ": "green" #radio loud quasar
}

# Plot Spatially Resolved Kinematics
for t, l, m, z in zip(T1, L1, M1, Z1):
    ax1.scatter(z, l, color=colours[t], linestyle='-')

# Plot Reverberation Mapping
for t, l, m, z in zip(T2, L2, M2, Z2):
    ax1.scatter(z, l, color=colours[t], linestyle='--')

# Plot Optical Luminosity
for t, l, m, z in zip(T3, L3, M3, Z3):
    ax1.scatter(z, l, color=colours[t], linestyle=':')

# Plot Stellar Velocity Dispersion
for t, l, m, Z in zip(T4, L4, M4, Z4):
    ax1.scatter(z, l, color=colours[t], linestyle='-.')

ax1.legend(fontsize=12)
ax1.set_ylabel('Mass', fontsize=16)
ax1.set_xlabel('Redshift', fontsize=16)
plt.show()





#mass vs accretion rate

Ac1 = ((10**L_bol1) * const.L_sun.to_value())/(0.1 * const.c.value**2)
Ac2 = ((10**log_L2) * const.L_sun.to_value())/(0.1 * const.c.value**2)


Accretion_rate = luminosities/(eta * const.c.to_value()**2)


fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)
ax1.plot(np.log10(AccR_Edd(Masses)), np.log10((Masses)/const.M_sun.to_value()), linestyle='-', color='black', label='Eddington limit')
#ax1.plot(np.log10((Masses)/const.M_sun.to_value()), np.log10(Accretion_rate), linestyle='-', color='black', label='Eddington limit')
ax1.errorbar(np.log10(Ac1), M_BH1, fmt='o', color='red', label='Observed Data')
ax1.errorbar(np.log10(Ac2), M_BH2, fmt='o', color='green', label='Observed Data')
#Line of best fit
#coefficients = np.polyfit(M_BH1, Ac, 1)
#polynomial = np.poly1d(coefficients)
#x_fit = np.linspace(min(M_BH1), max(M_BH1), 100)
#y_fit = polynomial(x_fit)
#ax1.plot(x_fit, y_fit, linestyle='--', color='green', label='Best Fit Line')
ax1.legend(fontsize=12)
ax1.set_ylabel('Log(Mass)', fontsize=16)
ax1.set_xlabel('Log(Accretion Rate)', fontsize=16)
'''


'''
#Luminosity vs accretion rate of black hole
luminosities2 = []
Accretion = np.linspace(0.00001, 2, 100) * AccR_Edd(MassSMBH2)
for i in Accretion:
    LumAGN = (Flux(Temp2(MassSMBH2, i, R_midpoints(MassSMBH2), Rin(MassSMBH2)))) * Area(R_midpoints(MassSMBH2))
    LumAGNsum = np.sum(LumAGN, axis=0)
    TotLumAGN = scipy.integrate.trapezoid(LumAGNsum, freq)
    luminosities2.append(TotLumAGN)
    L = L_edd(MassSMBH2)

print(luminosities2)

fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

ax1.plot(np.log10((luminosities2)/const.L_sun.to_value()), np.log10(LEdd), linestyle='-', color='black', label='')
ax1.errorbar(L_bol1, np.log10(Ac1), fmt='o', color='red', label='Observed Data')
ax1.errorbar(log_L2, np.log10(Ac2), fmt='o', color='green', label='Observed Data')
ax1.set_ylabel('Log(Accretion Rate)', fontsize=16)
ax1.set_xlabel('Log(Luminosity)', fontsize=16)
'''

#plt.show()







'''
#mag_absolute = mag- 5 * np.log10(dist - 1)

#Distance Array
dist1 = np.linspace(1, 1000, 1000) * 10e6

#flux_ref = #changes depending on filter used
#values given in erg/(s*cm^2*Hz)
flux_ref_U = 1.81e-20
flux_ref_B = 4.26e-20 
#* u.erg / (u.s * u.cm**2 * u.Hz)
flux_ref_V = 3.64e-20
flux_ref_R = 3.08e-20
flux_ref_I = 1.554e-20

#Max magnitude of possible observation
mag_faintest_Hubble = 31.5
mag_faintest_8to10m = 27
#faintest luminosity possible
def Lum_f_Hubble(dist):
    return Lum(mag_faintest_Hubble, dist)
def Lum_f_8to10m(dist):
    return Lum(mag_faintest_8to10m, dist)

#Hubble law
def Velocity(dist):
    return cosmo.H(0) * dist

#redshift equation (v << c)
def redshift(dist):
  return Velocity(dist)/const.c.to('km/s')

#Smallest black hole that can be observed
#Masses vs Magnitude
#would mean set distance

MassAGN = np.linspace(10**6, 10**10, 100) * const.M_sun.to_value()

luminosities = []

for i in MassAGN:
    LumAGN = (Flux(Temp2(i, AccR_Edd, R_midpoints, Rin))) * Area(R_midpoints)
    LumAGNsum = np.sum(LumAGN, axis=0)
    TotLumAGN = scipy.integrate.trapezoid(LumAGNsum, freq)
    luminosities.append(TotLumAGN)

#flux received from AGN
def f(Lum, dist):
  return Lum/(4*np.pi * dist**2)

#luminosity of zero magnitude star
L0 = 3*10e28

#flux to magnitude with no filter with respect to a zero point source
def mag(Lum, dist):
  return -2.5 * np.log10(f(Lum, dist)/(L0))


def efficiency():
  return (1/4)*(Rg/Rin)

def ACCR(dist):
  return Lum_f_Hubble(dist) / (efficiency() * const.c.to_value()**2)

'''






#start with black hole in a close galaxy or maybes even in milky way -> move further awawy until no longer visible -> increase mass 
#repeat until mass becomes too large for black hole -> change to AGN??????
#lighteset possible AGN observable if it were to be accreting at eddington accretion rate at the edge of the observalbe universe
#Can we find anything about the early universe and the first formed galaxies????
#AGN accretion rate??
#
