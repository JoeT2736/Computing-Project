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





#HOW DOES ACCRETION RATE OF AGN CHANGE WITH REDSHIFT??????????????????????????????????????????????????????????????????????




MassSMBH1 =  10**7 * const.M_sun.to_value()
MassSMBH2 =  10**9 * const.M_sun.to_value()
MassSMBH3 =  10**11 * const.M_sun.to_value()
MassBH =  10 * const.M_sun.to_value()  #Mass of Black hole
Mass2 = 10**4 * const.M_sun.to_value()
AccR = 10**15 #*u.kg/u.s  #Accretion rate
AccR2 = 10**10 #*u.kg/u.s
Mr2 = 10**14 #*u.kg/(u.s*u.m**3)

Fstart = 10
Fstop = 25
Fsteps = 1000
freq = np.logspace(Fstart, Fstop, Fsteps) #*u.Hz #Range of frequencies


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
  a = (3 * const.G.to_value()* M * Ar)
  b = (8 * np.pi * const.sigma_sb.to_value() * Radius**3)
  c = ( Rin / Radius )**(1/2)
  T = ((a / b) * (1 - c))**(1/4)
  return T
#print(Temp2(MassBH, AccR, R_midpoints, Rin))



#Temperatue equation for scaled units
#Temperatue equation for scaled units
def Temp(M, Ar, Radius, RIN):
  a = (3 * const.G.to_value()* M * Ar)
  b = (8 * np.pi * const.sigma_sb.to_value() * Radius**3 * Rg**3)
  c = ( RIN / Radius )**(1/2)
  T = ((a / b) * (1 - c[:, np.newaxis]))**(1/4)
  return T
#print(Temp(MassBH, AccR, r_midpoints, rin))

#Temp without viscous forces
def Temp3(M, Ar, Radius):
  a = (const.G.to_value()* M * Ar)
  b = (8 * np.pi * const.sigma_sb.to_value() * Radius**3)
  T = (a / b)**(1/4)
  return T

#T2 = Temp2(MassBH, AccR, R_midpoints, Rin)
#T3 = Temp3(MassBH, AccR, R_midpoints)

#print(T2[5000], T3[5000])



#GOOD TEMPERATURE PLOT MIGHT NEED THIS
'''
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 8)

fig.patch.set_facecolor('#D5D5D5') 
ax1.set_facecolor('#D5D5D5') 

ax1.plot(Rin(MassBH) + R_midpoints(MassBH), Temp2(MassBH, AccR, R_midpoints(MassBH), Rin(MassBH)), linestyle = '--', color = 'blue', label = 'Viscous')
ax1.plot(Rin(Mass2) + R_midpoints(Mass2), Temp2(Mass2, AccR, R_midpoints(Mass2), Rin(Mass2)), color = 'red', label = '')
ax1.set_xlim(0, 10e9)
ax1.set_ylabel('Temperature (K)', fontsize=20)
ax1.set_xlabel('Distance from centre of Black hole (m)', fontsize=20)
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.legend(fontsize=16)
plt.show()
'''






#Blackbody flux not using built in function
def Flux(T):
   temps = T.reshape(-1, 1)
   freqs = freq.reshape(1, -1)
   a = (2 * np.pi *const.h.to_value() * freqs**3)/(const.c.to_value()**2)
   b = (const.h.to_value() * freqs)/(const.k_B.to_value() * temps)
   c = np.exp(b)
   Fv = a/(c-1)
   return Fv





#print(Flux(Temp(MassBH, AccR, r_midpoints, rin)))

#Using function
#LL = F * Area(Log_rMid)
#Lsum = np.sum(LL, axis = 0)
#LogL = F * Area(Log_rMid)
#LogLsum = np.sum(LL, axis = 0)




#Milestone values
LumMilestone = (Flux(Temp2(MassBH, AccR, R_midpoints(MassBH), Rin(MassBH)))) * Area(R_midpoints(MassBH))  
LumMilestonesum = np.sum(LumMilestone, axis = 0)

'''
#milestone values, max spin black hole
L3 = (Flux(Temp2(MassBH, AccR, R_midpoints2, Rin2))) * Area(R_midpoints2) 
L3sum = np.sum(L3, axis = 0)
'''


#Object mass = mass of sun
LumSunMass = (Flux(Temp2(Mass2, AccR, R_midpoints(Mass2), Rin(Mass2)))) * Area(R_midpoints(Mass2))
LumSunMass2um = np.sum(LumSunMass, axis = 0)


'''
#mass of sun with spin
LumSunMass2pin = (Flux(Temp2(Mass2, AccR, R_midpoints2, Rin2))) * Area(R_midpoints2)
LumSunMass2umSpin = np.sum(LumSunMass2pin, axis = 0)
'''

'''
#Different accretion rate
LumdAccr = (Flux(Temp2(MassBH, AccR2, R_midpoints, Rin))) * Area(R_midpoints)
LumdAccrSum = np.sum(LumdAccr, axis = 0)
'''

'''
#D Accr with spin
LumdAccrSpin = (Flux(Temp2(MassBH, AccR2, R_midpoints2, Rin2))) * Area(R_midpoints2)
LumdAccrSumSpin = np.sum(LumdAccrSpin, axis = 0)
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
TotLumSunMass = scipy.integrate.trapezoid(LumSunMass2um, freq)*u.W  #.to(u.W)
#TotLumSunMass2pin = scipy.integrate.trapezoid(LumSunMass2umSpin, freq)*u.W  #.to(u.W)
#TotLumdAccr = scipy.integrate.trapezoid(LumdAccrSum, freq)*u.W  #.to(u.W)
#TotLumdAccrSpin = scipy.integrate.trapezoid(LumdAccrSumSpin, freq)*u.W  #.to(u.W)
#TotLumMilestone2 = scipy.integrate.trapezoid(Lsum_non_viscous, freq)*u.W  #.to(u.W)
TotLumSMBH1 = scipy.integrate.trapezoid(LumSMBH1sum, freq)*u.W
TotLumSMBH2 = scipy.integrate.trapezoid(LumSMBH2sum, freq)*u.W
TotLumSMBH3 = scipy.integrate.trapezoid(LumSMBH3sum, freq)*u.W
TotLumEdd = scipy.integrate.trapezoid(LumEddsum, freq)*u.W


'''
print(f'(2) L sum = {TotLumSunMass}')
print(f'(Milestone) L sum = {TotLumMilestone}')
print(f'(1) L sum = {TotLumSMBH1}')
print(f'(2) L sum = {TotLumSMBH2}')
print(f'(3) L sum = {TotLumSMBH3}')
print(f'(Eddington) L sum = {TotLumEdd}')
'''


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
e, = ax1.plot(np.log10(freq), np.log10(freq * LumEddsum), linestyle='-', color='purple', label='Eddington Limit (10)')
a, = ax1.plot(np.log10(freq), np.log10(freq * LumSMBH1sum), linestyle='-', color='green', label='10^7')
b, = ax1.plot(np.log10(freq), np.log10(freq * LumSMBH2sum), linestyle='-', color='blue', label='10^9')
c, = ax1.plot(np.log10(freq), np.log10(freq * LumSMBH3sum), linestyle='-', color='red', label='10^11')
f, = ax1.plot(np.log10(freq), np.log10(freq * LumSunMass2um), linestyle='-', color='brown', label='10^4')


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
''''
accretion_rates = [0.1, 1, 10]  # Different accretion rates (in units of Eddington accretion rate)

LLowEdd = (Flux(Temp2(MassBH, 0.1*AccR_Edd, R_midpoints, Rin))) * Area(R_midpoints)  
LLowEddSum = np.sum(LLowEdd, axis = 0)
LAboveEdd = (Flux(Temp2(MassBH, 3*AccR_Edd, R_midpoints, Rin))) * Area(R_midpoints)  
LAEddSum = np.sum(LAboveEdd, axis = 0)
LFarLessEdd = (Flux(Temp2(MassBH, 0.0001*AccR_Edd, R_midpoints, Rin))) * Area(R_midpoints)  
LFLEddSum = np.sum(LFarLessEdd, axis = 0)
LumEdd2 = (Flux(Temp2(30*MassBH, AccR_Edd, R_midpoints, Rin))) * Area(R_midpoints)  
LumEddsum2 = np.sum(LumEdd2, axis = 0)
LumEdd3 = (Flux(Temp2(60*MassBH, AccR_Edd, R_midpoints, Rin))) * Area(R_midpoints)  
LumEddsum3 = np.sum(LumEdd3, axis = 0)

fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

ax1.set_xlim(11, 20)
ax1.set_ylim(15, 38)

lowEdd, = ax1.plot(np.log10(freq), np.log10(freq * LLowEddSum), color='blue', linestyle='-', label='Accretion rate = 0.1 Edd')
Edd1, = ax1.plot(np.log10(freq), np.log10(freq * LumEddsum), color='red', linestyle='-', label='Accretion rate = 1 Edd')
Edd2, = ax1.plot(np.log10(freq), np.log10(freq * LumEddsum2), color='red', linestyle=':')
Edd3, = ax1.plot(np.log10(freq), np.log10(freq * LumEddsum3), color='red', linestyle='--')
AboveEdd, = ax1.plot(np.log10(freq), np.log10(freq * LAEddSum), color='purple', linestyle='-', label='Accretion rate = 3 Edd')
FarLessThanEdd, = ax1.plot(np.log10(freq), np.log10(freq * LFLEddSum), color='green', linestyle='-', label='Accretion rate = 0.001 Edd')

mass1_line = plt.Line2D([0], [0], color='black', linestyle='-', label='Mass = MassBH')
mass30_line = plt.Line2D([0], [0], color='black', linestyle=':', label='Mass = 30*MassBH')
mass60_line = plt.Line2D([0], [0], color='black', linestyle='--', label='Mass = 60*MassBH')

ax1.legend(handles=[FarLessThanEdd, lowEdd, Edd1, Edd2, Edd3, AboveEdd, 
mass1_line, mass30_line, mass60_line], fontsize=12)


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

print(luminosities)


#LogMass vs LogLuminosity

fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

ax1.plot((Masses)/const.M_sun.to_value(), luminosities, linestyle='-', color='black', label='')
ax1.plot((Masses)/const.M_sun.to_value(), LEdd, linestyle='-', color='black', label='Eddington limit')

ax1.set_xscale('log')
ax1.set_yscale('log')

#ax1.set_ylabel(r'$\log_{10}(Luminosity)$ W', fontsize=16)




luminosities2 = []
Accretion = np.linspace(0.001, 10, 100) * AccR_Edd(MassSMBH2)


for i in Accretion:
    LumAGN = (Flux(Temp2(MassSMBH2, i, R_midpoints(MassSMBH2), Rin(MassSMBH2)))) * Area(R_midpoints(MassSMBH2))
    LumAGNsum = np.sum(LumAGN, axis=0)
    TotLumAGN = scipy.integrate.trapezoid(LumAGNsum, freq)
    luminosities2.append(TotLumAGN)
    L = L_edd(MassSMBH2)

print(luminosities2)



fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

ax1.plot((Accretion)/AccR_Edd(MassSMBH2), (luminosities2)/L, linestyle='-', color='black', label='')


ax1.set_xscale('log')
ax1.set_yscale('log')



plt.show()








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




#Magnitude to flux equation (flux in watts per square kilometer)
def flux(mag):
    return 10**((-mag)/2.5) * flux_ref_B
    #.to(u.J / (u.km**2 * u.s * u.Hz)) #change this depending on filter used



#Luminosity of AGN from the apparent magnitude
def Lum(mag, dist):
    return flux(mag) * (4*np.pi * dist**2)


#Max magnitude of possible observation
mag_faintest_Hubble = 31.5
mag_faintest_8to10m = 27


#faintest luminosity possible
def Lum_f_Hubble(dist):
    return Lum(mag_faintest_Hubble, dist)

def Lum_f_8to10m(dist):
    return Lum(mag_faintest_8to10m, dist)



#Using luminosity to find Accretion rate 



#Hubble law
def Velocity(dist):
    return cosmo.H(0) * dist


#redshift equation (v << c)
def redshift(dist):
  return Velocity(dist)/const.c.to('km/s')






def Lum_f_Hubble(dist):
    return Lum(mag_faintest_Hubble, dist)

def Lum_f_8to10m(dist):
    return Lum(mag_faintest_8to10m, dist)







'''
#Accretion rate from temperature 
def Accretion_Rate(Mass, Radius, Rin, Temp):#i dont think this is going to work (would need to know size of black hole/AGN)
  a = 8*np.pi * Temp * const.sigma_sb.to_value() * Radius**3
  b = 1 - np.sqrt(Rin/Radius)
  c = 3 * const.G.to_value() * Mass
  return a / (b * c)


# accretion rate proportional to (np.ln(luminosity))**4



#plot of flux vs redshift
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

#ax1.set_xlim(11, 20)
#ax1.set_ylim(15, 38)


for i in dist:
  R = redshift(i)
  L1, = ax1.plot(R, flux(mag_faintest_Hubble), linestyle='-', color='blue', label='Hubble')
  L2, = ax1.plot(R, flux(mag_faintest_8to10m), linestyle='-', color='red', label='8-10m')
  
ax1.legend(fontsize=12)

plt.show()
'''








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

#in parsecs
distance = 1000*10e7

magnitudes = [mag(lum, distance) for lum in luminosities]



fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

#ax1.set_xlim(11, 20)
#ax1.set_ylim(15, 38)

ax1.plot(magnitudes, MassAGN/const.M_sun.to_value(), linestyle='-', color='black', label='')
ax1.axvline(x=mag_faintest_Hubble, color='red', linestyle='--', label=f'Faintest mag = {mag_faintest_Hubble}'
, ymin = 0, ymax = 2.2e40/const.M_sun.to_value())
ax1.axvline(x=mag_faintest_8to10m, color='red', linestyle='--', label=f'Faintest mag = {mag_faintest_8to10m}'
, ymin = 0, ymax = 2.2e40/const.M_sun.to_value())


closest_index = np.argmin(np.abs(np.array(magnitudes) - mag_faintest_Hubble))
intersection_mass = MassAGN[closest_index]
#print(intersection_mass/const.M_sun.to_value())

ax1.legend()
#plt.show()


def efficiency():
  return (1/4)*(Rg/Rin)

def ACCR(dist):
  return Lum_f_Hubble(dist) / (efficiency() * const.c.to_value()**2)

Lum(mag, dist1)
ACCR(dist1)



fig, ax1 = plt.subplots()
fig.set_size_inches(10, 6)

#ax1.set_xlim(11, 20)
#ax1.set_ylim(15, 38)

ax1.plot(redshift(dist1), ACCR(dist1), linestyle='-', color='black', label='')



ax1.legend()
plt.show()









#start with black hole in a close galaxy or maybes even in milky way -> move further awawy until no longer visible -> increase mass 
#repeat until mass becomes too large for black hole -> change to AGN??????
#lighteset possible AGN observable if it were to be accreting at eddington accretion rate at the edge of the observalbe universe
#Can we find anything about the early universe and the first formed galaxies????
#AGN accretion rate??
#


