from gpkit import Model, Variable, units, Vectorize, SignomialsEnabled

from path import h_alt, hdot_in, M, P, T, rho, mu, Cp, LDR, C_D_fuse, C_D_wing
from flight_state import FlightState

from mission import Mission
from aircraft import Aircraft
from fully_electric_bat import FullyElectricBatHybrid

import pandas as pd  

import matplotlib.pyplot as plt
import numpy as np

gamma = 1.4 # specific heat ratio for air
g = 9.81*units('m/s/s') # gravitational constant
R1 = 287*units('J/kg/K') # gas constant for air
R2 = 8.314*units('J/K') # universal gas constant
F = 96485.33*units('C') # Faraday's constant

BLI = 'on'
N_b = 1
m_bat = np.asarray([1, 10, 100, 1000]) # kg
m_bat_tot = N_b*m_bat

for ii in range(len(m_bat_tot)):
    m_bat_tot[ii] = int(m_bat_tot[ii])
    
title = ['No Battery', str(m_bat_tot[1]) + ' kg Battery', str(m_bat_tot[2]) + ' kg Battery', str(m_bat_tot[3]) + ' kg Battery']
F_HEX = np.zeros([len(m_bat),len(h_alt)])
D_core = np.zeros([len(m_bat),len(h_alt)])
D_nac = np.zeros([len(m_bat),len(h_alt)])
D_net = np.zeros([len(m_bat),len(h_alt)])
Q = np.zeros([len(m_bat),len(h_alt)])
V = np.zeros([len(m_bat),len(h_alt)])
eta_fc = np.zeros([len(m_bat),len(h_alt)])
eta_net = np.zeros([len(m_bat),len(h_alt)])
F_net = np.zeros([len(m_bat),len(h_alt)])
P_gross = np.zeros([len(m_bat),len(h_alt)])
PFEI = np.zeros([len(m_bat)])

for ii in range(len(m_bat)):
               
    # set up propulsion system
    feb = FullyElectricBatHybrid(9)
    feb.substitutions.update({
        # motor constants
        '(T/m)_{mot}': 70 * units('N*m/kg'),  # 70 for superconducting; 20 for traditional
        # ducted fan constants
        #'K_{fan}': 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
        'K_{fan}': (5216/2) * units('lb') / (61 * units('in'))**3, # half of CFM-56 weight is turbine and half is fan
        'HTR': 0.38, 
        'd_{fan}': 1.2*units('m'), # CHEETA fan diameter
        
        # fuel cell constants
        'P_{max,fc}':               850*units('kW'), # max power the model is accurate for. This could be replaced with a max power specified by Boeing
        'm_{comp}':                 36*units('kg'),
        #'N_s':                      6,
        # HEX constants
        '\psi':                    5.02*units('kg/ft^2'), # from Boeing (Chellappa)
        # battery constants
        '(P/m)_{bat}': 800*units('W/kg'), # conservative 2035 future estimate from LEARN report pg. 14
        '(E/m)_{\\rm bat}': 8*10**2*units('W*hr/kg'), # 2050 estimate of Li ion batteries based on 4% annual growth rate [Viswanathan, Nature, 2022]
        'N_b':  N_b,
    })

    if m_bat[ii] > 0:
        feb.substitutions.update({
        'm_{bat}': m_bat[ii]*units('kg'),
        })
    
    # set up aircraft system (to which propulsion system is passed as input)
    ac = Aircraft(feb)
    ac.substitutions.update({
        'f_{empty}': 0.554  , # assumed
        'b_{fuse}': 14*units('ft'),
        'b_{wing}': (136-14)*units('ft'),
        'S_{ref}': 2000*units('ft^2'),
        #'MTOW': 190000*units('lb'),
        })
    
    mission = Mission(ac, h_alt, M, P, T, rho, mu, Cp, hdot_in, LDR, C_D_fuse, C_D_wing, BLI)
    # all aircraft and propulsion system performance models (and thus their
    # variables) get set up inside Mission(), so we update their constants here
    
    mission.substitutions.update({
        '\rho_{fuel}':              71*units('kg/m^3'), # jet fuel density
        'h_{fuel}':                120 * units('MJ/kg'), # jet fuel lower heating value
        
    # mission parameters
        'f_{res}': 1.2, # assume 20% additional fuel for reserves
        'R_{des}': 5500 * units('km'), # design range
        'm_{pay}':  15875 * units('kg'), # design payload
    })
    
    mission.substitutions.update({
        'L/D':              LDR, # take input lift-to-drag ratios from Elias
        '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]; according to D.K. Hall, this number is generally representative         
        #'D_p':                       4830*units('lbf'), # profile drag ingested by propulsors
        '\phi_{fuse}':               .5, # assume half of fuselage dissipation occurs on top of fuselage 
        '\phi_{wing}':               .64, # assume most dissipation occurs on top of wing
        
        # fan constants
        #'D_corr':           0.5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'\\eta_{fan}':               0.95, 
        #'Cd':                        0.04, # [Robinson, 2017]
        
        # fuel cell constants
        '(P/m)_{fc}':               2425*units('W/kg'), # CHEETA
        'P_{t\\fuel}':               5*units('bar'), # fuel storage pressure
        'A':                         446*units('cm**2'), # membrane area
        'N_{cells}':                 1080, # number of cells in a stack
        '\mu_f':                     .96, # fuel utilization
        'alpha':                     .5, # charge transfer coefficient (empirical)
        'c':                         .21, # fraction of O2 in air
        'S':                         1.5, # stoichiometry (surplus air)
        'eta_{comp}':                .9, # fuel cell compressor efficiency
        'T_{t\\2}':                  (180+273.15)*units('K'), # stack operating temperature
        '\\DeltaG':                  224.24*units('kJ'),
        'k':                         5e-5*units('1/s'),   # electrochemical rate constant
        '\\ohm_n':                   .062*units('ohm*cm**2'), # fuel cell internal resistance
        'x_{parasitic}':             .05,        
        'm':                         4e-5*units('V'), # Elias
        'n':                         2.3e-3*units('cm**2/mA'), # Elias
        
        # HEX constants
        'C_d':                        0.04, # [Robinson, 2017]
        'Pr':                         0.72,
        #'K_f':                        2,
        #'K_h':                        0.93,
        #'\sigma':                     .464,
        #'\epsilon':                   0.85, # Kays and London
        'C_{p,h}': 2200*units('J/kg/K'),
        '\dot{m}_h': 15*units('kg/s'),
        'h': 75*units('W/m^2/K'), # based on lowest value given by Chellappa (.5 MW heat load during cruise conditions)
        
        # motor constants
        'omega_{motor}': 4000 * units('rpm'),
        '1-eta_{motor}': 1 - .999, #.999 for superconducting; .95 for traditional
        })
    
    mission.append(mission['E_{bat,max}'] >= mission['E_{bat}'][0] + mission['E_{bat}'][1] + mission['E_{bat}'][2] + mission['E_{bat}'][3] + mission['E_{bat}'][4] + mission['E_{bat}'][5] + mission['E_{bat}'][6] + mission['E_{bat}'][7] + mission['E_{bat}'][8] + mission['E_{bat}'][9] + mission['E_{bat}'][len(h_alt)-1],)
    
    
    # these constraints turn off the battery for all operating conditions except top-of-climb
   # for yy in range(len(h_alt)):
      #  if yy == len(h_alt)-2:
      #      continue
      #  if yy == 0:
      #      continue
      #  else:
        #    mission.append(.01*units('kJ') >= mission['E_{bat}'][yy]) # battery expends negligible energy at each flight state yy
            
    
    init = {
        mission['(i0/A)']:                                          1e-3*units('A/cm**2'),
        mission['P_{comp}']:                                        1e12*units('kW'),
        mission['V_{act}']:                                         .3*units('V'),
        mission['T_{t\\2}']:                                        (180+273.15)*units('K'),
        mission['\\dot{m}_{a}']:                                    360*units('kg/s'),
        mission['eta_{comp}']:                                      .9,
        mission['F_{net}']:                                         44*units('kN'),
        mission.fs.performance.prop.fan_fuse_perf['u_{jet}']:       100000*units('m/s'),
        mission.fs.performance.prop.fan_wing_perf['u_{jet}']:       100000*units('m/s'),
        mission['pi_c']:                                            5,
        mission['F_{HEX}']:                                         1*units('kN'),
        mission['P_{gross}']:                                       1*units('kW'),
        mission.fs['dt']:                                           1000*units('s'),
        mission.fs.performance.prop.fan_fuse_perf['F_{fan}']:       100000*units('kN'),
        mission.fs.performance.prop.fan_wing_perf['F_{fan}']:       100000*units('kN'),
        mission.fs.performance.prop.bat_perf['P_{out}'][len(h_alt)-2]: 300*units('kW'),
        mission['m_{bat}']:                                         100*units('kg'),
        mission['N_s']:                                             20,
        mission['E_{bat}']:                                         1000*units('kJ'),
        }
    
    model = Model(mission['m_{fuel,tot}'], [ac, mission])
    
    sol = model.localsolve(x0=init,solver='mosek_conif',verbosity=2,iteration_limit=1000)
    #print(sol.table())
    
    P_out_keys = list(sol['variables']['P_{out}'].keys())
    P_out_values = list(sol['variables']['P_{out}'].values())
    
    F_HEX[ii] = sol['variables']['F_{HEX}']
    D_core[ii] = sol['variables']['D_{core}']
    D_nac[ii] = sol['variables']['D_{Nacelle}']
    D_net[ii] = D_core[ii] + D_nac[ii]- F_HEX[ii]
    Q[ii] = sol['variables']['Q_{HEX}']
    V[ii] = sol['variables']['V_\\infty']
    eta_fc[ii] = sol['variables']['eta_{fc}']
    eta_net[ii] = sol['variables']['eta_{net}']
    F_net[ii] = sol['variables']['F_{net}']
    N_s = sol['variables']['N_s']
    N_p = sol['variables']['N_p']
    P_gross[ii] = sol['variables']['P_{gross}']
    PFEI[ii] = sol['variables']['PFEI']
    S_ref = sol['variables']['S_{ref}']
    A_m = sol['variables']['A']
    dt = sol['variables']['dt']
    
    # plot below is misleading; thrust is higher because more is needed to carry the battery NOT because the battery helps
#plt.figure(1)
# h_alt*0.0003048 converts h_alt from ft to km
#plt.plot(F_HEX[ii]*V[ii]/Q[ii], h_alt*0.0003048, 'r--', label='HEX Core Thrust')
#plt.plot(-D_core[ii]*V[ii]/Q[ii], h_alt*0.0003048, 'b--', label='HEX Core Drag')
#plt.plot(-D_nac[ii]*V[ii]/Q[ii], h_alt*0.0003048, 'b:', label='HEX Nacelle Drag')
#plt.plot(-D_net[ii]*V[ii]/Q[ii], h_alt*0.0003048, 'k-', label='HEX Drag')
#plt.plot(F_net[0]*V[0]/Q[0]/N_s, h_alt*0.0003048, 'k-', label=title[0])
#plt.plot(F_net[1]*V[1]/Q[1]/N_s, h_alt*0.0003048, 'b-.', label=title[1])
#plt.plot(F_net[2]*V[2]/Q[2]/N_s, h_alt*0.0003048, 'b:', label=title[2])
#plt.plot(F_net[3]*V[3]/Q[3]/N_s, h_alt*0.0003048, 'b--', label=title[3])
#plt.legend()
#plt.xlabel('$\dfrac{F_x V_\infty}{Q}$')
#plt.ylabel('Altitude (km)')
#plt.xticks(np.linspace(1.25, 1.75, 6))
#plt.yticks(np.linspace(0,12,7))
#plt.title(title[ii])
 
plt.figure(2)
# h_alt*0.0003048 converts h_alt from ft to km
plt.plot(eta_fc[0], h_alt*0.0003048, 'k-', label=title[0])
plt.plot(eta_fc[len(m_bat)-1], h_alt*0.0003048, 'b-', label='Battery')
plt.legend()
plt.xlabel('Fuel Cell Efficiency')
plt.ylabel('Altitude (km)')
plt.xticks(np.linspace(0.60,.70,6))
plt.yticks(np.linspace(0,12,7))
 
plt.figure(3)
# h_alt*0.0003048 converts h_alt from ft to km  
plt.plot(P_gross[0]*1000*N_s/rho/V[0]**3/A_m,h_alt*0.0003048, 'k-', label=title[0])    
#plt.plot(P_gross[1]*1000*N_s/rho/V[1]**3/A_m,h_alt*0.0003048, 'b--', alpha=1, label=title[1])
plt.plot(P_gross[2]*1000*N_s/rho/V[2]**3/A_m,h_alt*0.0003048, 'r--', alpha=.5, label=title[2])
plt.plot(P_gross[3]*1000*N_s/rho/V[3]**3/A_m,h_alt*0.0003048, 'b--', alpha=.5, label=title[3])    
plt.legend()
plt.xlabel('$\dfrac{N_s \dot{W}_{fc}}{\dfrac{1}{2} \\rho V_\infty^3 A_m}$')
plt.ylabel('Altitude (km)')
plt.xticks(np.linspace(0, .03, 4))
plt.yticks(np.linspace(0,12,7))
#plt.title(title[ii])

    

    
    

print(sol.summary())

print("PFEI (kJ/kg/km): " + str(sol['variables']['PFEI']))
print("MTOW (kg): " + str(sol['variables']['MTOW']))
print("d_{fan} (in): " + str(sol['variables']['d_{fan}']*39.37))
print("f_{wing,BLI}: " + str(sol['variables']['f_{wing,BLI}']))
print("f_{fuse,BLI}: " + str(sol['variables']['f_{fuse,BLI}']))
print(" ")
print("E_tot (kJ): " + str(sol['variables']['E_{bat,max}']))
print("Battery Weight (kg): " + str(sol['variables']['m_{bat}']))
#print(sol['variables']['\\dot{m}'])
print("# of batteries: " + str(sol['variables']['N_b']))
print("# of stacks: " + str(sol['variables']['N_s']))

print(" ")
print("Output Power")
print(sol['variables']['P_{out}'])
print("Net Thrust (kN)")
print(sol['variables']['F_{net}'])
print("Battery Energy (kJ)")
print(sol['variables']['E_{bat}'])
print('HEX Effectiveness')
print(sol['variables']['\epsilon'])
print('Heat Load (kW)')
print(sol['variables']['N_s']*sol['variables']['Q_{HEX}'])
print('Fuel Cell Efficiency')
print(sol['variables']['eta_{fc}'])
print('HEX core drag (kN)')
print(sol['variables']['D_{core}'])
print('HEX heat thrust (kN)')
print(sol['variables']['F_{HEX}'])
print('HEX nacelle drag (kN)')
print(sol['variables']['D_{Nacelle}'])

