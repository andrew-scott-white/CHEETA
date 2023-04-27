from gpkit import Model, Variable, units, Vectorize, SignomialsEnabled

from path import h_alt, hdot_in, M, P, T, rho, mu, Cp, LDR, C_D_fuse, C_D_wing
from flight_state import FlightState

from mission import Mission
from aircraft import Aircraft
from turbofan import Turbofan

import pandas as pd  

import matplotlib.pyplot as plt
import numpy as np

gamma = 1.4 # specific heat ratio for air
g = 9.81*units('m/s/s') # gravitational constant
R1 = 287*units('J/kg/K') # gas constant for air
R2 = 8.314*units('J/K') # universal gas constant
F = 96485.33*units('C') # Faraday's constant

BLI = 'on'
               
# set up propulsion system
tf = Turbofan(2)
tf.substitutions.update({
    # gas generator constants
    #'(P/m)_{gg}': 7687 * units('hp') / 3065 / units('lb'), # half of LEAP-1B
    '(P/m)_{gg}': 6380 * units('hp') / (5216/2) / units('lb'), # assuming half of CFM-56 weight is gas turbine and half is fan
    # ducted fan constants
    #'K_{fan}': 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
    'K_{fan}': (5216/2) * units('lb') / (61 * units('in'))**3, # half of CFM-56 weight is turbine and half is fan
    'HTR': 0.3,
})

# set up aircraft system (to which propulsion system is passed as input)
ac = Aircraft(tf)
ac.substitutions.update({
    'f_{empty}': 0.47  , # 737-800 at 3000 mi mission for 35000 lb payload
    'b_{fuse}': 14*units('ft'),
    'b_{wing}': (136-14)*units('ft'),
    'S_{ref}': 2000*units('ft^2'),
    'd_{fan}': 61*units('in'), # CFM-56 engine fan diameter
    })

mission = Mission(ac, h_alt, M, P, T, rho, mu, Cp, hdot_in, LDR, C_D_fuse, C_D_wing, BLI)
# all aircraft and propulsion system performance models (and thus their
# variables) get set up inside Mission(), so we update their constants here

mission.substitutions.update({
    '\rho_{fuel}':              840*units('kg/m^3'), # jet fuel density
    'h_{fuel}':                43 * units('MJ/kg'), # jet fuel lower heating value
    
# mission parameters
    'f_{res}': 1.2, # assume 20% additional fuel for reserves
    'R_{des}': 5500 * units('km'), # design range
    'm_{pay}':  15875 * units('kg'), # design payload
})

mission.substitutions.update({
   # 'D_corr':           0.5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
    'L/D':              LDR, # take input lift-to-drag ratios from Elias
    #'f_{BLI}':                   0.0001, # assuming no BLI
        '\phi_{fuse}':               .5, # assume half of fuselage dissipation occurs on top of fuselage 
        '\phi_{wing}':               .64, # assume most dissipation occurs on top of wing
    '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]; according to D.K. Hall, this number is generally representative         
    #'D_p':                       4830*units('lbf'), # profile drag ingested by propulsors
    '\\eta_{th}':                0.55,
    #'\\eta_{fan}':               0.95, 
    #'Cd':                        0.04, # [Robinson, 2017]
    
    })

model = Model(mission['m_{fuel,tot}'], [ac, mission])

init = {
    mission['F_{net}']:                                         44*units('kN'),
    mission.fs.performance.prop.fan_perf['u_{jet}']:       100*units('m/s'),
    mission.fs.performance.prop.fan_perf['\\dot{m}']:       1000000*units('kg/s'),
    mission['dt']:                                           1000*units('s'),
    mission.fs.performance.prop.fan_perf['F_{fan}']:       20*units('kN'),
    mission['F_{net}']:                                         44*units('kN'),
    mission['K_{inl}']:                                         300*units('kW'),
    
    }

if BLI == 'on':
    sol = model.localsolve(solver='mosek_conif', x0=init)
else:
    sol = model.solve(solver='mosek_conif', x0=init)
print(sol.summary())
#print(sol.table())
print("PFEI (kJ/kg/km): " + str(sol['variables']['PFEI']))
print("MTOW (kg): " + str(sol['variables']['MTOW']))
print("OEW (kg): " + str(sol['variables']['m_{OE}']))
print("Fuel Weight (kg): " + str(sol['variables']['m_{fuel,tot}']))
print("d_{fan} (in): " + str(sol['variables']['d_{fan}']*39.37))
print("f_{wing,BLI}: " + str(sol['variables']['f_{wing,BLI}']))
print("f_{fuse,BLI}: " + str(sol['variables']['f_{fuse,BLI}']))