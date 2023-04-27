from gpkit import Model, Variable, units, Vectorize, SignomialsEnabled

from path import h_alt, hdot_in, M, P, T, rho, mu, Cp, LDR, C_D_fuse, C_D_wing
from flight_state import FlightState

from mission import Mission
from aircraft import Aircraft
from hybrid_electric import HybridElectric

import numpy as np

gamma = 1.4 # specific heat ratio for air
g = 9.81*units('m/s/s') # gravitational constant
R1 = 287*units('J/kg/K') # gas constant for air
R2 = 8.314*units('J/K') # universal gas constant
F = 96485.33*units('C') # Faraday's constant

BLI = 'on'

# set up propulsion system
he = HybridElectric(9)
he.substitutions.update({
    # gas generator constants
    #'(P/m)_{gg}': 7687 * units('hp') / 3065 / units('lb'), # half of LEAP-1B
    '(P/m)_{gg}': 6380 * units('hp') / (5216/2) / units('lb'), # assuming half of CFM-56 weight is gas turbine and half is fan
    # generator constants
    '(T/m)_{gen}': 20 * units('N*m/kg'),
    # motor constants
    '(T/m)_{mot}': 70 * units('N*m/kg'),  # 70 for superconducting; 20 for traditional
    # ducted fan constants
    #'K_{fan}': 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
    'K_{fan}': (5216/2) * units('lb') / (61 * units('in'))**3, # half of CFM-56 weight is turbine and half is fan
    'HTR': 0.38, 
    #'d_{fan}': 56*units('in'),
    'd_{fan}': 1.2*units('m'), # CHEETA fan diameter
    # fuel cell constants
    'P_{max,fc}':               850*units('kW'), # max power the model is accurate for. This could be replaced with a max power specified by Boeing
    'm_{comp}':                 36*units('kg'),
    'N_s':                      6,
    'N_{eng}':                  2,
    # HEX constants
    '\psi':                    5.02*units('kg/ft^2'), # from Boeing (Chellappa)
})

# set up aircraft system (to which propulsion system is passed as input)
ac = Aircraft(he)
ac.substitutions.update({
    'f_{empty}': 0.554  , # tuned to Elias' weights
    
    'b_{fuse}': 14*units('ft'),
    'b_{wing}': (136-14)*units('ft'),
    'S_{ref}': 2000*units('ft^2'),
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
    'm_{pay}':  15875* units('kg'), # design payload
})

mission.substitutions.update({
    'L/D':              LDR, # take input lift-to-drag ratios from Elias
    #'f_{BLI,max}':                   0.3388, 
    '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]; according to D.K. Hall, this number is generally representative         
    #'D_p':                       4830*units('lbf'), # profile drag ingested by propulsors
    
    '\\eta_{th}':                0.55,
    
    # fan constants
   # 'D_corr':           0.5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
   # '\\eta_{fan}':               0.95, 
   # 'Cd':                        0.04, # [Robinson, 2017]
        '\phi_{fuse}':               .5, # assume half of fuselage dissipation occurs on top of fuselage 
        '\phi_{wing}':               .64, # assume most dissipation occurs on top of wing
    
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
    
    # generator constants
    'omega_{gen}': 16000 * units('rpm'),
    '1-eta_{gen}': 1 - .98,
    
    # motor constants
    'omega_{motor}': 4000 * units('rpm'),
    '1-eta_{motor}': 1 - .999, #.999 for superconducting; .95 for traditional
    })

init = {
    mission['(i0/A)']:                                              1e-3*units('A/cm**2'),
    mission['P_{comp}']:                                            1e12*units('kW'),
    mission['V_{act}']:                                             .3*units('V'),
    mission['T_{t\\2}']:                                            (180+273.15)*units('K'),
    mission['\\dot{m}_{a}']:                                        360*units('kg/s'),
    mission['eta_{comp}']:                                          .9,
    mission['F_{net}']:                                             44*units('kN'),
    mission.fs.performance.prop.fan_fuse_perf['u_{jet}']:           10*units('m/s'),
    mission.fs.performance.prop.fan_wing_perf['u_{jet}']:           10*units('m/s'),
    mission.fs.performance.prop.fan_fuse_perf['\\dot{m}']:          1000000*units('kg/s'),
    mission.fs.performance.prop.fan_wing_perf['\\dot{m}']:          1000000*units('kg/s'),
    mission['f_{fuse,BLI}']:                                        .1,
    mission['f_{wing,BLI}']:                                        .01,
    mission['pi_c']:                                                5,
    mission['F_{HEX}']:                                             1*units('kN'),
    mission['D_{core}']:                                            1*units('kN'),
    mission['P_{gross}']:                                           200*units('kW'),
    mission['dt']:                                                  1000*units('s'),
    mission.fs.performance.prop.fan_fuse_perf['F_{fan}']:           10*units('kN'),
    mission.fs.performance.prop.fan_wing_perf['F_{fan}']:           10*units('kN'),
    mission.fs.performance.prop.gen_perf['P_{out}']:                200*units('kW'),
    mission.fs.performance.prop.motor_perf['P_{out}']:              200*units('kW'),
    mission.fs.performance.prop.stack_perf['P_{out}']:              200*units('kW'),
    mission.fs.performance.prop.core_perf['P_{out}']:               200*units('kW'),
    mission['d_{fan}']:                                             56*units('in'),
    mission['N_s']:                                                 1e-12,    
    }

model = Model(mission['PFEI'], [ac, mission])

sol = model.localsolve(x0=init,solver='mosek_conif',verbosity=2, iteration_limit=1000)
print(sol.table())
print(sol.summary())

print("PFEI (kJ/kg/km): " + str(sol['variables']['PFEI']))
print("MTOW (kg): " + str(sol['variables']['MTOW']))
print("OEW (kg): " + str(sol['variables']['m_{OE}']))
print("Fuel Weight (kg): " + str(sol['variables']['m_{fuel,tot}']))
print("Fuel Cell Weight (kg): " + str(sol['variables']['W_{stack}']))
print("d_{fan} (in): " + str(sol['variables']['d_{fan}']*39.37))
print("f_{wing,BLI}: " + str(sol['variables']['f_{wing,BLI}']))
print("f_{fuse,BLI}: " + str(sol['variables']['f_{fuse,BLI}']))
print(sol['variables']['\\dot{m}'])

print("# of stacks: " + str(sol['variables']['N_s']))
