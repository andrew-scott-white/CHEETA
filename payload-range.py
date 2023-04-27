import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from gpkit import Model, units

from flight_state import FlightState

from aircraft import Aircraft
from mission_OptRange import Mission

from turbofan import Turbofan
from turboelectric import Turboelectric
from hybrid_electric import HybridElectric
from fully_electric import FullyElectric

# GLOBAL CONSTANTS
gamma = 1.4
g = 9.81*units('m/s/s')
R2 = 8.314*units('J/K')
F = 96485.33*units('C')


# PARAMETER SWEEPS

# Number of engines
N_eng = 2

# Number of motor-driven propulsors
N_p = 9

# design payload weight
m_pay = 15875
#m_pay = np.linspace(500,30000,num=20)

# design range
R_des = np.linspace(8000,1000, num=60) # used a small number of iterations due to the instability of hybrid model. Too many signomials likely the cause.

# PEMFC stack specific power 
P_m = 2700 # W/kg

    
#   PROPULSION SYSTEMS
tf = Turbofan(N_eng)
tf.substitutions.update({
    # gas generator constants
    '(P/m)_{gg}': 7687 * units('hp') / 3065 / units('lb'), # half of LEAP-1B
    # ducted fan constants
    'K_{fan}': 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
    'HTR': 0.3,
    })

te_non = Turboelectric(N_eng, N_p) # non-superconducting turboelectric
te_non.substitutions.update({
    # gas generator constants
    '(P/m)_{gg}': 7687 * units('hp') / 3065 / units('lb'), # half of LEAP-1B
    # generator constants
    '(T/m)_{gen}': 20 * units('N*m/kg'),
    # motor constants
    '(T/m)_{mot}': 20 * units('N*m/kg'),  # 70 for superconducting; 20 for traditional
    # ducted fan constants
    'K_{fan}': 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
    'HTR': 0.38, 
    })


fe = FullyElectric(N_p, P_m)
fe.substitutions.update({
    # ducted fan constants
    'K_{fan}':                 0.5 * 6130 * units('lb') / ((69.4 * units('in'))**3), # half of LEAP-1B
    'HTR':                     0.38,
    #'d_{fan}':                 56*units('in'), # current CHEETA result
    
     # motor constants
    '(T/m)_{mot}':             70*units('N*m/kg'),  # anticipated value for superconducting machines [] | CHEETA table 47.7465
    
    # fuel cell constants
    'N_{cells}':                1080, # number of cells in a stack
    'A':                        446*units('cm**2'), # membrane area
    'P_{max,fc}':               850*units('kW'), 
    'm_{comp}':                 36*units('kg'),
    #'N_s':                      29,
    
    # HEX constants
    '\psi':                    5.02*units('kg/ft^2'), # from Boeing (Chellappa)
})
    
te_sc = Turboelectric(N_eng, N_p)
te_sc.substitutions.update({
    # gas generator constants
    '(P/m)_{gg}': 7687 * units('hp') / 3065 / units('lb'), # LEAP-1B
    # generator constants
    '(T/m)_{gen}': 20 * units('N*m/kg'),
    # motor constants
    '(T/m)_{mot}': 70 * units('N*m/kg'),  # 70 for superconducting; 20 for traditional
    # ducted fan constants
    'K_{fan}': 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
    'HTR': 0.38, 
    })


he = HybridElectric(N_p, P_m)
he.substitutions.update({
    
    # ducted fan constants
    'K_{fan}':                 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
    'HTR':                     0.38,
     # gas generator constants
    '(P/m)_{gg}':               7687 * units('hp') / 3065 / units('lb'), # LEAP-1B
    # generator constants
    '(T/m)_{gen}':              20 * units('N*m/kg'),
     # motor constants
    '(T/m)_{mot}':             70*units('N*m/kg'),  # anticipated value for superconducting machines []
    # fuel cell constants
    #'(P/m)_{fc}':              2700*units('W/kg'), # CHEETA Data Sheet
    'P_{max,fc}':               850*units('kW'), # max power the model is accurate for. This could be replaced with a max power specified by Boeing
    'm_{comp}':                 36*units('kg'),
    'N_s':                      18,
    'N_{eng}':                  2,
    # HEX constants
    '\psi':                    5.02*units('kg/ft^2'), # from Boeing (Chellappa)
    })







#   AIRCRAFT SYSTEMS
ac_737 = Aircraft(tf)
ac_737.substitutions.update({
    'f_{empty}': 0.5,  # sans propulsion system
    })

ac_737_te = Aircraft(te_non)
ac_737_te.substitutions.update({
    'f_{empty}': 0.5,  # sans propulsion system
    })

ac_ZeroE = Aircraft(tf)
ac_ZeroE.substitutions.update({
    'f_{empty}': 0.54,  # sans propulsion system
    })

ac_CHEETA = Aircraft(fe)
ac_CHEETA.substitutions.update({
    'f_{empty}': 0.54, # sans propulsion system
    })

ac_CHEETA_te = Aircraft(te_sc)
ac_CHEETA_te.substitutions.update({
    'f_{empty}': 0.54, # sans propulsion system
    })

ac_CHEETA_he = Aircraft(he)
ac_CHEETA_he.substitutions.update({
    'f_{empty}': 0.54, # sans propulsion system
    })





# FLIGHT STATE PARAMETERS
flight_state = FlightState()
flight_state.substitutions.update({ 
    # 0 kft
    'P_{t\\infty,TO}':       101325 * (1 + 0.5*(gamma-1)*0.25**2)**(gamma/(gamma-1)) * units('Pa'),
    'T_{t\\infty,TO}':       311.15 * (1 + 0.5*(gamma-1)*0.25**2) * units('K'),
    'V_\\infty,TO':          88.4*units('m/s'),
    'M_\infty,TO':           0.25,
    'rho_\infty,TO':         1.135*units('kg/m**3'),
    
    # 37.5 kft
    'M_\infty':           0.773, 
    'P_{t\\infty}':       21139 * (1 + 0.5*(gamma-1)*0.773**2)**(gamma/(gamma-1)) * units('Pa'),
    'T_{t\\infty}':       217.786 * (1 + 0.5*(gamma-1)*0.773**2) * units('K'),
    'V_\\infty':          443.369 * units('kts'),

    'rho_\infty':         .3509518*units('kg/m**3'),
    'C_p_\infty':         1006*units('J/kg/K'),
    '\\mu_\infty':        0.00001447704*units('Pa*s'),
})

# CONVENTIONAL 737 MISSION
mission_737 = Mission(ac_737, flight_state, R_des, m_pay)
# all aircraft and propulsion system performance models (and thus their
# variables) get set up inside Mission(), so we update their constants here

mission_737.substitutions.update({
    'F_{net,max}':               125*units('kN'), # thrust required at take-off
    '\rho_{fuel}':              840*units('kg/m^3'),
    'f_{TO}':                   0.10, # fuel fraction required for take-off/climb
    
    # Boundary-layer ingestion (BLI) parameters
       'f_{BLI}':                   0.17, # from D8, [Hall et al, 2017]
       '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
       'D_p':                       4830*units('lbf'),
    
    # gas generator performance constants
        'HV_{gg}': 43 * units('MJ/kg'),
        '\\eta_{th}': .55,
              
    # ducted fan performance constant
        '\\eta_{fan}': 0.95, # [Chandel, 2021]
        'Cd': .05,
        'D_{cruise}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{cruise}': .4870425385, # D(M_fan=0.6) corrected flow per unit area
        #'D_{cruise}': .4319187993, # D(M_fan=0.5) corrected flow per unit area
        'D_{takeoff}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{takeoff}': .2408545822, # D(M_fan=0.25) corrected flow per unit area

    # aircraft performance constant
        'L/D': 18.4,
    # mission parameters
        'f_{res}': 1.15,
    # flight state parameters
        'T_{t\\2}': 242.541*units('K'),
})

model = Model(mission_737['m_{fuel}'], [ac_737, mission_737])

sol_737 = model.localsolve(skipsweepfailures=True)
#print(sol_737.table())
#print(sol_737.summary())    

MTOW_737 = sol_737['variables']['MTOW']
dryWeight_737 = sol_737['variables']['m_{ZF}']
payload_737 = sol_737['variables']['m_{pay}']
R_des_737 = sol_737['variables']['R_{des}']
PFEI_737 = sol_737['variables']['PFEI']
V_fuel_737 = sol_737['variables']['V_{fuel}']



# TURBOELECTRIC 737 MISSION
mission_737_te = Mission(ac_737_te, flight_state, R_des, m_pay)
# all aircraft and propulsion system performance models (and thus their
# variables) get set up inside Mission(), so we update their constants here
mission_737_te.substitutions.update({
    'F_{net,max}':               125*units('kN'), # thrust required at take-off
    '\rho_{fuel}':              840*units('kg/m^3'),
    'f_{TO}':                   0.10, # fuel fraction required for take-off/climb
    
        # Boundary-layer ingestion (BLI) parameters
       'f_{BLI,max}':               0.3388, # weighted average of wing and fuselage BLI based on profile drag values from Elias 
       '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
       'D_p':                       4830*units('lbf'), # total profile drag for Cd_wing = .005784, Cd_fuse = .00804 for Sref = 1840 ft^2
    
    # turboelectric substitutions
        # generator constants
        'omega_{gen}': 22000 * units('rpm'),
        '1-eta_{gen}': 1 - .98,
        # motor constants
        'omega_{motor}': 4000 * units('rpm'),
        '1-eta_{motor}': 1 - .95, #.999 for superconducting; .95 for traditional
       })

mission_737_te.substitutions.update({
    # gas generator performance constants
        'HV_{gg}': 43 * units('MJ/kg'),
        '\\eta_{th}': 0.55,
       
    # ducted fan performance constant
        '\\eta_{fan}': 0.95, # [Chandel, 2021]
        'Cd': .05,
        'D_{cruise}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{cruise}': .4870425385, # D(M_fan=0.6) corrected flow per unit area
        #'D_{cruise}': .4319187993, # D(M_fan=0.5) corrected flow per unit area
        'D_{takeoff}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{takeoff}': .2408545822, # D(M_fan=0.25) corrected flow per unit area

    # aircraft performance constant
        'L/D': 18.4, # 16.8475292 for CHEETA, 18.4 for turboE
    # mission parameters
        'f_{res}': 1.15,
})

model = Model(mission_737_te['m_{fuel}'], [ac_737_te, mission_737_te])

sol_737_te = model.localsolve(skipsweepfailures=True)
print(sol_737_te.table())

MTOW_737_te = sol_737_te['variables']['MTOW']
dryWeight_737_te = sol_737_te['variables']['m_{ZF}']
payload_737_te = sol_737_te['variables']['m_{pay}']
R_des_737_te = sol_737_te['variables']['R_{des}']
PFEI_737_te = sol_737_te['variables']['PFEI']
V_fuel_737_te = sol_737_te['variables']['V_{fuel}']



# H2-BURNING TURBOFAN (A LA AIRBUS ZERO-E)
mission_ZeroE = Mission(ac_ZeroE, flight_state, R_des, m_pay)
# all aircraft and propulsion system performance models (and thus their
# variables) get set up inside Mission(), so we update their constants here

mission_ZeroE.substitutions.update({
    'F_{net,max}':               125*units('kN'), # thrust required at take-off
    '\rho_{fuel}':              71*units('kg/m^3'),
    'f_{TO}':                   0.10, # fuel fraction required for take-off/climb
    
    # Boundary-layer ingestion (BLI) parameters
       'f_{BLI}':                   0.17, # from D8, [Hall et al, 2017]
       '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
       'D_p':                       4830*units('lbf'), # from CHEETA, [Waddington, 2021]
    
    # gas generator performance constants
        'HV_{gg}': 120 * units('MJ/kg'),
        '\\eta_{th}': .55,
              
    # ducted fan performance constant
        '\\eta_{fan}': 0.95, # [Chandel, 2021]
        'Cd': .05,
        #'D_{cruise}': .4870425385, # D(M_fan=0.6) corrected flow per unit area
        'D_{cruise}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{cruise}': .4870425385, # D(M_fan=0.6) corrected flow per unit area
        #'D_{cruise}': .4319187993, # D(M_fan=0.5) corrected flow per unit area
        'D_{takeoff}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{takeoff}': .2408545822, # D(M_fan=0.25) corrected flow per unit area

    # aircraft performance constant
        'L/D': 18.4,
    # mission parameters
        'f_{res}': 1.15,
    # flight state parameters
        'T_{t\\2}': 242.541*units('K'),
})

model = Model(mission_ZeroE['m_{fuel}'], [ac_ZeroE, mission_ZeroE])

sol_ZeroE = model.localsolve(skipsweepfailures=True)
#print(sol_ZeroE.table())
 
MTOW_ZeroE = sol_ZeroE['variables']['MTOW']
dryWeight_ZeroE = sol_ZeroE['variables']['m_{ZF}']
payload_ZeroE = sol_ZeroE['variables']['m_{pay}']
R_des_ZeroE = sol_ZeroE['variables']['R_{des}']
ZeroE_PFEI = sol_ZeroE['variables']['PFEI']
V_fuel_ZeroE = sol_ZeroE['variables']['V_{fuel}']


# FULLY-ELECTRIC H2 PEMFC (A LA CHEETA)
mission_CHEETA = Mission(ac_CHEETA, flight_state, R_des, m_pay)    
mission_CHEETA.substitutions.update({
    'F_{net,max}':               125*units('kN'), # thrust required at take-off
    '\rho_{fuel}':              71*units('kg/m^3'),
    'f_{TO}':                   0.10, # fuel fraction required for take-off/climb
    
    # Boundary-layer ingestion (BLI) parameters
       'f_{BLI,max}':               0.3388, # weighted average of wing and fuselage BLI based on profile drag values from Elias 
       '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
       'D_p':                       4830*units('lbf'), # total profile drag for Cd_wing = .005784, Cd_fuse = .00804 for Sref = 1840 ft^2
        
    # high-temperature parameters            
        'HV_{fc}':              120.0*units('MJ/kg'), # fuel lower heating value
        'HV_{gg}':              120.0*units('MJ/kg'), # (no GG, but needed to satisfy PFEI constraint in mission)
        'T_{t\\2}':             (180+273.15)*units('K'), # stack operating temperature
        '\\DeltaG':             224.24*units('kJ'),
        'k':                    5e-5*units('1/s'),   # electrochemical rate constant
        '\\ohm_n':              .062*units('ohm*cm**2'), # fuel cell internal resistance
        'x_{parasitic}':        .05,        
        'm':                    4e-5*units('V'), # Elias
        'n':                    2.3e-3*units('cm**2/mA'), # Elias
        
    # fuel cell design inputs
        'P_{t\\fuel}':          5*units('bar'), # fuel storage pressure
        '\mu_f':                .95, # fuel utilization factor
        'alpha':                .5, # charge transfer coefficient
        'c':                    .21, # fraction of O2 in air
        'S':                    1.5, # stoichiometry (surplus air)
        'eta_{comp}':           .9, # fuel cell compressor efficiency
        
    # HEX parameters
        # constants
        'Pr':                   0.72,
        'C_d':                  .05,
        
        # HEX design
        'K_f':                  2,
        'K_h':                  0.93,
        '\sigma':               .464,
        '\epsilon':             0.85, # Kays and London
        '\epsilon_{to}':        0.60, # Kays and London
                  
    # motor constants
        'omega_{motor}':        4000 * units('rpm'),
        '1-eta_{motor}':        1 - .999, # reduced efficiency to accomodate cooling. Was 95% but could be 99.9% if we model cooling power draw
        
    # ducted fan performance constant
        '\\eta_{fan}': 0.95, # [Chandel, 2021]
        'Cd':                   0.05,
        'D_{cruise}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{cruise}': .4870425385, # D(M_fan=0.6) corrected flow per unit area
        #'D_{cruise}': .4319187993, # D(M_fan=0.5) corrected flow per unit area
        'D_{takeoff}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{takeoff}': .2408545822, # D(M_fan=0.25) corrected flow per unit area
        
    # aircraft performance constant
        'L/D':                  18.4, # 18.4 w/out HEX | from Elias
        
    # mission parameters
        'f_{res}':              1.15,
})

init = {
        mission_CHEETA['(i0/A)']:                        1e-3*units('A/cm**2'),
        mission_CHEETA['(i0/A,max)']:                    1e-3*units('A/cm**2'),
        mission_CHEETA['P_{comp}']:                      1e12*units('kW'),
        mission_CHEETA['P_{comp,max}']:                  1e12*units('kW'),
        mission_CHEETA['V_{act}']:                       .3*units('V'),
        mission_CHEETA['V_{act,max}']:                   .3*units('V'),
        mission_CHEETA['T_{t\\2}']:                      (180+273.15)*units('K'),
        mission_CHEETA['\\dot{m}_{a}']:                   360*units('kg/s'),
        mission_CHEETA['\\dot{m}_{a,max}']:               360*units('kg/s'),
        mission_CHEETA['eta_{comp}']:                    .9,
        #mission_CHEETA['Q_{HEX}']:                      300*units('W'),
        mission_CHEETA['F_{net}']:                       44*units('kN'),
        mission_CHEETA['u_{jet}']:                       100*units('m/s'),
        }

model = Model(mission_CHEETA['m_{fuel}'], [ac_CHEETA, mission_CHEETA])  

sol_CHEETA = model.localsolve(x0=init, solver='mosek_conif', verbosity=2, skipsweepfailures=True, pccp_penalty=200)
print(sol_CHEETA.table())

MTOW_CHEETA = sol_CHEETA['variables']['MTOW']
dryWeight_CHEETA = sol_CHEETA['variables']['m_{ZF}']
payload_CHEETA = sol_CHEETA['variables']['m_{pay}']
R_des_CHEETA = sol_CHEETA['variables']['R_{des}']
CHEETA_PFEI = sol_CHEETA['variables']['PFEI']
V_fuel_CHEETA = sol_CHEETA['variables']['V_{fuel}']

    


# CHEETA TURBO-ELECTRIC
mission_CHEETA_te = Mission(ac_CHEETA_te, flight_state, R_des, m_pay)
# all aircraft and propulsion system performance models (and thus their
# variables) get set up inside Mission(), so we update their constants here
mission_CHEETA_te.substitutions.update({
    'F_{net,max}':               125*units('kN'), # thrust required at take-off
    '\rho_{fuel}':              71*units('kg/m^3'),
    'f_{TO}':                   0.10, # fuel fraction required for take-off/climb
    
        # Boundary-layer ingestion (BLI) parameters
       'f_{BLI,max}':               0.3388, # weighted average of wing and fuselage BLI based on profile drag values from Elias 
       '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
       'D_p':                       4830*units('lbf'), # total profile drag for Cd_wing = .005784, Cd_fuse = .00804 for Sref = 1840 ft^2

    
    # turboelectric substitutions
        # generator constants
        'omega_{gen}':  22000 * units('rpm'),
        '1-eta_{gen}':  1 - .98,
        # motor constants
        'omega_{motor}': 4000 * units('rpm'),
        '1-eta_{motor}': 1 - .999, #.999 for superconducting; .95 for traditional
       })

mission_CHEETA_te.substitutions.update({
    # gas generator performance constants
        'HV_{gg}':      120 * units('MJ/kg'),
        '\\eta_{th}':   0.55,
       
    # ducted fan performance constant
        '\\eta_{fan}': 0.95, # [Chandel, 2021]
        'Cd': .05,
        'D_{cruise}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{cruise}': .4870425385, # D(M_fan=0.6) corrected flow per unit area
        #'D_{cruise}': .4319187993, # D(M_fan=0.5) corrected flow per unit area
        'D_{takeoff}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{takeoff}': .2408545822, # D(M_fan=0.25) corrected flow per unit area
        
    # aircraft performance constant
        'L/D':          18.4, # 16.8475292 for CHEETA, 18.4 for turboE
    # mission parameters
        'f_{res}':      1.15,

})

model = Model(mission_CHEETA_te['m_{fuel}'], [ac_CHEETA_te, mission_CHEETA_te])

sol_CHEETA_te = model.localsolve(skipsweepfailures=True)
#print(sol_CHEETA_te.table())

MTOW_CHEETA_te = sol_CHEETA_te['variables']['MTOW']
dryWeight_CHEETA_te = sol_CHEETA_te['variables']['m_{ZF}']
payload_CHEETA_te = sol_CHEETA_te['variables']['m_{pay}']
R_des_CHEETA_te = sol_CHEETA_te['variables']['R_{des}']
CHEETA_te_PFEI = sol_CHEETA_te['variables']['PFEI']
V_fuel_CHEETA_te = sol_CHEETA_te['variables']['V_{fuel}']




# CHEETA HYBRID-ELECTRIC
mission_CHEETA_he = Mission(ac_CHEETA_he, flight_state, R_des, m_pay)
    # all aircraft and propulsion system performance models (and thus their
    # variables) get set up inside Mission(), so we update their constants here
    
mission_CHEETA_he.substitutions.update({
        'F_{net,max}':               125*units('kN'), # thrust required at take-off
        '\rho_{fuel}':                  71*units('kg/m^3'),
        'f_{TO}':                   0.10, # fuel fraction required for take-off/climb
        
        # Boundary-layer ingestion (BLI) parameters
       'f_{BLI,max}':               0.3388, # weighted average of wing and fuselage BLI based on profile drag values from Elias 
       '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
       'D_p':                       4830*units('lbf'), # total profile drag for Cd_wing = .005784, Cd_fuse = .00804 for Sref = 1840 ft^2
        
    # gas generator performance constants
        'HV_{gg}':                  120 * units('MJ/kg'),
        '\\eta_{th}':               .55,
        
    # high-temperature parameters            
        'T_{t\\2}':             (180+273.15)*units('K'), # stack operating temperature
        '\\DeltaG':             224.24*units('kJ'),
        'k':                    5e-5*units('1/s'),   # electrochemical rate constant
        '\\ohm_n':              .062*units('ohm*cm**2'), # fuel cell internal resistance
        'x_{parasitic}':        .05,
        'm':                    4e-5*units('V'), # Elias
        'n':                    2.3e-3*units('cm**2/mA'), # Elias

        
    # fuel cell design inputs
        'HV_{fc}':              120 * units('MJ/kg'),
        'P_{t\\fuel}':          5*units('bar'), # fuel storage pressure
        'A':                    446*units('cm**2'), # membrane area
        'N_{cells}':            1080, # number of cells in a stack
        '\mu_f':                .95, # fuel utilization
        'alpha':                .5, # charge transfer coefficient (empirical)
        'c':                    .21, # fraction of O2 in air
        'S':                    1.5, # stoichiometry (surplus air)
        'eta_{comp}':           .9, # fuel cell compressor efficiency
        
    # HEX parameters
        # assumptions
        'Pr':                   0.72,
        'C_d':                  .05,

        # operating conditions (for matrix) 
        'K_f':                  2,
        'K_h':                  0.93,
        '\sigma':               .464,
        '\epsilon':             0.85, # Kays and London
        '\epsilon_{to}':        0.60, # Kays and London
        
    # generator constants
        'omega_{gen}':          22000 * units('rpm'),
        '1-eta_{gen}':          1 - .98,
        
    # motor constants
        'omega_{motor}':        4000 * units('rpm'),
        '1-eta_{motor}':        1 - .999, # reduced efficiency to accomodate cooling. Was 95% but could be 99.9% if we model cooling power draw
        
    # ducted fan performance constant
        '\\eta_{fan}': 0.95, # [Chandel, 2021]
        'Cd':                   0.05,
        'D_{cruise}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{cruise}': .4870425385, # D(M_fan=0.6) corrected flow per unit area
        #'D_{cruise}': .4319187993, # D(M_fan=0.5) corrected flow per unit area
        'D_{takeoff}': .5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'D_{takeoff}': .2408545822, # D(M_fan=0.25) corrected flow per unit area
        
    # aircraft performance constant
        'L/D':                  18.4, # 18.4 w/out HEX | from Elias
        
    # mission parameters
        'f_{res}':              1.15
})

init = {
    mission_CHEETA_he['(i0/A)']:                        1e-3*units('A/cm**2'),
    mission_CHEETA_he['(i0/A,max)']:                    1e-3*units('A/cm**2'),
    mission_CHEETA_he['P_{comp}']:                      1e12*units('kW'),
    mission_CHEETA_he['P_{comp,max}']:                  1e12*units('kW'),
    mission_CHEETA_he['V_{act}']:                       .3*units('V'),
    mission_CHEETA_he['V_{act,max}']:                   .3*units('V'),
    mission_CHEETA_he['T_{t\\2}']:                      (180+273.15)*units('K'),
    mission_CHEETA_he['\\dot{m}_{a}']:                   360*units('kg/s'),
    mission_CHEETA_he['\\dot{m}_{a,max}']:               360*units('kg/s'),
    mission_CHEETA_he['eta_{comp}']:                    .9,
    #mission_CHEETA_he['Q_{HEX}']:                      300*units('W'),
    mission_CHEETA_he['F_{net}']:                       44*units('kN'),
    mission_CHEETA_he['u_{jet}']:                       100*units('m/s'),
    mission_CHEETA_he['P_{max,gg}']:                    8*units('MW'),
    #mission_CHEETA_he['F_{HEX,to}']:              1e-12*units('N'),
    }



model = Model(mission_CHEETA_he['m_{fuel}'], [ac_CHEETA_he, mission_CHEETA_he])  
#model = Model(mission_CHEETA_he['m_{\\rm fuel}'], Bounded([ac, mission]))

sol_CHEETA_he = model.localsolve(x0=init, solver='mosek_conif', verbosity=2, iteration_limit=1000, skipsweepfailures=True)
print(sol_CHEETA_he.table())

MTOW_CHEETA_he = sol_CHEETA_he['variables']['MTOW']
dryWeight_CHEETA_he = sol_CHEETA_he['variables']['m_{ZF}']
payload_CHEETA_he = sol_CHEETA_he['variables']['m_{pay}']
R_des_CHEETA_he = sol_CHEETA_he['variables']['R_{des}']
CHEETA_he_PFEI = sol_CHEETA_he['variables']['PFEI']
V_fuel_CHEETA_he = sol_CHEETA_he['variables']['V_{fuel}']


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 11}

matplotlib.rc('font', **font)

        
# PLOT PAYLOAD-RANGE DIAGRAMS
fig1 = plt.figure(1)
line11, = plt.plot(R_des_737, dryWeight_737, 'k:', label='ATF Turbo-fan')
line21, = plt.plot(R_des_737_te, dryWeight_737_te, 'k--', label='ATF Turbo-electric')
line31, = plt.plot(R_des_ZeroE, dryWeight_ZeroE, 'b:', label='Hydrogen Turbo-fan')
line41, = plt.plot(R_des_CHEETA_te, dryWeight_CHEETA_te, 'b--', label='Hydrogen Turbo-electric')
line51, = plt.plot(R_des_CHEETA, dryWeight_CHEETA, 'b', label='Hydrogen Fully-electric')
line61, = plt.plot(R_des_CHEETA_he, dryWeight_CHEETA_he, 'b-.', label='Hydrogen Hybrid-electric')
plt.xlabel("Range (km)")
plt.ylabel("Zero-fuel Weight (kg)")
plt.ylim([0, 100000])
plt.xlim([1000, 8000])
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend(handles=[line11, line21, line31, line41, line51, line61])

fig2 = plt.figure(2)
line12, = plt.plot(R_des_737, PFEI_737, 'k:', label='ATF Turbo-fan')
line22, = plt.plot(R_des_737_te, PFEI_737_te, 'k--', label='ATF Turbo-electric')
line32, = plt.plot(R_des_ZeroE, ZeroE_PFEI, 'b:', label='Hydrogen Turbo-fan')
line42, = plt.plot(R_des_CHEETA_te, CHEETA_te_PFEI, 'b--', label='Hydrogen Turbo-electric')
line52, = plt.plot(R_des_CHEETA, CHEETA_PFEI, 'b', label='Hydrogen Fully-electric')
line62, = plt.plot(R_des_CHEETA_he, CHEETA_he_PFEI, 'b-.', label='Hydrogen Hybrid-electric')
plt.xlabel("Range (km)")
plt.ylabel("PFEI (kJ/kg-km)")
plt.ylim([0, 10])
plt.xlim([1000, 8000])
plt.legend(handles=[line12, line22, line32, line42, line52, line62])

fig2 = plt.figure(3)
line13, = plt.plot(R_des_737, MTOW_737, 'k:', label='ATF Turbo-fan')
line23, = plt.plot(R_des_737_te, MTOW_737_te, 'k--', label='ATF Turbo-electric')
line33, = plt.plot(R_des_ZeroE, MTOW_ZeroE, 'b:', label='Hydrogen Turbo-fan')
line43, = plt.plot(R_des_CHEETA_te, MTOW_CHEETA_te, 'b--', label='Hydrogen Turbo-electric')
line53, = plt.plot(R_des_CHEETA, MTOW_CHEETA, 'b', label='Hydrogen Fully-electric')
line63, = plt.plot(R_des_CHEETA_he, MTOW_CHEETA_he, 'b-.', label='Hydrogen Hybrid-electric')
plt.xlabel("Range (km)")
plt.ylabel("Max Take-off Weight (kg)")
plt.ylim([0, 125000])
plt.xlim([1000, 8000])
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend(handles=[line12, line22, line32, line42, line52, line62])

fig2 = plt.figure(4)
line13, = plt.plot(R_des_737, V_fuel_737, 'k:', label='ATF Turbo-fan')
line23, = plt.plot(R_des_737_te, V_fuel_737_te, 'k--', label='ATF Turbo-electric')
line33, = plt.plot(R_des_ZeroE, V_fuel_ZeroE, 'b:', label='Hydrogen Turbo-fan')
line43, = plt.plot(R_des_CHEETA_te, V_fuel_CHEETA_te, 'b--', label='Hydrogen Turbo-electric')
line53, = plt.plot(R_des_CHEETA, V_fuel_CHEETA, 'b', label='Hydrogen Fully-electric')
line63, = plt.plot(R_des_CHEETA_he, V_fuel_CHEETA_he, 'b-.', label='Hydrogen Hybrid-electric')
plt.xlabel("Range (km)")
plt.ylabel("Required Fuel Volume ($m^3$)")
plt.ylim([0, 150])
plt.xlim([1000, 8000])
plt.legend(handles=[line12, line22, line32, line42, line52, line62])