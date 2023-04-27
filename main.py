import numpy as np
import matplotlib.pyplot as plt

from gpkit import Model, Variable, units, Vectorize

from flight_state import FlightState

from aircraft import Aircraft
from mission import Mission

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

# design range
R_des = 2935 # nmi
#R_des = [1500, 2000, 2500, 2935, 3500] # nmi

# design payload weight
m_pay = 35000 # lb
#m_pay = [20000 25000 30000 35000 40000 45000] # lb

# PEMFC stack specific power
#P_m = 2700 # W/kg
#P_m = [2000, 2700] # W/kg
P_m = [2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3000] # W/kg

#for x in R_des:
for ii in P_m:

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
        'HTR': 0.3, # .47 for superconducting, .3 for traditional
        })

    
    fe = FullyElectric(N_p, P_m)
    fe.substitutions.update({
        # ducted fan constants
        'K_{fan}':                 0.5 * 6130 * units('lb') / ((69.4 * units('in'))**3), # half of LEAP-1B
        'HTR':                     0.47,
        #'d_{fan}':                 56*units('in'), # current CHEETA result
        
         # motor constants
        '(T/m)_{mot}':             70*units('N*m/kg'),  # anticipated value for superconducting machines [] | CHEETA table 47.7465
        
        # fuel cell constants
        #'(P/m)_{fc}':               2700*units('W/kg'), # Hypoint
        'N_{cells}':                1080, # number of cells in a stack
        'A':                        446*units('cm**2'), # membrane area
        'P_{max,fc}':               850*units('kW'), 
        #'vol_{fc}':                 37.5*units('cm**3'), # volume of one cell [Wilson, 1996]
        'm_{comp}':                 36*units('kg'),
        #'N_s':                      29,
        
        # HEX constants
        '\psi':                    5.02*units('kg/ft^2'), # from Boeing (Chellappa)
       # 'A_{face}':                700*units('in**2'),
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
        'HTR': 0.47, 
        })
    
    
    he = HybridElectric(N_p, P_m)
    he.substitutions.update({
        
        # ducted fan constants
        'K_{fan}':                 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
        'HTR':                     0.47,
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
        'N_s':                      12,
        'N_{eng}':                  2,
        #'//nu_{conv}':              .8,
        # HEX constants
        '\psi':                    5.02*units('kg/ft^2'), # from Boeing (Chellappa)
        })

    
    
    
    
    
    
    #   AIRCRAFT SYSTEMS
    ac_737 = Aircraft(tf)
    ac_737.substitutions.update({
        'f_{empty}': 0.473, 
        'd_{fan}': 69.4*units('in'), # for 737 baseline with LEAP-1B twin engines
        })
    
    ac_737_te = Aircraft(te_non)
    ac_737_te.substitutions.update({
        'f_{empty}': 0.473, 
        })
    
    ac_ZeroE = Aircraft(tf)
    ac_ZeroE.substitutions.update({
            'f_{empty}': 0.53442, 
        })
    
    ac_CHEETA = Aircraft(fe)
    ac_CHEETA.substitutions.update({
        'f_{empty}':                0.53442, # no propulsion system mass (.518?? for Cheeta, .473 for 737)
        })

    ac_CHEETA_te = Aircraft(te_sc)
    ac_CHEETA_te.substitutions.update({
        'f_{empty}': 0.53442, # no propulsion system mass (.511 for Cheeta, .473 for 737)
        })
    
    ac_CHEETA_he = Aircraft(he)
    ac_CHEETA_he.substitutions.update({
        'f_{empty}':                0.53442, # no propulsion system mass (.511 for Cheeta, .473 for 737)
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
        'F_{net,max}':               23300*units('lbf'), # thrust required per propulsor at take-off
        
        # Boundary-layer ingestion (BLI) parameters
           'f_{BLI}':                   0.000000001, 
           '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
           #'f_{BLI}':                   0.0000001, # non-BLI
           'D_p':                       4830*units('lbf'),
        
        # gas generator performance constants
            'HV_{gg}': 43 * units('MJ/kg'),
            '\\eta_{th}': .45,
                  
        # ducted fan performance constant
            '\\eta_{fan}': 0.95, # 0.93 w/ BLI, 0.95 w/out BLI
            'Cd': .04,
            'D_{cruise}': .4870425385, # D(M_fan=0.6) corrected flow per unit area
            'D_{takeoff}': .2408545822, # D(M_fan=0.25) corrected flow per unit area
    
        # aircraft performance constant
            'L/D': 18.4,
        # mission parameters
            'f_{res}': 1.199921,
            #'R_{des}': 2935 * units('nmi'),
            #'m_{pay}':  35000 * units('lb'),
        # flight state parameters
            'T_{t\\2}': 242.541*units('K'),
    })
    
    model = Model(mission_737['m_{fuel}'], [ac_737, mission_737])
    
    sol_737 = model.localsolve()
    #print(sol_737.table())
    #print(sol_737.summary())    
    
    PFEI_737 = sol_737['variables']['PFEI']
    
    
    
    # TURBOELECTRIC 737 MISSION
    mission_737_te = Mission(ac_737_te, flight_state, R_des, m_pay)
    # all aircraft and propulsion system performance models (and thus their
    # variables) get set up inside Mission(), so we update their constants here
    mission_737_te.substitutions.update({
        'F_{net,max}':               2600*units('lbf'), # thrust required per propulsor at take-off   
        
            # Boundary-layer ingestion (BLI) parameters
           'f_{BLI,max}':               0.3388, # weighted average of wing and fuselage BLI based on profile drag values from Elias 
           '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
           'D_p':                       4830*units('lbf'), # total profile drag for Cd_wing = .005784, Cd_fuse = .00804 for Sref = 1840 ft^2
        
        # turboelectric substitutions
            # generator constants
            'omega_{gen}': 16000 * units('rpm'),
            '1-eta_{gen}': 1 - .98,
            # motor constants
            'omega_{motor}': 4000 * units('rpm'),
            '1-eta_{motor}': 1 - .95, #.999 for superconducting; .95 for traditional
           })
    
    mission_737_te.substitutions.update({
        # gas generator performance constants
            'HV_{gg}': 43 * units('MJ/kg'),
            '\\eta_{th}': 0.45,
           
        # ducted fan performance constant
            '\\eta_{fan}': 0.93, # 0.93 w/ BLI, 0.95 w/out BLI
            'Cd': .04,
            'D_{cruise}': .4319187993, # D(M_fan=0.5) corrected flow per unit area
            'D_{takeoff}': .2408545822, # D(M_fan=0.25) corrected flow per unit area
    
        # aircraft performance constant
            'L/D': 18.4, # 16.8475292 for CHEETA, 18.4 for turboE
        # mission parameters
            'f_{res}': 1.199921,
            #'R_{des}': 2935 * units('nmi'),
            #'m_{pay}': 35000 * units('lb'),
    })
    
    model = Model(mission_737_te['m_{fuel}'], [ac_737_te, mission_737_te])
    
    sol_737_te = model.localsolve()
    #print(sol_737_te.table())
    
    PFEI_737_te = sol_737_te['variables']['PFEI']
    
    
    
    
    # H2-BURNING TURBOFAN (A LA AIRBUS ZERO-E)
    mission_ZeroE = Mission(ac_ZeroE, flight_state, R_des, m_pay)
    # all aircraft and propulsion system performance models (and thus their
    # variables) get set up inside Mission(), so we update their constants here
    
    mission_ZeroE.substitutions.update({
        #'P_{des}':                 24*units('MW'),
        'F_{net,max}':               23300*units('lbf'), # thrust required per propulsor at take-off
        
        # Boundary-layer ingestion (BLI) parameters
           'f_{BLI}':                   0.000000001, 
           '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
           'D_p':                       4830*units('lbf'),
        
        # gas generator performance constants
            'HV_{gg}': 120 * units('MJ/kg'),
            '\\eta_{th}': .45,
                  
        # ducted fan performance constant
            '\\eta_{fan}': 0.95, # 0.93 w/ BLI, 0.95 w/out BLI
            'Cd': .04,
            'D_{cruise}': .4870425385, # D(M_fan=0.6) corrected flow per unit area
            'D_{takeoff}': .2408545822, # D(M_fan=0.25) corrected flow per unit area
    
        # aircraft performance constant
            'L/D': 18.4,
        # mission parameters
            'f_{res}': 1.199921,
            #'R_{des}': 2935 * units('nmi'),
            #'m_{pay}':  35000 * units('lb'),
        # flight state parameters
            'T_{t\\2}': 242.541*units('K'),
    })
    
    model = Model(mission_ZeroE['m_{fuel}'], [ac_ZeroE, mission_ZeroE])
    
    sol_ZeroE = model.localsolve()
    #print(sol_ZeroE.table())
    
    ZeroE_PFEI = sol_ZeroE['variables']['PFEI']
    
    
    
    # FULLY-ELECTRIC H2 PEMFC (A LA CHEETA)
    mission_CHEETA = Mission(ac_CHEETA, flight_state, R_des, m_pay)    
    mission_CHEETA.substitutions.update({
        'F_{net,max}':               2600*units('lbf'), # thrust required per propulsor at take-off
        
        # Boundary-layer ingestion (BLI) parameters
           'f_{BLI,max}':               0.3388, # weighted average of wing and fuselage BLI based on profile drag values from Elias 
           '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
           'D_p':                       4830*units('lbf'), # total profile drag for Cd_wing = .005784, Cd_fuse = .00804 for Sref = 1840 ft^2
            
        # high-temperature parameters            
            'HV_{fc}':              120.0*units('MJ/kg'), # fuel lower heating value
            'HV_{gg}':              120.0*units('MJ/kg'), # (no GG, but needed to satisfy PFEI constraint in mission)
            'T_{t\\2}':             (180+273.15)*units('K'), # stack operating temperature
            '\\DeltaG':             220.4*units('kJ'),
            'k':                    6e-6*units('1/s'),   # electrochemical rate constant
            '\\ohm_n':              .062*units('ohm*cm**2'), # fuel cell internal resistance
            'x_{parasitic}':        .05,        
            'm':                    8e-5*units('V'), # Elias
            'n':                    2e-4*units('cm**2/mA'), # Elias
            
        # fuel cell design inputs
            'P_{t\\fuel}':          5*units('bar'), # fuel storage pressure
            '\mu_f':                .96, # fuel utilization factor
            'alpha':                .5, # charge transfer coefficient
            'c':                    .21, # fraction of O2 in air
            'S':                    1.5, # stoichiometry (surplus air)
            'eta_{comp}':           .9, # fuel cell compressor efficiency
            
        # HEX parameters
            # constants
            'Pr':                   0.72,
            'C_d':                  .04,
            
            # HEX design
            'K_f':                  2,
            'K_h':                  0.93,
            '\sigma':               .464,
            '\epsilon':             0.60,
                      
        # motor constants
            'omega_{motor}':        4000 * units('rpm'),
            '1-eta_{motor}':        1 - .999, # reduced efficiency to accomodate cooling. Was 95% but could be 99.9% if we model cooling power draw
            
        # ducted fan performance constant
            '\\eta_{fan}':          0.93, # 0.93 w/ BLI, 0.95 w/out BLI
            'Cd':                   0.04,
            'D_{cruise}':           .4319187993, # D(M_fan=0.5) corrected flow per unit area
            'D_{takeoff}':          .2408545822, # D(M_fan=0.25) corrected flow per unit area
            
        # aircraft performance constant
            'L/D':                  18.4, # 18.4 w/out HEX | from Elias
            
        # mission parameters
            'f_{res}':              1.199921,
            #'R_{des}':              2935 * units('nmi'),
            #'m_{pay}':              35000 * units('lb'),
    })
    
    init = {
        mission_CHEETA['(i0/A)']:                   1e-3*units('A/cm**2'),
        mission_CHEETA['P_{comp}']:                 1e12*units('kW'),
        mission_CHEETA['V_{act}']:                  .3*units('V'),
        mission_CHEETA['T_{t\\2}']:                 (180+273.15)*units('K'),
        mission_CHEETA['eta_{fc}']:                 .6,
        mission_CHEETA['\\dot{m}_{a}']:             360*units('kg/s'),
        }
    
    model = Model(mission_CHEETA['m_{fuel}'], [ac_CHEETA, mission_CHEETA])  
    
    sol_CHEETA = model.localsolve(x0=init, solver='mosek_conif', verbosity=2)
    #print(sol_CHEETA.table())
    
    CHEETA_PFEI = sol_CHEETA['variables']['PFEI']
    
    
    
        
    
    
    # CHEETA TURBO-ELECTRIC
    mission_CHEETA_te = Mission(ac_CHEETA_te, flight_state, R_des, m_pay)
    # all aircraft and propulsion system performance models (and thus their
    # variables) get set up inside Mission(), so we update their constants here
    mission_CHEETA_te.substitutions.update({
        'F_{net,max}':               2600*units('lbf'), # thrust required per propulsor at take-off
        #'P_{des}':                 24*units('MW'),
        
            # Boundary-layer ingestion (BLI) parameters
           'f_{BLI,max}':               0.3388, # weighted average of wing and fuselage BLI based on profile drag values from Elias 
           '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
           #'f_{BLI}':                   0.0000001, # non-BLI
           'D_p':                       4830*units('lbf'), # total profile drag for Cd_wing = .005784, Cd_fuse = .00804 for Sref = 1840 ft^2
           #'D_{p,TO}':                  3086.564*units('lbf'), # total profile drag for Cd_wing = .004754, Cd_fuse = .008412 for Sref = 1840 ft^2
           #'l_{span}':                  135.65*units('ft'),
           #'\phi_{inj}':                0.9,
        
        # turboelectric substitutions
            # generator constants
            'omega_{gen}': 16000 * units('rpm'),
            '1-eta_{gen}': 1 - .98,
            # motor constants
            'omega_{motor}': 4000 * units('rpm'),
            '1-eta_{motor}': 1 - .999, #.999 for superconducting; .95 for traditional
           })
    
    mission_CHEETA_te.substitutions.update({
        # gas generator performance constants
            'HV_{gg}': 120 * units('MJ/kg'),
            '\\eta_{th}': 0.45,
           
        # ducted fan performance constant
            '\\eta_{fan}': 0.93, # 0.93 w/ BLI, 0.95 w/out BLI
            'Cd': .04,
            'D_{cruise}': .4319187993, # D(M_fan=0.5) corrected flow per unit area
            'D_{takeoff}': .2408545822, # D(M_fan=0.25) corrected flow per unit area
            
        # aircraft performance constant
            'L/D': 18.4, # 16.8475292 for CHEETA, 18.4 for turboE
        # mission parameters
            'f_{res}': 1.199921,
            #'R_{des}': 2935 * units('nmi'),
            #'m_{pay}': 35000 * units('lb'),
        # flight state parameters
            #'T_{t\\2}': 242.541*units('K'),
    })
    
    model = Model(mission_CHEETA_te['m_{fuel}'], [ac_CHEETA_te, mission_CHEETA_te])
    
    sol_CHEETA_te = model.localsolve()
    #print(sol_CHEETA_te.table())
    
    CHEETA_te_PFEI = sol_CHEETA_te['variables']['PFEI']
    
    
    
    
    
    
    
    # CHEETA HYBRID-ELECTRIC
    mission_CHEETA_he = Mission(ac_CHEETA_he, flight_state, R_des, m_pay)
        # all aircraft and propulsion system performance models (and thus their
        # variables) get set up inside Mission(), so we update their constants here
        
    mission_CHEETA_he.substitutions.update({
            'F_{net,max}':               2600*units('lbf'), # thrust required per propulsor at take-off
            #'P_{des}':                 24*units('MW'),
            
            # Boundary-layer ingestion (BLI) parameters
           'f_{BLI,max}':               0.3388, # weighted average of wing and fuselage BLI based on profile drag values from Elias 
           '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]
           #'f_{BLI}':                   0.0000001, # non-BLI
           'D_p':                       4830*units('lbf'), # total profile drag for Cd_wing = .005784, Cd_fuse = .00804 for Sref = 1840 ft^2
           #'D_{p,TO}':                  3086.564*units('lbf'), # total profile drag for Cd_wing = .004754, Cd_fuse = .008412 for Sref = 1840 ft^2
           #'l_{span}':                  135.65*units('ft'),
           #'\phi_{inj}':                0.9,
        
        # flight state parameters            
            #'rho_\infty':           .35762*units('kg/m**3'), # 37 kft
            #'C_p_\infty':           1003*units('J/kg/K'), 
            #'\\mu_\infty':          0.00001461*units('Pa*s'),
            
        # gas generator performance constants
            'HV_{gg}': 120 * units('MJ/kg'),
             #'h_{\\rm fuel}': 43 * units('MJ/kg'),
            '\\eta_{th}': .45,
            
        # low-temperature parameters
            #'T_{t\\2}':             (90+273.15)     *units('K'), # stack operating temperature
            #'\\Delta G':            228.2*units('kJ'),
            #'k':                    1.75e-7          *units('1/s'), # electrochemical rate constant
            #'\\ohm_n':              .038            *units('ohm*cm**2'), # fuel cell internal resistance
            #'x_{parasitic}':          .02,
            #'c_0':                  .015895*units('in**2/W'),
            
        # high-temperature parameters            
            'T_{t\\2}':             (180+273.15)*units('K'), # stack operating temperature
            '\\DeltaG':             220.4*units('kJ'),
            'k':                    6e-6*units('1/s'),   # electrochemical rate constant
            '\\ohm_n':              .062*units('ohm*cm**2'), # fuel cell internal resistance
            'x_{parasitic}':        .05,
            #'c_0':                  .007259*units('in**2/W'), # from Chellappa
            'm':                    8e-5*units('V'), # Elias
            'n':                    2e-4*units('cm**2/mA'), # Elias
            #'m':                    3e-5*units('V'), # Larminie
            #'n':                    8e-3*units('cm**2/mA'), # Larminie
            
        # fuel cell operating conditions (for matrix)
            #'P_{t\\2}':             .967*units('bar'), # air pressure at cathode
            #'e':                    1.2*units('V'),
            #'V':                    899.9*units('V'),
            #'i_{\\rm n}':           500*units('A'),
            
        # fuel cell design inputs
            'HV_{fc}':              120 * units('MJ/kg'),
            'P_{t\\fuel}':          5*units('bar'), # fuel storage pressure
            'A':                    446*units('cm**2'), # membrane area
            'N_{cells}':            1080, # number of cells in a stack
            '\mu_f':                .96, # fuel utilization
            'alpha':                .5, # charge transfer coefficient (empirical)
            'c':                    .21, # fraction of O2 in air
            'S':                    1.5, # stoichiometry (surplus air)
            #'pi_c':                 2.5,
            'eta_{comp}':           .9, # fuel cell compressor efficiency
            
        # HEX parameters
            # assumptions
            'Pr':                   0.72,
            'C_d':                  .04,
    
            # operating conditions (for matrix) 
            'K_f':                  2,
            'K_h':                  0.93,
            '\sigma':               .464,
            
             # Radiator exit temperature estimates
           # 'T_e': 360*units('K'), # conservative
            #'T_e': 396*units('K'), # optimistic
            
        # generator constants
            'omega_{gen}':          16000 * units('rpm'),
            '1-eta_{gen}':          1 - .98,
            
        # motor constants
            'omega_{motor}':        4000 * units('rpm'),
            '1-eta_{motor}':        1 - .999, # reduced efficiency to accomodate cooling. Was 95% but could be 99.9% if we model cooling power draw
            
        # ducted fan performance constant
            '\\eta_{fan}':          0.93, # 0.93 w/ BLI, 0.95 w/out BLI
            'Cd':                   0.04,
            'D_{cruise}':           .4319187993, # D(M_fan=0.6) corrected flow per unit area
            'D_{takeoff}':          .2408545822, # D(M_fan=0.25) corrected flow per unit area
            
        # aircraft performance constant
            'L/D':                  18.4, # 18.4 w/out HEX | from Elias
            
        # mission parameters
            'f_{res}':              1.199921,
            #'R_{des}':              2935 * units('nmi'),
            #'m_{pay}':              35000 * units('lb'),
            #'P_{\\rm out}':       700*units('kW'),
            #'m_{\\rm OE}':          100000*units('lb'),
    })
    
    init = {
        mission_CHEETA_he['(i0/A)']:                   1e-3*units('A/cm**2'),
        mission_CHEETA_he['P_{comp}']:                 1e12*units('kW'),
        mission_CHEETA_he['V_{act}']:                  .3*units('V'),
        mission_CHEETA_he['T_{t\\2}']:                 (180+273.15)*units('K'),
        mission_CHEETA_he['\\dot{m}_{a}']:             360*units('kg/s'),
        }
    
    model = Model(mission_CHEETA_he['m_{fuel}'], [ac_CHEETA_he, mission_CHEETA_he])  
    #model = Model(mission['m_{\\rm fuel}'], Bounded([ac, mission]))
    
    sol_CHEETA_he = model.localsolve(x0=init, solver='mosek_conif', verbosity=2, iteration_limit=1000)
    #print(sol_CHEETA_he.table())
    
    CHEETA_he_PFEI = sol_CHEETA_he['variables']['PFEI']
    
plt.figure(1)
line1, = plt.plot(P_m, PFEI_737*np.ones(len(P_m)), 'k:', label='Kerosene Turbo-fan')
line2, = plt.plot(P_m, PFEI_737_te*np.ones(len(P_m)), 'k--', label='Kerosene Turbo-electric')
line3, = plt.plot(P_m, ZeroE_PFEI*np.ones(len(P_m)), 'b:', label='Hydrogen Turbo-fan')
line4, = plt.plot(P_m, CHEETA_te_PFEI*np.ones(len(P_m)), 'b--', label='Hydrogen Turbo-electric')
line5, = plt.plot(P_m, CHEETA_PFEI, 'b', label='Hydrogen Fully-electric')
line6, = plt.plot(P_m, CHEETA_he_PFEI, 'b-.', label='Hydrogen Hybrid-electric')
plt.xlabel("PEMFC Specific Power (W/kg)")
plt.ylabel("PFEI (kJ/kg-km)")
plt.ylim([0,8])
plt.xlim([2000, 3000])
plt.legend(handles=[line1, line2, line3, line4, line5, line6])