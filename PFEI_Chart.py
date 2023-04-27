import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import re

from gpkit import Model, units
from SensSort import SensSort

from path import h_alt, hdot_in, M, P, T, rho, mu, Cp, LDR, C_D_fuse, C_D_wing
from flight_state import FlightState

from aircraft import Aircraft
from mission import Mission

from turbofan import Turbofan
from turboelectric import Turboelectric
from hybrid_electric import HybridElectric
from fully_electric import FullyElectric
from fully_electric_bat import FullyElectricBatHybrid

# GLOBAL CONSTANTS
gamma = 1.4
g = 9.81*units('m/s/s')
R2 = 8.314*units('J/K')
F = 96485.33*units('C')

# Is BLI on or off?
BLI_737 = 'on'
BLI_737_te = 'on'
BLI_ZeroE = 'on'
BLI_CHEETA = 'on'
BLI_CHEETA_te = 'on'
BLI_CHEETA_he = 'on'
BLI_CHEETA_bat = 'on'

# Initialization
performance = np.zeros([3,7])
MTOW = np.zeros([3,7])

# PARAMETER SWEEPS

# Number of engines
N_eng = 2

# Number of motor-driven propulsors
N_p = 9

# design range
R_des = 5500 # km

# design payload weight
m_pay = 15875 # kg

# PEMFC stack specific power
P_m = np.asarray([2000, 2425, 3000]) # W/kg

# GG thermal efficiency
eta_th = np.asarray([0.5, 0.55, 0.6])

b_fuse = 14 # ft
b_wing = 136 - b_fuse # ft
S_ref = 1840 # sq ft
d_fan = 47.24 # inches


for oo in range(len(P_m)):

    #   PROPULSION SYSTEMS
    tf = Turbofan(N_eng)
    tf.substitutions.update({
        # gas generator constants
        #'(P/m)_{gg}': 7687 * units('hp') / (6130/2) / units('lb'), # assuming half of LEAP-1B weight is gas turbine and half is fan
        '(P/m)_{gg}': 6380 * units('hp') / (5216/2) / units('lb'), # assuming half of CFM-56 weight is gas turbine and half is fan
        # ducted fan constants
        #'K_{fan}': 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B weight is turbine and half is fan
        'K_{fan}': (5216/2) * units('lb') / (61 * units('in'))**3, # half of CFM-56 weight is turbine and half is fan
        'HTR': 0.3,
        #'d_{fan}': 69.4*units('in'), # LEAP-1B Fan Diameter
        'd_{fan}': 61*units('in'), # CFM-56 Fan Diameter
        })
    
    te_non = Turboelectric(N_eng, N_p) # non-superconducting turboelectric
    te_non.substitutions.update({
        # gas generator constants
        #'(P/m)_{gg}': 7687 * units('hp') / 3065 / units('lb'), # half of LEAP-1B
        '(P/m)_{gg}': 6380 * units('hp') / (5216/2) / units('lb'), # assuming half of CFM-56 weight is gas turbine and half is fan
        # generator constants
        '(T/m)_{gen}': 20 * units('N*m/kg'),
        # motor constants
        '(T/m)_{mot}': 20 * units('N*m/kg'),  # 70 for superconducting; 20 for traditional
        # ducted fan constants
        #'K_{fan}': 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
        'K_{fan}': (5216/2) * units('lb') / (61 * units('in'))**3, # half of CFM-56 weight is turbine and half is fan
        'HTR': 0.38, 
        'd_{fan}': d_fan*units('in'),
        })
    
    
    fe = FullyElectric(N_p)
    fe.substitutions.update({
        # ducted fan constants
        #'K_{fan}':                 0.5 * 6130 * units('lb') / ((69.4 * units('in'))**3), # half of LEAP-1B
        'K_{fan}': (5216/2) * units('lb') / (61 * units('in'))**3, # half of CFM-56 weight is turbine and half is fan
        'HTR':                     0.38,
        'd_{fan}':                 d_fan*units('in'), # current CHEETA result
        
         # motor constants
        '(T/m)_{mot}':             70*units('N*m/kg'),  # anticipated value for superconducting machines [] | CHEETA table 47.7465
        
        # fuel cell constants
        #'(P/m)_{fc}':               2700*units('W/kg'), # Hypoint
        'N_{cells}':                1080, # number of cells in a stack
        'A':                        446*units('cm**2'), # membrane area
        'P_{max,fc}':               850*units('kW'), 
        #'vol_{fc}':                 37.5*units('cm**3'), # volume of one cell [Wilson, 1996]
        'm_{comp}':                 .00001*units('kg'),
        #'N_s':                      18,
        
        # HEX constants
        '\psi':                    5.02*units('kg/ft^2'), # from Boeing (Chellappa)
       # 'A_{face}':                700*units('in**2'),
    })
        
    te_sc = Turboelectric(N_eng, N_p)
    te_sc.substitutions.update({
        # gas generator constants
        #'(P/m)_{gg}': 7687 * units('hp') / 3065 / units('lb'), # LEAP-1B
        '(P/m)_{gg}': 6380 * units('hp') / (5216/2) / units('lb'), # assuming half of CFM-56 weight is gas turbine and half is fan
        # generator constants
        '(T/m)_{gen}': 20 * units('N*m/kg'),
        # motor constants
        '(T/m)_{mot}': 70 * units('N*m/kg'),  # 70 for superconducting; 20 for traditional
        # ducted fan constants
        #'K_{fan}': 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
        'K_{fan}': (5216/2) * units('lb') / (61 * units('in'))**3, # half of CFM-56 weight is turbine and half is fan
        'HTR': 0.38, 
        'd_{fan}': d_fan*units('in'),
        })
    
    
    he = HybridElectric(N_p)
    he.substitutions.update({
        
        # ducted fan constants
        #'K_{fan}':                 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
        'K_{fan}': (5216/2) * units('lb') / (61 * units('in'))**3, # half of CFM-56 weight is turbine and half is fan
        'HTR':                     0.38,
        'd_{fan}':                  d_fan*units('in'),
         # gas generator constants
        #'(P/m)_{gg}':               7687 * units('hp') / 3065 / units('lb'), # LEAP-1B
        '(P/m)_{gg}': 6380 * units('hp') / (5216/2) / units('lb'), # assuming half of CFM-56 weight is gas turbine and half is fan
        # generator constants
        '(T/m)_{gen}':              20 * units('N*m/kg'),
         # motor constants
        '(T/m)_{mot}':             70*units('N*m/kg'),  # anticipated value for superconducting machines []
        # fuel cell constants
        'P_{max,fc}':               850*units('kW'), # max power the model is accurate for. This could be replaced with a max power specified by Boeing
        'm_{comp}':                 .00001*units('kg'),
        'N_s':                      6,
        'N_{eng}':                  2,
        #'//nu_{conv}':              .8,
        # HEX constants
        '\psi':                    5.02*units('kg/ft^2'), # from Boeing (Chellappa)
        })
    
    feb = FullyElectricBatHybrid(N_p)
    feb.substitutions.update({
        # motor constants
        '(T/m)_{mot}': 70 * units('N*m/kg'),  # 70 for superconducting; 20 for traditional
        # ducted fan constants
        #'K_{fan}': 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
        'K_{fan}': (5216/2) * units('lb') / (61 * units('in'))**3, # half of CFM-56 weight is turbine and half is fan
        'HTR': 0.38, 
        'd_{fan}': d_fan*units('in'),
        
        # fuel cell constants
        'P_{max,fc}':               850*units('kW'), # max power the model is accurate for. This could be replaced with a max power specified by Boeing
        'm_{comp}':                 .00001*units('kg'),
        #'N_s':                      6,
        # HEX constants
        '\psi':                    5.02*units('kg/ft^2'), # from Boeing (Chellappa)
        # battery constants
        '(P/m)_{bat}': 800*units('W/kg'), # conservative 2035 future estimate from LEARN report pg. 14
        '(E/m)_{\\rm bat}': 8*10**2*units('W*hr/kg'), # 2050 estimate of Li ion batteries based on 4% annual growth rate [Viswanathan, Nature, 2022]
        'N_b':  4,
        'm_{bat}': 100*units('kg'),
    })
    
    
    ########################
    
    
    #   AIRCRAFT SYSTEMS
    ac_737 = Aircraft(tf)
    ac_737.substitutions.update({
        'f_{empty}': 0.47,  # 737-800 value for 3000 nmi mission
        #'f_{empty}': 0.51,  # 737 MAX 8 value for 3000 nmi mission
        #'d_{fan}': 69.4*units('in'), # for 737 baseline with LEAP-1B twin engines
        
        'b_{fuse}': b_fuse*units('ft'),
        'b_{wing}': b_wing*units('ft'),
        'S_{ref}': S_ref*units('ft^2'),
        })
    
    ac_737_te = Aircraft(te_non)
    ac_737_te.substitutions.update({
        'f_{empty}': 0.47,  # (CHEETA empty weight minus propulsion system weight) divided by MTOW
        
        'b_{fuse}': b_fuse*units('ft'),
        'b_{wing}': b_wing*units('ft'),
        'S_{ref}': S_ref*units('ft^2'),
        })
    
    ac_ZeroE = Aircraft(tf)
    ac_ZeroE.substitutions.update({
            'f_{empty}':            0.554, # (CHEETA empty weight minus propulsion system weight) divided by MTOW
            # tuned to Elias' weights
            'b_{fuse}': b_fuse*units('ft'),
            'b_{wing}': b_wing*units('ft'),
            'S_{ref}': S_ref*units('ft^2'),
        })
    
    ac_CHEETA = Aircraft(fe)
    ac_CHEETA.substitutions.update({
        'f_{empty}':                0.554, # (CHEETA empty weight minus propulsion system weight) divided by MTOW
        
        'b_{fuse}': b_fuse*units('ft'),
        'b_{wing}': b_wing*units('ft'),
        'S_{ref}': S_ref*units('ft^2'),
        })
    
    ac_CHEETA_te = Aircraft(te_sc)
    ac_CHEETA_te.substitutions.update({
        'f_{empty}':                0.554, # (CHEETA empty weight minus propulsion system weight) divided by MTOW
        
        'b_{fuse}': b_fuse*units('ft'),
        'b_{wing}': b_wing*units('ft'),
        'S_{ref}': S_ref*units('ft^2'),
        })
    
    ac_CHEETA_he = Aircraft(he)
    ac_CHEETA_he.substitutions.update({
        'f_{empty}':                0.554, # (CHEETA empty weight minus propulsion system weight) divided by MTOW
        
        'b_{fuse}': b_fuse*units('ft'),
        'b_{wing}': b_wing*units('ft'),
        'S_{ref}': S_ref*units('ft^2'),
        })
    
    ac_CHEETA_bat = Aircraft(feb)
    ac_CHEETA_bat.substitutions.update({
        'f_{empty}':                0.554, # (CHEETA empty weight minus propulsion system weight) divided by MTOW
        
        'b_{fuse}': b_fuse*units('ft'),
        'b_{wing}': b_wing*units('ft'),
        'S_{ref}': S_ref*units('ft^2'),
        })
    
    #%%
    # CONVENTIONAL 737 MISSION

    mission_737 = Mission(ac_737, h_alt, M, P, T, rho, mu, Cp, hdot_in, LDR, C_D_fuse, C_D_wing, BLI_737)
    # all aircraft and propulsion system performance models (and thus their
    # variables) get set up inside Mission(), so we update their constants here
    
    mission_737.substitutions.update({
        '\rho_{fuel}':              840*units('kg/m^3'), # jet fuel density
        'h_{fuel}':                43 * units('MJ/kg'), # jet fuel lower heating value
        
        # mission parameters
        'f_{res}': 1.2, # assume 15% additional fuel for reserves
        'R_{des}': R_des*units('km'), # design range
        'm_{pay}':  m_pay*units('kg'), # design payload
    
        # aircraft performance parameters
        #'D_corr':           0.5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        'L/D':              LDR, # take input lift-to-drag ratios from Elias

        '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]; according to D.K. Hall, this number is generally representative         
        '\phi_{fuse}':               .5, # assume half of fuselage dissipation occurs on top of fuselage 
        '\phi_{wing}':               .64, # assume most dissipation occurs on top of wing
        
        # gg parameters
        '\\eta_{th}':                eta_th[oo],
        #'\\eta_{fan}':               0.95, 
        #'Cd':                        0.04, # [Robinson, 2017], nacelle drag coefficient
        
        })
    
    init = {
        mission_737['F_{net}']:                       44*units('kN'),
        mission_737['u_{jet}']:                       10000*units('m/s'),
        mission_737['dt']:                                           1000*units('s'),
        mission_737.fs.performance.prop.fan_perf['F_{fan}']:       100000*units('kN'),
        mission_737['F_{net}']:                                         44*units('kN'),
        }
    
    model = Model(mission_737['PFEI'], [ac_737, mission_737])
    
    if BLI_737 == 'on':
        sol_737 = model.localsolve(x0=init, solver='mosek_conif', verbosity=2)
    else:
        sol_737 = model.solve(x0=init, solver='mosek_conif', verbosity=2)
        
    #(sol_737.table())
    #print(sol_737.table())    
    
    PFEI_737 = sol_737['variables']['PFEI']
    MTOW_737 = sol_737['variables']['MTOW']
    sens_737 = list(sol_737['sensitivities']['variables'].values())
    
    sens_737 = SensSort(sens_737)
    
    list_variables_737 = list(sol_737['sensitivities']['variables'].keys())
    
    variables_737 = ['']*len(list_variables_737)
    jj = 0
    
    while jj < len(list_variables_737):
        variables_737[jj] = str(list_variables_737[jj])
        variables_737[jj] = re.sub("^.*\\.","", variables_737[jj])
        
        jj = jj + 1
    
    
    #%%
    # TURBOELECTRIC 737 MISSION

    mission_737_te = Mission(ac_737_te, h_alt, M, P, T, rho, mu, Cp, hdot_in, LDR, C_D_fuse, C_D_wing, BLI_737_te)
    # all aircraft and propulsion system performance models (and thus their
    # variables) get set up inside Mission(), so we update their constants here
    mission_737_te.substitutions.update({
        '\rho_{fuel}':              840*units('kg/m^3'), # jet fuel density
        'h_{fuel}':                43 * units('MJ/kg'), # jet fuel lower heating value
        
        # mission parameters
        'f_{res}': 1.2, # assume 15% additional fuel for reserves
        'R_{des}': R_des*units('km'), # design range
        'm_{pay}':  m_pay*units('kg'), # design payload
    
        # aircraft performance parameters
        'L/D':              LDR, # take input lift-to-drag ratios from Elias
        #'f_{BLI,max}':                   0.3388, 
        '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]; according to D.K. Hall, this number is generally representative         
        #'D_p':                       4830*units('lbf'), # profile drag ingested by propulsors
        '\phi_{fuse}':               .5, # assume half of fuselage dissipation occurs on top of fuselage 
        '\phi_{wing}':               .64, # assume most dissipation occurs on top of wing
    
        # gg parameters
        '\\eta_{th}':                eta_th[oo],
        
        # fan parameters
        #'D_corr':           0.5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'\\eta_{fan}':               0.95, 
        #'Cd':                        0.04, # [Robinson, 2017]
        
        # generator constants
        'omega_{gen}': 16000 * units('rpm'),
        '1-eta_{gen}': 1 - .98,
        
        # motor constants
        'omega_{motor}': 4000 * units('rpm'),
        '1-eta_{motor}': 1 - .95, #.999 for superconducting; .95 for traditional
        
        })
    
    init = {
        mission_737_te['F_{net}']:                       44*units('kN'),
        mission_737_te.fs.performance.prop.fan_fuse_perf['u_{jet}']:       100000*units('m/s'),
        mission_737_te.fs.performance.prop.fan_wing_perf['u_{jet}']:       100000*units('m/s'),
        mission_737_te['dt']:                                           1000*units('s'),
        mission_737_te.fs.performance.prop.fan_fuse_perf['F_{fan}']:       100000*units('kN'),
        mission_737_te.fs.performance.prop.fan_wing_perf['F_{fan}']:       100000*units('kN'),
        }
    
    model = Model(mission_737_te['PFEI'], [ac_737_te, mission_737_te])
    
    sol_737_te = model.localsolve(x0=init, solver='mosek_conif', verbosity=2)
    #print(sol_737_te.table())
    
    PFEI_737_te = sol_737_te['variables']['PFEI']
    MTOW_737_te = sol_737_te['variables']['MTOW']
    sens_737_te = list(sol_737_te['sensitivities']['variables'].values())
    sens_737_te = SensSort(sens_737_te)
    list_variables_737_te = list(sol_737_te['sensitivities']['variables'].keys())
    
    variables_737_te = ['']*len(list_variables_737_te)
    jj = 0
    
    while jj < len(list_variables_737_te):
        variables_737_te[jj] = str(list_variables_737_te[jj])
        variables_737_te[jj] = re.sub("^.*\\.","", variables_737_te[jj])
        
        jj = jj + 1
    
    
    #%%
    # H2-BURNING TURBOFAN (A LA AIRBUS ZERO-E)
    mission_ZeroE = Mission(ac_ZeroE, h_alt, M, P, T, rho, mu, Cp, hdot_in, LDR, C_D_fuse, C_D_wing, BLI_ZeroE)
    # all aircraft and propulsion system performance models (and thus their
    # variables) get set up inside Mission(), so we update their constants here
    
    mission_ZeroE.substitutions.update({
        '\rho_{fuel}':              71*units('kg/m^3'), # jet fuel density
        'h_{fuel}':                120 * units('MJ/kg'), # jet fuel lower heating value
        
        # mission parameters
        'f_{res}': 1.2, # assume 15% additional fuel for reserves
        'R_{des}': R_des*units('km'), # design range
        'm_{pay}':  m_pay*units('kg'), # design payload
    
    # aircraft performance parameters
        #'D_corr':           0.5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        'L/D':              LDR, # take input lift-to-drag ratios from Elias
        #'f_{BLI,max}':                   0.0001, # assuming no BLI
        #'f_{BLI,max}':                   0.17, # D8 BLI
        '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]; according to D.K. Hall, this number is generally representative         
        #'D_p':                       4830*units('lbf'), # profile drag ingested by propulsors
        '\\eta_{th}':               eta_th[oo],
        #'\\eta_{fan}':               0.95, 
        #'Cd':                        0.04, # [Robinson, 2017]
        '\phi_{fuse}':               .5, # assume half of fuselage dissipation occurs on top of fuselage 
        '\phi_{wing}':               .64, # assume most dissipation occurs on top of wing
        
        })
    
    init = {
        mission_ZeroE['F_{net}']:                       44*units('kN'),
        mission_ZeroE['u_{jet}']:                       100*units('m/s'),
        mission_ZeroE['dt']:                                           1000*units('s'),
        mission_ZeroE.fs.performance.prop.fan_perf['F_{fan}']:       100000*units('kN'),
        mission_ZeroE.fs.performance.prop.fan_perf['\\dot{m}']:       10000*units('kg/s'),
        }
    
    model = Model(mission_ZeroE['PFEI'], [ac_ZeroE, mission_ZeroE])
    
    if BLI_ZeroE == 'on':
        sol_ZeroE = model.localsolve(x0=init, solver='mosek_conif', verbosity=2)
    else:
        sol_ZeroE = model.solve(x0=init, solver='mosek_conif', verbosity=2)
    #print(sol_ZeroE.table())
    
    ZeroE_PFEI = sol_ZeroE['variables']['PFEI']
    ZeroE_MTOW = sol_ZeroE['variables']['MTOW']
    sens_ZeroE = list(sol_ZeroE['sensitivities']['variables'].values())
    sens_ZeroE = SensSort(sens_ZeroE)
    list_variables_ZeroE = list(sol_ZeroE['sensitivities']['variables'].keys())
    
    variables_ZeroE = ['']*len(list_variables_ZeroE)
    jj = 0
    
    while jj < len(list_variables_ZeroE):
        variables_ZeroE[jj] = str(list_variables_ZeroE[jj])
        variables_ZeroE[jj] = re.sub("^.*\\.","", variables_ZeroE[jj])
        
        jj = jj + 1
    
    
   #%% 
    # FULLY-ELECTRIC H2 PEMFC (A LA CHEETA)
    mission_CHEETA = Mission(ac_CHEETA, h_alt, M, P, T, rho, mu, Cp, hdot_in, LDR, C_D_fuse, C_D_wing, BLI_CHEETA)    
    mission_CHEETA.substitutions.update({
        '\rho_{fuel}':              71*units('kg/m^3'), # jet fuel density
        'h_{fuel}':                120 * units('MJ/kg'), # jet fuel lower heating value
        
        # mission parameters
        'f_{res}': 1.2, # assume 15% additional fuel for reserves
        'R_{des}': R_des*units('km'), # design range
        'm_{pay}':  m_pay*units('kg'), # design payload
    
    # aircraft performance parameters
        'L/D':              LDR, # take input lift-to-drag ratios from Elias
        #'f_{BLI,max}':                   0.3388, 
        '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]; according to D.K. Hall, this number is generally representative         
        #'D_p':                       4830*units('lbf'), # profile drag ingested by propulsors
        '\phi_{fuse}':               .5, # assume half of fuselage dissipation occurs on top of fuselage 
        '\phi_{wing}':               .64, # assume half of dissipation occurs on top of wing
        
        # fan constants
        #'D_corr':           0.5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'\\eta_{fan}':               0.95, 
        #'Cd':                        0.04, # [Robinson, 2017]
        
        # fuel cell constants
        '(P/m)_{fc}':               P_m[oo]*units('W/kg'), # Hypoint
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
       # 'K_f':                        2,
       # 'K_h':                        0.93,
       # '\sigma':                     .464,
       # '\epsilon':                   0.85, # Kays and London
       'C_{p,h}': 2200*units('J/kg/K'),
       '\dot{m}_h': 15*units('kg/s'),
       'h': 75*units('W/m^2/K'), # based on lowest value given by Chellappa (.5 MW heat load during cruise conditions)
        
        # motor constants
        'omega_{motor}': 4000 * units('rpm'),
        '1-eta_{motor}': 1 - .999, #.999 for superconducting; .95 for traditional
        })
    
    init = {
        mission_CHEETA['(i0/A)']:                                          1e-3*units('A/cm**2'),
        mission_CHEETA['P_{comp}']:                                        1e12*units('kW'),
        mission_CHEETA['V_{act}']:                                         .3*units('V'),
        mission_CHEETA['T_{t\\2}']:                                        (180+273.15)*units('K'),
        mission_CHEETA['\\dot{m}_{a}']:                                    360*units('kg/s'),
        mission_CHEETA['eta_{comp}']:                                      .9,
        mission_CHEETA['F_{net}']:                                         44*units('kN'),
        mission_CHEETA.fs.performance.prop.fan_fuse_perf['u_{jet}']:       100*units('m/s'),
        mission_CHEETA.fs.performance.prop.fan_wing_perf['u_{jet}']:       100*units('m/s'),
        mission_CHEETA.fs.performance.prop.fan_fuse_perf['\\dot{m}']:              1000000*units('kg/s'),
        mission_CHEETA.fs.performance.prop.fan_wing_perf['\\dot{m}']:              1000000*units('kg/s'),
        mission_CHEETA['pi_c']:                                            5,
        mission_CHEETA['F_{HEX}']:                                         1*units('kN'),
        mission_CHEETA['P_{gross}']:                                       1*units('kW'),
        mission_CHEETA['dt']:                                           1000*units('s'),
        mission_CHEETA.fs.performance.prop.fan_fuse_perf['F_{fan}']:       100000*units('kN'),
        mission_CHEETA.fs.performance.prop.fan_wing_perf['F_{fan}']:       100000*units('kN'),
        }
    
    model = Model(mission_CHEETA['PFEI'], [ac_CHEETA, mission_CHEETA])  
    
    sol_CHEETA = model.localsolve(x0=init, solver='mosek_conif', verbosity=2, iteration_limit=1000)
    #print(sol_CHEETA.table())
    
    CHEETA_PFEI = sol_CHEETA['variables']['PFEI']
    CHEETA_MTOW = sol_CHEETA['variables']['MTOW']
    CHEETA_sens = list(sol_CHEETA['sensitivities']['variables'].values())
    CHEETA_sens = SensSort(CHEETA_sens)
    list_CHEETA_variables = list(sol_CHEETA['sensitivities']['variables'].keys())
    
    CHEETA_variables = ['']*len(list_CHEETA_variables)
    jj = 0
    
    while jj < len(list_CHEETA_variables):
        CHEETA_variables[jj] = str(list_CHEETA_variables[jj])
        CHEETA_variables[jj] = re.sub("^.*\\.","", CHEETA_variables[jj])
        
        jj = jj + 1
    
        
    
    
    # CHEETA TURBO-ELECTRIC
    mission_CHEETA_te = Mission(ac_CHEETA_te, h_alt, M, P, T, rho, mu, Cp, hdot_in, LDR, C_D_fuse, C_D_wing, BLI_CHEETA_te)
    # all aircraft and propulsion system performance models (and thus their
    # variables) get set up inside Mission(), so we update their constants here
    
    mission_CHEETA_te.substitutions.update({
        '\rho_{fuel}':              71*units('kg/m^3'), # jet fuel density
        'h_{fuel}':                120 * units('MJ/kg'), # jet fuel lower heating value
        
        # mission parameters
        'f_{res}': 1.2, # assume 15% additional fuel for reserves
        'R_{des}': R_des*units('km'), # design range
        'm_{pay}':  m_pay*units('kg'), # design payload
    
    # aircraft performance parameters
        'L/D':              LDR, # take input lift-to-drag ratios from Elias
        #'f_{BLI,max}':                   0.3388, 
        '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]; according to D.K. Hall, this number is generally representative         
        #'D_p':                       4830*units('lbf'), # profile drag ingested by propulsors
        '\phi_{fuse}':               .5, # assume half of fuselage dissipation occurs on top of fuselage 
        '\phi_{wing}':               .64, # assume most dissipation occurs on top of wing
        
        '\\eta_{th}':                eta_th[oo],
        
       # 'D_corr':           0.5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
       # '\\eta_{fan}':               0.95, 
       # 'Cd':                        0.04, # [Robinson, 2017]
        
        # generator constants
        'omega_{gen}': 16000 * units('rpm'),
        '1-eta_{gen}': 1 - .98,
        
        # motor constants
        'omega_{motor}': 4000 * units('rpm'),
        '1-eta_{motor}': 1 - .999, #.999 for superconducting; .95 for traditional
        
        })
    
    
    init = {
        mission_CHEETA_te['F_{net}']:                       44*units('kN'),
        mission_CHEETA_te.fs.performance.prop.fan_fuse_perf['u_{jet}']:       100000*units('m/s'),
        mission_CHEETA_te.fs.performance.prop.fan_wing_perf['u_{jet}']:       100000*units('m/s'),
        mission_CHEETA_te['dt']:                                           1000*units('s'),
        mission_CHEETA_te.fs.performance.prop.fan_fuse_perf['F_{fan}']:       100000*units('kN'),
        mission_CHEETA_te.fs.performance.prop.fan_wing_perf['F_{fan}']:       100000*units('kN'),
        }
    
    model = Model(mission_CHEETA_te['PFEI'], [ac_CHEETA_te, mission_CHEETA_te])
    
    sol_CHEETA_te = model.localsolve(x0=init,solver='mosek_conif',verbosity=2)
    #print(sol_CHEETA_te.table())
    
    CHEETA_te_PFEI = sol_CHEETA_te['variables']['PFEI']
    CHEETA_te_MTOW = sol_CHEETA_te['variables']['MTOW']
    CHEETA_te_sens = list(sol_CHEETA_te['sensitivities']['variables'].values())
    CHEETA_te_sens = SensSort(CHEETA_te_sens)
    list_CHEETA_te_variables = list(sol_CHEETA_te['sensitivities']['variables'].keys())
    
    CHEETA_te_variables = ['']*len(list_CHEETA_te_variables)
    jj = 0
    
    while jj < len(list_CHEETA_te_variables):
        CHEETA_te_variables[jj] = str(list_CHEETA_te_variables[jj])
        CHEETA_te_variables[jj] = re.sub("^.*\\.","", CHEETA_te_variables[jj])
        
        jj = jj + 1
    
        
    #%%
    
    
    # CHEETA HYBRID-ELECTRIC
    mission_CHEETA_he = Mission(ac_CHEETA_he, h_alt, M, P, T, rho, mu, Cp, hdot_in, LDR, C_D_fuse, C_D_wing, BLI_CHEETA_he)
        # all aircraft and propulsion system performance models (and thus their
        # variables) get set up inside Mission(), so we update their constants here
        
    mission_CHEETA_he.substitutions.update({
        '\rho_{fuel}':              71*units('kg/m^3'), # jet fuel density
        'h_{fuel}':                120 * units('MJ/kg'), # jet fuel lower heating value
        
        # mission parameters
        'f_{res}': 1.2, # assume 15% additional fuel for reserves
        'R_{des}': R_des*units('km'), # design range
        'm_{pay}':  m_pay*units('kg'), # design payload
    
    # aircraft performance parameters
        'L/D':              LDR, # take input lift-to-drag ratios from Elias
       # 'f_{BLI,max}':                   0.3388, 
        '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]; according to D.K. Hall, this number is generally representative         
       # 'D_p':                       4830*units('lbf'), # profile drag ingested by propulsors
        
        '\\eta_{th}':                eta_th[oo],
        
        # fan constants
       # 'D_corr':           0.5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
       # '\\eta_{fan}':               0.95, 
       ## 'Cd':                        0.04, # [Robinson, 2017]
        '\phi_{fuse}':               .5, # assume half of fuselage dissipation occurs on top of fuselage 
        '\phi_{wing}':               .64, # assume most dissipation occurs on top of wing
        
        # fuel cell constants
        '(P/m)_{fc}':               P_m[oo]*units('W/kg'), # Hypoint
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
      #  'K_f':                        2,
      #  'K_h':                        0.93,
       # '\sigma':                     .464,
      #  '\epsilon':                   0.85, # Kays and London
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
        mission_CHEETA_he['(i0/A)']:                                              1e-3*units('A/cm**2'),
        mission_CHEETA_he['P_{comp}']:                                            1e12*units('kW'),
        mission_CHEETA_he['V_{act}']:                                             .3*units('V'),
        mission_CHEETA_he['T_{t\\2}']:                                            (180+273.15)*units('K'),
        mission_CHEETA_he['\\dot{m}_{a}']:                                        360*units('kg/s'),
        mission_CHEETA_he['eta_{comp}']:                                          .9,
        mission_CHEETA_he['F_{net}']:                                             44*units('kN'),
        mission_CHEETA_he.fs.performance.prop.fan_fuse_perf['u_{jet}']:           10*units('m/s'),
        mission_CHEETA_he.fs.performance.prop.fan_wing_perf['u_{jet}']:           10*units('m/s'),
        mission_CHEETA_he.fs.performance.prop.fan_fuse_perf['\\dot{m}']:          1000000*units('kg/s'),
        mission_CHEETA_he.fs.performance.prop.fan_wing_perf['\\dot{m}']:          1000000*units('kg/s'),
        mission_CHEETA_he['f_{fuse,BLI}']:                                        .1,
        mission_CHEETA_he['f_{wing,BLI}']:                                        .01,
        mission_CHEETA_he['pi_c']:                                                5,
        mission_CHEETA_he['F_{HEX}']:                                             1*units('kN'),
        mission_CHEETA_he['D_{core}']:                                            1*units('kN'),
        mission_CHEETA_he['P_{gross}']:                                           200*units('kW'),
        mission_CHEETA_he['dt']:                                               1000*units('s'),
        mission_CHEETA_he.fs.performance.prop.fan_fuse_perf['F_{fan}']:           10*units('kN'),
        mission_CHEETA_he.fs.performance.prop.fan_wing_perf['F_{fan}']:           10*units('kN'),
        mission_CHEETA_he.fs.performance.prop.gen_perf['P_{out}']:                200*units('kW'),
        mission_CHEETA_he.fs.performance.prop.motor_perf['P_{out}']:              200*units('kW'),
        mission_CHEETA_he.fs.performance.prop.stack_perf['P_{out}']:              200*units('kW'),
        mission_CHEETA_he.fs.performance.prop.core_perf['P_{out}']:               200*units('kW'),
        mission_CHEETA_he['d_{fan}']:                                             56*units('in'),
        mission_CHEETA_he['N_s']:                                               1e-12,
       # mission_CHEETA_he['N_s']:                                                 600,
        }
    
    
    model = Model(mission_CHEETA_he['PFEI'], [ac_CHEETA_he, mission_CHEETA_he])  
    #model = Model(mission['m_{\\rm fuel}'], Bounded([ac, mission]))
    
    sol_CHEETA_he = model.localsolve(x0=init, solver='mosek_conif', verbosity=2, iteration_limit=1000)
    #print(sol_CHEETA_he.table())
    
    CHEETA_he_PFEI = sol_CHEETA_he['variables']['PFEI']
    CHEETA_he_MTOW = sol_CHEETA_he['variables']['MTOW']
    CHEETA_he_sens = list(sol_CHEETA_he['sensitivities']['variables'].values())
    CHEETA_he_sens = SensSort(CHEETA_he_sens)
    list_CHEETA_he_variables = list(sol_CHEETA_he['sensitivities']['variables'].keys())
    
    CHEETA_he_variables = ['']*len(list_CHEETA_he_variables)
    jj = 0
    
    while jj < len(list_CHEETA_he_variables):
        CHEETA_he_variables[jj] = str(list_CHEETA_he_variables[jj])
        CHEETA_he_variables[jj] = re.sub("^.*\\.","", CHEETA_he_variables[jj])
        
        jj = jj + 1
        
        #%%
    # FULLY-ELECTRIC H2 PEMFC with Battery 
    mission_CHEETA_bat = Mission(ac_CHEETA_bat, h_alt, M, P, T, rho, mu, Cp, hdot_in, LDR, C_D_fuse, C_D_wing, BLI_CHEETA_bat)    
    mission_CHEETA_bat.substitutions.update({
        '\rho_{fuel}':              71*units('kg/m^3'), # jet fuel density
        'h_{fuel}':                120 * units('MJ/kg'), # jet fuel lower heating value
        
        # mission parameters
        'f_{res}': 1.2, # assume 15% additional fuel for reserves
        'R_{des}': R_des*units('km'), # design range
        'm_{pay}':  m_pay*units('kg'), # design payload
    
    # aircraft performance parameters
        'L/D':              LDR, # take input lift-to-drag ratios from Elias
        #'f_{BLI,max}':                   0.3388, 
        '1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]; according to D.K. Hall, this number is generally representative         
        #'D_p':                       4830*units('lbf'), # profile drag ingested by propulsors
        '\phi_{fuse}':               .5, # assume half of fuselage dissipation occurs on top of fuselage 
        '\phi_{wing}':               .64, # assume most dissipation occurs on top of wing
        
        # fan constants
        #'D_corr':           0.5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        #'\\eta_{fan}':               0.95, 
        #'Cd':                        0.04, # [Robinson, 2017]
        
        # fuel cell constants
        '(P/m)_{fc}':               P_m[oo]*units('W/kg'), # Hypoint
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
       # '\sigma':                     .464,
        #'\epsilon':                   0.85, # Kays and London
        'C_{p,h}': 2200*units('J/kg/K'),
         '\dot{m}_h': 15*units('kg/s'),
         'h': 75*units('W/m^2/K'), # based on lowest value given by Chellappa (.5 MW heat load during cruise conditions)
        
        # motor constants
        'omega_{motor}': 4000 * units('rpm'),
        '1-eta_{motor}': 1 - .999, #.999 for superconducting; .95 for traditional
        })
    
    mission_CHEETA_bat.append(mission_CHEETA_bat['E_{bat,max}'] >= mission_CHEETA_bat['E_{bat}'][0] + mission_CHEETA_bat['E_{bat}'][1] + mission_CHEETA_bat['E_{bat}'][2] + mission_CHEETA_bat['E_{bat}'][3] + mission_CHEETA_bat['E_{bat}'][4] + mission_CHEETA_bat['E_{bat}'][5] + mission_CHEETA_bat['E_{bat}'][6] + mission_CHEETA_bat['E_{bat}'][7] + mission_CHEETA_bat['E_{bat}'][8] + mission_CHEETA_bat['E_{bat}'][9] + mission_CHEETA_bat['E_{bat}'][len(h_alt)-1],)
    
    # these constraints turn off the battery for all operating conditions except top-of-climb
   # for yy in range(len(h_alt)):
    #    if yy == len(h_alt)-2:
     #       mission_CHEETA_bat.append(mission_CHEETA_bat['E_{bat,max}'] >= mission_CHEETA_bat['E_{bat}'][yy])
      #  else:
       #     mission_CHEETA_bat.append(.0001*units('kJ') >= mission_CHEETA_bat['E_{bat}'][yy])

    init = {
        mission_CHEETA_bat['(i0/A)']:                                          1e-3*units('A/cm**2'),
        mission_CHEETA_bat['P_{comp}']:                                        1e12*units('kW'),
        mission_CHEETA_bat['V_{act}']:                                         .3*units('V'),
        mission_CHEETA_bat['T_{t\\2}']:                                        (180+273.15)*units('K'),
        mission_CHEETA_bat['\\dot{m}_{a}']:                                    360*units('kg/s'),
        mission_CHEETA_bat['eta_{comp}']:                                      .9,
        mission_CHEETA_bat['F_{net}']:                                         44*units('kN'),
        mission_CHEETA_bat.fs.performance.prop.fan_fuse_perf['u_{jet}']:       100000*units('m/s'),
        mission_CHEETA_bat.fs.performance.prop.fan_wing_perf['u_{jet}']:       100000*units('m/s'),
        mission_CHEETA_bat['pi_c']:                                            5,
        mission_CHEETA_bat['F_{HEX}']:                                         1*units('kN'),
        mission_CHEETA_bat['P_{gross}']:                                       1*units('kW'),
        mission_CHEETA_bat['dt']:                                           1000*units('s'),
        mission_CHEETA_bat.fs.performance.prop.fan_fuse_perf['F_{fan}']:       100000*units('kN'),
        mission_CHEETA_bat.fs.performance.prop.fan_wing_perf['F_{fan}']:       100000*units('kN'),
        mission_CHEETA_bat.fs.performance.prop.bat_perf['P_{out}']:            1e-12*units('kW'),
        mission_CHEETA_bat['m_{bat}']:                                         100*units('kg'),
        }
    
    model = Model(mission_CHEETA_bat['PFEI'], [ac_CHEETA_bat, mission_CHEETA_bat])  
    
    sol_CHEETA_bat = model.localsolve(x0=init, solver='mosek_conif', verbosity=2)
    #print(sol_CHEETA.table())
    
    CHEETA_bat_PFEI = sol_CHEETA_bat['variables']['PFEI']
    CHEETA_bat_MTOW = sol_CHEETA_bat['variables']['MTOW']
    CHEETA_bat_sens = list(sol_CHEETA_bat['sensitivities']['variables'].values())
    CHEETA_bat_sens = SensSort(CHEETA_bat_sens)
    list_CHEETA_bat_variables = list(sol_CHEETA_bat['sensitivities']['variables'].keys())
    
    CHEETA_bat_variables = ['']*len(list_CHEETA_bat_variables)
    jj = 0
    
    while jj < len(list_CHEETA_bat_variables):
        CHEETA_bat_variables[jj] = str(list_CHEETA_bat_variables[jj])
        CHEETA_bat_variables[jj] = re.sub("^.*\\.","", CHEETA_bat_variables[jj])
        
        jj = jj + 1
        
    font = {'family' : 'normal',
                'weight' : 'normal',
                'size'   : 12}
    
    matplotlib.rc('font', **font)
        
    systems = ('ATF Turbo-fan', 'ATF Turbo-electric', 'Hydrogen Turbo-fan', 'Hydrogen Turbo-electric', 'Hydrogen Hybrid-electric','Hydrogen Fully-electric', 'Hydrogen Fully-electric + Battery')
    
    performance[oo] = [PFEI_737, PFEI_737_te, ZeroE_PFEI, CHEETA_te_PFEI, CHEETA_he_PFEI, CHEETA_PFEI, CHEETA_bat_PFEI]
    MTOW[oo] = np.asarray([MTOW_737, MTOW_737_te, ZeroE_MTOW, CHEETA_te_MTOW, CHEETA_he_MTOW, CHEETA_MTOW, CHEETA_bat_MTOW])/1000 # converted from kg to tonnes
    
    if oo == 1:
        variables = [variables_737, variables_737_te, variables_ZeroE, CHEETA_te_variables, CHEETA_he_variables, CHEETA_variables, CHEETA_bat_variables]
        sensitivities = [sens_737, sens_737_te, sens_ZeroE, CHEETA_te_sens, CHEETA_he_sens, CHEETA_sens, CHEETA_bat_sens]
    
    #%%
# PFEI Bar Chart
plt.figure(1)
y_pos = np.arange(len(systems))

X1 = performance[2]
X2 = np.asarray(performance[1]) - X1
X3 = .062*np.ones(len(systems))
X4 = performance[0] - np.asarray(performance[1]) - X3

X12 = np.add(X1, X2).tolist()
X13 = np.add(X12, X3).tolist()

plt.barh(y_pos, X1, capsize=3, alpha=.8, color='w')
plt.barh(y_pos, X2, capsize=3, alpha=0.5, label='Optimistic', color='b', left=X1)
plt.barh(y_pos, X3, capsize=3, label='Projected', color='k', left=X12)
plt.barh(y_pos, X4, capsize=3, alpha=0.65, label='Conservative', color='#808080', left=X13)

plt.yticks(y_pos, systems)
plt.xlabel("PFEI (kJ/kg-km)")
plt.xlim([0,12])
plt.legend(loc='upper left')

# MTOW Chart
plt.figure(2)
y_pos = np.arange(len(systems))

plt.barh(y_pos, MTOW[1], capsize=3, alpha=.5, color='k')

plt.yticks(y_pos, systems)
plt.xlabel("MTOW (tonnes)")
plt.xlim([0,100])

#%%
# Bar chart of fuel weight sensitivity to a 1% change of each parameter

# Parameter Dictionary
Dict = {'f_{empty}':        'Aircraft Empty Weight Fraction', 
        'f_{res}':          'Fuel Reserve Fraction',
        'R_{des}':          'Flight Range',
        'P_{max,fc}':       'Fuel Cell Stack Maximum Power',
        'm_{pay}':          'Payload Mass',
        'V_\infty[:]':         'Flight Velocity',
        'T_{t\infty}[:]':      'Cruise Ambient Stagnation Temperature',
        'F_{net,max}':      'Take-off Thrust Required',
        'A':                'Fuel Cell Membrane Area',
        'alpha[:]':            'Fuel Cell Charge Transfer Coefficient',
        'N_{cells}':        'Number of Fuel Cells per Stack',
        'L/D[:]':              'Aircraft Lift-to-Drag Ratio',
        '\eta_{fan}[:]':       'Fan Efficiency',
        '\mu_f[:]':            'Fuel Cell Fuel Utilization Coefficient',
        '\eta_{th}[:]':        'Gas Turbine Thermal Efficiency',
        'h_{fuel}':          'Fuel Heating Value',
        'N_p':              'Number of Propulsors',
        'h_{alt}[:]':         'Altitude',
        '(P/m)_{fc}':       'Fuel Cell Specific Power'}

ii = 0 # initializing accumulator
N = 3

while ii < len(sensitivities):
    var_all = variables[ii]
    sens_all = sensitivities[ii]
    
    #sens_sign = np.asarray(sens)/abs(np.asarray(sens))
    
    # sort sens and var based on the sens
    Z = sorted(zip(sens_all, var_all))
    
    var = ['']
    sens = np.array([])
    sens_positive = ['']
    
    kk = 0
    
    while kk < len(Z):
        if Z[kk][0] <= -0.5:
            sens = np.append(sens, (Z[kk][0]))
            var.append(Z[kk][1])
        elif Z[kk][0] >= 0.5:
            sens = np.append(sens, Z[kk][0])
            var.append(Z[kk][1])
            
        kk = kk + 1
    
    sens = np.asarray(sens[sens != 0])
    del var[0]
    
    rr = 0 # initialize counter
    
    Dict_list = list(Dict.items())
    for ll in var:
        for (oo,uu) in Dict_list:
                if var[rr] == oo:
                    var[rr] = uu
        rr = rr + 1 # accumulator
        
    sens_abs = abs(sens)
        
    # remove parameters with sens near zero; separate positive and negative sens
    
    def bar_color(df,color1,color2):
        return np.where(sens>0,color1,color2).T
    
    plt.figure(N)
    plt.barh(np.arange(len(var)), sens_abs, align='center', color=bar_color(sens,'b','r'))
    plt.yticks(np.arange(len(var)), var)
    plt.xlabel('% change in PFEI')
    plt.xlim([0,3.5])
    #plt.title(systems[ii] + ' Parameter Sensitivities')
    
    ii = ii + 1 # accumulator
    N = N + 1
    
# Mission Sweep Chart


    
#print(sol_CHEETA.table())
#print(sol_CHEETA_bat.table())
#print(sol_CHEETA_te.table())
#print(sol_CHEETA_he.table())
#print(sol_737_te.table())
#print(sol_737.table())

#%%
print("Systems")
print(systems)
print("PFEI (kj/kg/km)")
print(performance)
print("MTOW (tonnes)")
print(MTOW)

print('Morph')
print(performance[1])

