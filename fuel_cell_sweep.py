import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from path import C_D_fuse, C_D_wing

from gpkit import Model, Variable, units, SignomialsEnabled
from flight_state import FlightState

R2 = 8.314*units('J/K')
F = 96485.33*units('C')
g = 9.81*units('m/s/s')
gamma = 1.4

class FuelCell(Model):
    """A component sizing model of the fuel cell.
    """
    def setup(self):
        m_fc = Variable('m_{fc}', 'lb', 'stack mass')
        P_max = Variable('P_{max,fc}', 'kW', 'fuel cell stack max power')
        P_m = Variable('(P/m)_{fc}', 'W/kg', 'fuel cell stack specific power')
        m_comp = Variable('m_{comp}', 'lb', 'stack compression system weight')
        N_c = Variable('N_{cells}', '-', 'number of cells in fuel cell stack')
        A = Variable('A', 'cm**2', 'cell membrane area')
        
        constraints = [
            m_fc == P_max/P_m,
            #1.5*units('MW') >= P_max,
            #P_max >= .75*units('MW')
            ]

        return constraints
    
    def performance(self, flight_state):
        return FuelCellPerformance(self, flight_state)

class FuelCellPerformance(Model):
    """A component performance model of the fuel cell.
    """
    def setup(self, fuelcell, flight_state_fc, i_d):
        
        c = Variable('c', '-', 'oxygen concentration in air')
        
        # Device design variables
        i_n = Variable('i_{n}', 'A', 'fuel cell internal current')
        i = Variable('i', 'A', 'fuel cell output current')
        e = Variable('e', 'V', 'fuel cell (Nernst) electromotive force')
        e0 = Variable('e0', 'V', 'reversible reference Nernst voltage')
        Vc = Variable('Vc', 'V', 'real cell voltage')
        V = Variable('V', 'V', 'fuel cell stack output voltage')
        ohm_n = Variable('\\ohm_n', 'ohm*cm**2', 'cell internal resistance-area')
        
        mdot_fuel = Variable('\\dot{m}_{f}', 'g/s','anode-side fuel flow rate')
        mdot_air = Variable('\\dot{m}_{a}', 'g/s', 'cathode-side air flow rate')
        Ptf = Variable('P_{t\\fuel}', 'kPa', 'fuel pressure')
        HV = Variable('HV_{fc}', 'MJ/kg', 'fuel heating value of fuel')
        mu_f = Variable('\mu_f', '-', 'fuel utilization')
               
        # Performance variables
        P_out = Variable('P_{out}', 'kW', 'output power of fuel cell system')
        P_gross = Variable('P_{gross}', 'kW', 'gross power of fuel cell stack')
        Q_ohm = Variable('Q_{ohm}', 'kW', 'ohmic power losses')
        Q_act = Variable('Q_{act}', 'kW', 'activation power losses')
        Q_mass = Variable('Q_{mass}', 'kW', 'mass transport power losses')
        delG = Variable('\\DeltaG', 'kJ', 'Gibbs energy of formation')
        Q = Variable('Q_{HEX}', 'kW', 'heat rejected by the fuel cell')
       
        S = Variable('S', '-', 'air stoichiometric constant')
        pi_c = Variable('pi_c', '-', 'stack compressor pressure ratio')
        P_comp = Variable('P_{comp}', 'kW', 'power required to compress air')
        eta_c = Variable('eta_{comp}', '-', 'fuel cell system compressor efficiency')
        P_para = Variable('P_parasitic', 'kW', 'parasitic losses for coolant pump and H2 recycler')
        x_para = Variable('x_{parasitic}', '-', 'parasitic fraction of power')
        
        eta_fc = Variable('eta_{fc}', '-', 'fuel cell efficiency at cruise')
        eta_net = Variable('eta_{net}', '-', 'net power system efficiency at cruise')
        V_ohm = Variable('V_{\ohm}', 'V', 'ohmic voltage drop per cell')
        V_mass = Variable('V_{mass}', 'V', 'mass transport voltage drop per cell')
        m = Variable('m','V', 'empirical constant')
        n = Variable('n','cm**2/mA', 'empirical constant')
    
        # Activation loss variables
        V_act = Variable('V_{act}', 'V', 'activation voltage drop per cell') 
        alpha = Variable('alpha', '-', 'charge transfer coefficient (empirical)')
        i_d = Variable('(i/A)', i_d, 'A/cm**2', 'fuel cell internal current density')
        i0_d = Variable('(i0/A)', 'A/cm**2', 'cathode exchange current density')
        k_r = Variable('k', '1/s', 'electrochemical rate constant')
                  
        # Taylor expansion dummy variable for activation losses
        z = 2*alpha*F*(V_act)/R2/flight_state_fc['T_{t\\2}']

     
        constraints = [
            fuelcell, 
                        
           # Sum of voltages per cell 
            e >= Vc + V_ohm + V_act + V_mass, # sum of voltages
            V == P_gross/i, 
            V == Vc * fuelcell['N_{cells}'], 
            
            i_d == i_n/fuelcell['A'], # define internal current density
                
           # ideal gas usage current
               i_n == 2*F*mdot_fuel/2/.0010079/units('kg')/fuelcell['N_{cells}'], # for testing with Pratt paper (flow rates given at STP)
               i == i_n*mu_f,
            
            # Low-temp GPfit for Nernst (80°C)
                #(flight_state_fc['P_{t\\2}']/2000/units('kPa'))**0.1 == 0.0385288 * (e/1.48/units('V'))**(-15.2633) * (Ptf/2000/units('kPa'))**0.2,
                
            # High-temp GPfit for Nernst (170°C)
                #(flight_state_fc['P_{t\\2}']/2000/units('kPa'))**0.1 == 0.0541579 * (e/1.48/units('V'))**(-12.1878) * (Ptf/2000/units('kPa'))**0.2,
                
            # Higher-temp GPfit for Nernst (180°C) - Mean-squared error of 1.6%; MAX error 4%
                (flight_state_fc['P_{t\\2}']/2000/units('kPa'))**0.1 == 0.0572938 * (e/1.48/units('V'))**-11.8829 * (Ptf/2000/units('kPa'))**0.200026,
                
            # Ohmic voltage drop
            V_ohm == i_d*ohm_n, # Ohm's law
            
            # Mass transport voltage drop
            V_mass >= m*(1 + n*i_d + (n*i_d)**2/2 + (n*i_d)**3/6 + (n*i_d)**4/24 + (n*i_d)**5/120 + (n*i_d)**6/720 + (n*i_d)**7/5040 + (n*i_d)**8/40320),
            ]
        
        with SignomialsEnabled():
            constraints += [
            # Activation voltage drop
                i_d/i0_d <= 1 + z + z**2/2 + z**3/6 + z**4/24 + z**5/120 + z**6/720 + z**7/5040 + z**8/40320 + z**9/362880 + z**10/3628800 + z**11/39916800 + z**12/479001600 \
                    + z**13/6227020800 + z**14/8.71782912e10 + z**15/1.307674368e12 + z**16/2.092278989e13 + z**17/3.556874281e14 + z**18/6.402373706e15 + z**19/1.216451004e17 + z**20/2.432902008e18, # (9) Tafel equation
                     
            # Stack compression system
                    P_comp >= mdot_air*flight_state_fc['C_p_\infty']*flight_state_fc['T_{t\\infty}']*(pi_c**(.4/1.4)-1)/eta_c,
                    ]
                
        constraints += [
                
            # Activation voltage drop cont...
                i0_d == 2*F*k_r*c*flight_state_fc['P_{t\\2}']/(101.325*units('kPa'))/fuelcell['A'], # (10) exchange current estimate (Pratt, JPP 2007)
            
            # heat loss mechanisms
                Q_ohm == i_d**2*ohm_n*fuelcell['N_{cells}']*fuelcell['A'], # stack ohmic power losses
                Q_act == i_d*fuelcell['A']*V_act*fuelcell['N_{cells}'], # stack activation power losses
                Q_act >= i0_d*fuelcell['A']*V_act*fuelcell['N_{cells}'], # minimum activation power losses
                Q_mass == i_d*fuelcell['A']*V_mass*fuelcell['N_{cells}'], # mass transport power losses
                    
            # air flow rate
                mdot_air == i_n*fuelcell['N_{cells}']/4/F * 28.9647*units('g')/c * S,
                
                pi_c == flight_state_fc['P_{t\\2}']/flight_state_fc['P_{t\\infty}'],
                
                # Fit from Boeing
                 flight_state_fc['P_{t\\2}'] >= 1.125*units('bar') + .175*units('cm**4*bar/A**2') * i_d**2,
                
       
            # fuel cell performance
                eta_net == P_out/(HV*mdot_fuel),
                eta_fc == P_gross/(HV*mdot_fuel),
                e0 == delG/2/F,            
                
            # "balance of plant" 
                mdot_fuel*HV >= Q + P_gross,
                P_gross >= P_out + P_comp + P_para, # Assume all unused power is heat generated
                P_para == x_para*P_gross, 
                Q >= Q_ohm + Q_act + Q_mass,
                fuelcell['P_{max,fc}'] >= P_out,
                
                #3.5*units('A/cm**2') >= i_d, # limiting current density
                     ]
        
        return constraints

         
if __name__ == '__main__':
    
    i_d_low = 0.01 # A/cm^2
    i_d_high = 4.3 # A/cm^2
    #i_d = np.linspace(i_d_low, i_d_high, num=200) # A/cm^2
    i_d_1 = np.asarray([.02,.03,.04,.05,.075,.1,.125,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9])
    #i_d_1 = np.asarray([.05,.5])
    i_d_2 = np.linspace(1,i_d_high,num=40)
    i_d = np.concatenate((i_d_1, i_d_2))
    #i_d = [2.5,3]

    
    fuelcell = FuelCell()
    fuelcell.substitutions.update({
        'P_{max,fc}': 850*units('kW')
        })
    flight_state_fc = FlightState(C_D_fuse, C_D_wing)
    fcPerf = FuelCellPerformance(fuelcell, flight_state_fc, i_d)
    fcPerf.substitutions.update({
    
    # Inputs       
    
         #'P_{gross}':       850        *units('kW'),
         #'P_{out}':             700         *units('kW'),
        
        # Take-off flight state
         #'P_{t\\infty}':       101325 * (1 + 0.5*(gamma-1)*0.25**2)**(gamma/(gamma-1)) * units('Pa'),
         #'T_{t\\infty}':       311.15 * (1 + 0.5*(gamma-1)*0.25**2) * units('K'),
         #'C_p_\infty':         1006*units('J/kg/K'),
          
        # Cruise flight state (altitude = 37.5 kft)
            'P_{t\\infty}':       21876 * (1 + 0.5*(gamma-1)*0.773**2)**(gamma/(gamma-1)) * units('Pa'),
            'T_{t\\infty}':       217.786 * (1 + 0.5*(gamma-1)*0.773**2) * units('K'),
            'C_p_\infty':         1006*units('J/kg/K'),
            
        # Cruise flight state (altitude = 20 kft)
              #  'P_{t\\infty}':       46540 * (1 + 0.5*(gamma-1)*0.773**2)**(gamma/(gamma-1)) * units('Pa'),
              #  'T_{t\\infty}':       248.15 * (1 + 0.5*(gamma-1)*0.773**2) * units('K'),
               # 'C_p_\infty':         1006*units('J/kg/K'),

    # low-temperature fuel cell parameters
        #'HV':                   141.8*units('MJ/kg'),
        #'T_{t\\2}':             (90+273.15)     *units('K'),
        #'\\Delta G':            228.2*units('kJ'),
        #'k':                    1.75e-7          *units('1/s'),  
        #'\\ohm_n':              .038            *units('ohm*cm**2'),
        #'x_{parasitic}':          .02,
        #'\mu_f':                .95,
        
    # high-temperature fuel cell parameters
        'HV_{fc}':              120.0*units('MJ/kg'),
        'T_{t\\2}':             (180+273.15)     *units('K'),
        '\\DeltaG':            224.24*units('kJ'),
        'k':                    5e-5          *units('1/s'),  
        '\\ohm_n':              .062            *units('ohm*cm**2'),
        'x_{parasitic}':        .05,
        '\mu_f':                .96,
        'm':                    4e-5*units('V'), # [Larminie and Dicks, Fuel Cell Systems Explained]
        'n':                    2.3e-3*units('cm**2/mA'), # Based on Larminie, but relaxed
        #'n':                    8e-3*units('cm**2/mA'), # Larminie
        
    # Design inputs
    
        # Fuel Cells
        'P_{t\\fuel}':          5        *units('bar'),
        #'P_{t\\2}':          3        *units('bar'),
        'alpha':                .5, 
        'c':                    .21,  # assumed mole ratio of O2 in dry air
        '(P/m)_{fc}':           2700*units('W/kg'), # CHEETA Data Sheet
        'A':                    446*units('cm**2'),
        'N_{cells}':            1080,
        'S':                    1.5,
        #'k_m':                  5.12e-11*units('1/bar/cm/s'), # permeability of H2 in Nafion at 80°C
        'P_{max,fc}':           850*units('kW'), # max power the model is accurate for. This could be replaced with a max power specified by Boeing
        
        # Compressor
        #   'pi_c':                 2.5,
        'eta_{comp}':           .9,
             
        })
        
           
    #model = Model(fcPerf['\\dot{m}_{f}'], [fcPerf])
    model = Model(1/fcPerf['eta_{net}'], [fcPerf])
    sol = model.localsolve(solver='mosek_conif', verbosity=1)
    print(sol.table())
    
i_d = sol['variables']['(i/A)']
V = sol['variables']['Vc']
eta_fc = sol['variables']['eta_{fc}']
eta_net = sol['variables']['eta_{net}']
P_out = sol['variables']['P_{out}']
P_gross = sol['variables']['P_{gross}']
P_max = sol['variables']['P_{max,fc}']
P_sp = P_gross/1000
    
plt.figure(1)
line1, = plt.plot(i_d, V, 'k-', label='Fuel Cell Voltage (V)')
line2, = plt.plot(i_d, eta_fc, 'k--', label='Fuel Cell Efficiency')
#line3, = plt.plot(i_d, eta_net, 'k:', label='Stack Efficiency')
line4, = plt.plot(i_d, P_sp, 'k:', label='Output Power Fraction')
plt.xlabel("Internal Current Density (A/$cm^2$)")
plt.ylabel("")
plt.ylim([0,1.2])
plt.xlim([i_d_low, 4.5])
plt.legend(handles=[line1, line2, line4], loc='lower right')
    
