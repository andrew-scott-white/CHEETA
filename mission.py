from gpkit import Model, Variable, units, Vectorize, SignomialsEnabled

from path import h_alt, hdot_in, M, P, T, rho, mu, Cp, LDR, C_D_fuse, C_D_wing
from flight_state import FlightState

from aircraft import Aircraft

import numpy as np

gamma = 1.4 # specific heat ratio for air
g = 9.81*units('m/s/s') # gravitational constant
R1 = 287*units('J/kg/K') # gas constant for air
R2 = 8.314*units('J/K') # universal gas constant
F = 96485.33*units('C') # Faraday's constant

class Mission(Model):
    """Breguet Range equation model of an aircraft mission (dynamic).
    """
    def setup(self, aircraft, h_alt, M, P, T, rho, mu, Cp, hdot_in, LDR, C_D_fuse, C_D_wing, BLI):
        f_res = Variable('f_{res}', '-', 'reserve range fraction')
        m_fuel_tot = Variable('m_{fuel,tot}', 'kg', 'total mission fuel mass')
        f_fuel = Variable('f_{fuel}', '-', 'fuel weight fraction')
        h_fuel = Variable('h_{fuel}', 'MJ/kg', 'fuel heating value')
        V_fuel = Variable('V_{fuel}', 'm^3', 'maximum volume of fuel')
        rho_fuel = Variable('\rho_{fuel}', 'kg/m^3', 'density of fuel')
        
        m_pay = Variable('m_{pay}', 'kg', 'payload mass')
        #aircraft['MTOW'] = Variable('MTOW', 'kg', 'aircraft max takeoff mass')
        m_ZF = Variable('m_{ZF}', 'kg', 'aircraft zero-fuel mass')
        m_OE = aircraft['m_{OE}']
        
        R_des = Variable('R_{des}', 'km', 'design range') 
        R_tot = Variable('R_{tot}', 'km', 'total range incl. reserve')
    
        PFEI = Variable('PFEI', 'kJ/kg/km', 'payload-fuel energy intensity')

        # define each flight segment based on the number of altitude points specified
        with Vectorize(len(h_alt)):               
            flight_state = FlightState(C_D_fuse, C_D_wing)
            segment = self.fs = FlightSegment(aircraft, flight_state, h_fuel, BLI)
            
            alt = Variable('h_{alt}', h_alt, 'ft', 'flight altitude')
            hdot = Variable('\dot{h}_{climb}', hdot_in, 'm/s', 'climb rate')
            m_fuel = Variable('m_{fuel}', 'kg', 'fuel consumed per segment')
            
        # define the flight state conditions for each flight segment    
        for ii in range(len(h_alt)):
            flight_state.substitutions.update({
                flight_state['M_\infty'][ii]:           M[ii], 
                flight_state['P_{t\\infty}'][ii]:       P[ii] * (1 + 0.5*(gamma-1)*M[ii]**2)**(gamma/(gamma-1)) * units('Pa'),
                flight_state['P_\\infty'][ii]:          P[ii]*units('Pa'),
                flight_state['T_{t\\infty}'][ii]:       T[ii] * (1 + 0.5*(gamma-1)*M[ii]**2) * units('K'),
                flight_state['T_\\infty'][ii]:          T[ii]*units('K'),
                flight_state['V_\\infty'][ii]:          M[ii]*(gamma*R1*T[ii]*units('K'))**0.5,
                flight_state['rho_\infty'][ii]:         rho[ii]*units('kg/m**3'),
                flight_state['C_p_\infty'][ii]:         Cp[ii]*units('J/kg/K'),
                flight_state['C_{d,wing}'][ii]:         C_D_wing[ii],
                flight_state['C_{d,fuse}'][ii]:         C_D_fuse[ii],
                flight_state['\mu_\infty'][ii]:         mu[ii],
                })
                    
        constraints = [
            flight_state,
            segment,
            h_fuel == segment['HV_{gg}'], # fuel heating value defined in mission model then sent to each segment and subsystem
            segment['L'] == aircraft['MTOW'] * g,  # small angle approximation, assuming cosx = 1 where x is climb angle
            ]
        
        # off-design
        for jj in range(len(h_alt)-1):
            constraints += [
                # the following constraint assumes aircraft weight is constant during climb
                aircraft['MTOW'] >= segment['m_{ac}'][jj] + m_fuel[jj],
                segment['m_{ac}'][jj] >= segment['m_{ac}'][jj+1] + m_fuel[jj], 
                segment['F_{net}'][jj] >= aircraft['MTOW']*g*(hdot[jj]/segment['V_\\infty'][jj] + 1/segment['L/D'][jj]), # small angle approximation of sin(theta) = theta = hdot/V_\infty
                m_fuel[jj] == segment['-\\dot{m}'][jj]*flight_state['dt'][jj], # fuel consumed is equal to fuel consumption rate times duration of segment
                segment['D'][jj] == segment['L'][jj]/segment['L/D'][jj], # definition of lift-to-drag ratio
                ]
            
            with SignomialsEnabled():
                constraints += [
                    flight_state['dt'][jj]*hdot[jj] + alt[jj] >= alt[jj+1], # duration of flight segment is equal to change in altitude divided by climb rate
                    ]
        
        # design point: cruise (the following two lines approximate the Breguet range equation via Taylor expansion)
        x = segment['-\\dot{m}'][len(h_alt)-1] / segment['m_{ac}'][len(h_alt)-1] * R_tot / segment['V_\\infty'][len(h_alt)-1] # neglecting changes in aircraft weight during climb
        exp_x_minus_1 = x + x**2/2 + x**3/6 + x**4/24 + x**5/120 
        
        constraints += [
            # cruise constraints
            segment['m_{ac}'][len(h_alt)-1] >= m_ZF + m_fuel[len(h_alt)-1],
            m_fuel[len(h_alt)-1] >= m_ZF * exp_x_minus_1,
            segment['D'][len(h_alt)-1] == segment['F_{net}'][len(h_alt)-1], # during cruise, thrust is equal to drag
            segment['D'][len(h_alt)-1] == segment['L'][len(h_alt)-1]/segment['L/D'][len(h_alt)-1],
            flight_state['dt'][len(h_alt)-1] == R_des/flight_state['V_\infty'][len(h_alt)-1],
            ]
        
        # neglect fuel consumption duing descent, landing, and taxi
        
       # mission and performance constraints
        constraints += [
             # aircraft performance constraints
             V_fuel == m_fuel_tot / rho_fuel,
             f_fuel == m_fuel_tot/aircraft['MTOW'],
             m_ZF >= m_OE + m_pay, # zero-fuel mass
             m_fuel_tot >= m_fuel[0] + m_fuel[1] + m_fuel[2] + m_fuel[3] + m_fuel[4] + m_fuel[5] + m_fuel[6] + m_fuel[7] + m_fuel[8] + m_fuel[9] + m_fuel[len(h_alt)-1],
            
             # mission-related sizing constraints 
             aircraft['MTOW'] >= m_fuel_tot + m_ZF, # take-off mass
             R_tot == R_des * f_res, # range times 15% for fuel reserves
             
             PFEI == m_fuel_tot*h_fuel/R_des/m_pay, # performance parameter
            ]

        return constraints
    

         
class FlightSegment(Model):
     """flight segment model that defines aircraft performance throughout flight
     """
     def setup(self, aircraft, flight_state, h_fuel, BLI):
         
        m_fuel = Variable('m_{fuel}', 'kg', 'fuel mass')
        #dt = Variable('dt', 's', 'time')
        m_aircraft = Variable('m_{ac}', 'kg', 'aircraft weight')
        
        dt = flight_state['dt']
        
        
        self.aircraft = aircraft
        self.performance = ac_perf = aircraft.performance(flight_state, h_fuel, BLI)
        
        constraints = [           
            ac_perf
        ]
        
        
        return constraints
    


# Turbofan configuration --------------------------------------------------------------------------------------------------------------                    

if __name__ == '__main__':
    from turbofan import Turbofan
               
    BLI = 'on'
    
    # set up propulsion system
    tf = Turbofan(2)
    tf.substitutions.update({
        # gas generator constants
        '(P/m)_{gg}': 7687 * units('hp') / 3065 / units('lb'), # half of LEAP-1B
        # ducted fan constants
        'K_{fan}': 0.5 * 6130 * units('lb') / (69.4 * units('in'))**3, # half of LEAP-1B
        'HTR': 0.3,
    })

    # set up aircraft system (to which propulsion system is passed as input)
    ac = Aircraft(tf)
    ac.substitutions.update({
        'f_{empty}': 0.5  , # assumed
        'd_{fan}': 69.4*units('in'), # LEAP-1B engine fan diameter
        })
    
    mission = Mission(ac, h_alt, M, P, T, rho, mu, Cp, hdot_in, LDR, C_D_fuse, C_D_wing, BLI)
    # all aircraft and propulsion system performance models (and thus their
    # variables) get set up inside Mission(), so we update their constants here
    
    mission.substitutions.update({
        '\rho_{fuel}':              840*units('kg/m^3'), # jet fuel density
        'h_{fuel}':                43 * units('MJ/kg'), # jet fuel lower heating value
        
    # mission parameters
        'f_{res}': 1.15, # assume 15% additional fuel for reserves
        'R_{des}': 2935 * units('nmi'), # design range
        'm_{pay}':  35000 * units('lb'), # design payload
    })
    
    mission.substitutions.update({
        'D_corr':           0.5137022627, # D(M_fan=0.66) corrected flow per unit area, https://asmedigitalcollection-asme-org.libproxy.mit.edu/turbomachinery/article/129/1/184/474859/Preliminary-Fan-Design-for-a-Silent-Aircraft
        'L/D':              LDR, # take input lift-to-drag ratios from Elias
        #'f_{BLI}':                   0.0001, # assuming no BLI
        #'1 - f_{wake}':              0.868, # from D8, [Hall et al, 2017]; according to D.K. Hall, this number is generally representative         
        #'D_p':                       4830*units('lbf'), # profile drag ingested by propulsors
        '\\eta_{th}':                0.55,
        #'\\eta_{fan}':               0.95, 
        #'Cd':                        0.04, # [Robinson, 2017]
        
        })
    
    model = Model(mission['m_{fuel,tot}'], [ac, mission])

    sol = model.localsolve()
    print(sol.summary())
    
