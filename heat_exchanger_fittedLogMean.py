from path import C_D_fuse, C_D_wing

from gpkit import Model, Variable, units, SignomialsEnabled

from flight_state import FlightState

import numpy as np
import matplotlib.pyplot as plt

gamma = 1.4
g = 9.81*units('m/s^2')

class HeatExchanger(Model):
    """A component sizing model of the heat exchanger (assumed to be counterflow, cylindrical tube with fins)
    """
    def setup(self):
       # the HEX should be sized for take-off conditions on a hot day.
        m = Variable('m_{HEX}', 'kg', 'mass of HEX + plumbing + coolant + pumps')
        psi = Variable('\psi', 'kg/m^2', 'thermal system mass scaling')
        sigma = Variable('\sigma', '-', 'HEX blockage factor (ratio of freeflow to frontal area)')
      #  beta = Variable('\\beta', '-', '1 - \sigma')
      #  rho_HEX = Variable('\rho_{HEX}', 'kg/m^3', 'density of HEX material')
      #  Vol_HEX = Variable('Vol_{HEX}', 'm^3', 'volume of HEX material')
        
        # HEX dimensions
        A_r = Variable('A_{face}', 'm^2', 'HEX frontal area')
        A_w = Variable('A_{surf}', 'm^2', 'HEX (wetted) surface area')
        
        #N_plate = Variable('N_{plates}', 60, '-', 'number of plates on cold side of HEX')
        #d_plate = Variable('d_{plate}', .25, 'mm', 'thickness of plate')
        #w_plate = Variable('w_{plate}', 'cm', 'width of plate')
        
        
        r_h = Variable('r_h', 'cm', 'HEX hydraulic radius')
        L_c = Variable('L_c', 7, 'cm', 'length of cold-side flow path') # value from Chellappa
        
        
        constraints = [ 
           m == psi*A_r,
           r_h == L_c*sigma*A_r/A_w,
            ]

        return constraints
    
    def performance(self, flight_state_HEX):
        return HeatExchangerPerformance(self, flight_state_HEX)
        
    
class HeatExchangerPerformance(Model):
    """A component performance model of the heat exchanger.
    """
    def setup(self, HEX, flight_state):
        
    # Constants
        Cp_c = Variable('C_{p,c}', 'J/kg/K', 'Specific heat of cold-side fluid (air)')
        Cp_h = Variable('C_{p,h}', 'J/kg/K', 'Specific heat of hot-side fluid (coolant)')
            
        # performance parameters
        Q = Variable('Q_{HEX}', 'kW', 'Primary heat dissipated')
        K_f = Variable('K_f', '-', 'non-dimensional pressure drop')
        K_h = Variable('K_h', '-', 'non-dimensional heat transfer')
        h = Variable('h', 'W/m^2/K', 'average heat transfer coefficient')
        C_f = Variable('C_f', '-', 'radiator skin friction factor')
        P = Variable('\phi', '-', 'radiator pressure drop parameter')
        z = Variable('\zeta^2', '-', 'HEX non-dimensional length parameter')
        St = Variable('St', '-', 'Stanton number')
        Pr = Variable('Pr', '-', 'Prandtl number')
        eff = Variable('\epsilon', '-', 'HEX effectiveness')
        
        mdot_c = Variable('\dot{m}_c', 'kg/s', 'cold-side (air) flow rate')
        mdot_h = Variable('\dot{m}_h', 'kg/s', 'hot-side (coolant) flow rate')
       
        Tc_in = Variable('T_{c,in}', 'K', 'Cold-side (air) inlet temperature')
        Tc_out = Variable('T_{c,out}', 'K', 'Cold-side (air) exit temperature')
        Th_in = Variable('T_{h,in}', 'K', 'Hot-side (coolant) inlet temperature')
        Th_out = Variable('T_{h,out}', 'K', 'Hot-side (coolant) exit temperature') 
       
        delT = Variable('\DelT', 'K', 'HEX max temperature differential')
        delT_1 = Variable('\DelT_1', 'K', 'Temperature difference between cold-side inlet and hot side outlet')
        delT_2 = Variable('\DelT_2', 'K', 'Temperature difference between cold-side outlet and hot side inlet')
        delT_c = Variable('\DelT_c', 'K', 'temperature difference between cold-side outlet and inlet')
        delT_h = Variable('\DelT_h', 'K', 'temperature difference between hot-side inlet and outlet')
        LMdelT = Variable('\DelT_{LM}', 'K', 'HEX log-mean temperature differential')
                
        # drag parameters
        delP = Variable('\DelP', 'Pa', 'HEX Pressure Drop')
        D_core = Variable('D_{core}', 'kN', 'HEX core friction drag')
        F_HEX = Variable('F_{HEX}', 'kN', 'HEX core heating-induced thrust')
        D_N = Variable('D_{Nacelle}', 'kN', 'nacelle drag')
        Cd = Variable('C_d', '-', 'nacelle drag coefficient') # [Robinson et al., 2019]          
        D_tot = Variable('D_{HEX}', 'kN', 'sum of HEX drag sources')
        
        # HEX drag power parameters
        
        VR = Variable('VR', '-', 'HEX velocity ratio')
        
        V_1 = Variable('v_1', 'm/s', 'radiator entrance velocity')
        rho_1 = Variable('\\rho_1', 'kg/m**3', 'air density at radiator face')
        
        rho_1 = flight_state['rho_\infty'] # assume incompressible flow
        Cp_c = flight_state['C_p_\infty'] # cold-side specific heat is that of freestream air
        Tc_in = flight_state['T_{t\\infty}'] # cold-side inlet temperature is equal to ambient stagnation temperature
        mu_1 = flight_state['\mu_\infty'] # cold-side dynamic viscosity
                
        constraints = [
         Th_in >= (delT + Tc_in), # define max temperature difference
         
         # Temperature differences
         Th_in >= delT_2 + Tc_out,
         Th_out >= delT_1 + Tc_in,
         Tc_out >= delT_c + Tc_in,
         Th_in >= delT_h + Th_out,
         
         # Temperature limits based on second law of thermodynamics
         Th_in >= Tc_out, # cold-side cannot reach a higher temperature than the incoming heat
         Th_out >= Tc_in, # hot-side cannot reach a lower temperature than incoming air
         #Tc_out >= Tc_in, # cold-side will heat up during exchange
         #Th_in >= Th_out, # hot-side will cool down during exchange
    
        # Log-mean temperature difference fit (RMS error of 2.2%; max error of 3.9%)
         (LMdelT/253/units('K'))**0.1 == 1.07788*(delT_1/253/units('K'))**0.0370743*(delT_2/253/units('K'))**0.0599609,
         #eff == delT_c/delT, # assuming mdot_c*Cp_c < mdot_h*Cp_h
         
        # Log-mean temperature difference fit (RMS error of 0.0%; max error of 0.0%)
        #(LMdelT/253/units('K'))**0.549219 >= 0.777845*(delT_1/253/units('K'))**0.492314*(delT_2/253/units('K'))**0.056878 + 0.803171*(delT_1/253/units('K'))**0.0638169*(delT_2/253/units('K'))**0.485421,

         #mdot_c >= Q/Cp_c/delT/eff, # required airflow for cooling
         #Q == z*K_h*Cp_c*delT*Pr**(-1/6)*mdot_c,
         mdot_c == rho_1*HEX['A_{face}']*V_1,
         Q == h*HEX['A_{surf}']*LMdelT,
         Q == mdot_c*Cp_c*delT_c,
         Q == mdot_h*Cp_h*delT_h,
                  
         VR == V_1/flight_state['V_\infty'], # define velocity ratio
         
         P == (C_f/HEX['\sigma']**2)*HEX['L_c']/HEX['r_h'],
         C_f == HEX['\sigma']*HEX['A_{face}']/HEX['A_{surf}']*2*delP/rho_1/(V_1/HEX['\sigma'])**2,

         St == Q/HEX['A_{surf}']/Cp_c/delT/rho_1/V_1*HEX['\sigma'],
         delP == rho_1*V_1**2*P/2,
         
         K_f == C_f/2*(rho_1*V_1*HEX['L_c']/HEX['\sigma']/mu_1)**(1/2),
         K_h == St*Pr**(2/3)*(rho_1*V_1*HEX['L_c']/HEX['\sigma']/mu_1)**(1/2),
         
         10 == K_f/K_h/HEX['\sigma']**2,
         
     # Fit to Drela data
         # definition of non-dimensional HEX length
         z == HEX['L_c']/HEX['r_h']*HEX['\sigma']/Pr*mu_1/rho_1/V_1/HEX['r_h'],
         
         # limits based on [Drela, 1996] plots
         1 >= z,
         z >= 0.001,
         K_f >= 0.664, # flat plate
         3 >= K_f,
         K_h >= 0.664, # flat plate
         1 >= K_h,
         
         
         D_core*flight_state['V_\infty'] == Q*flight_state['V_\infty']**2 * Pr**(2/3) / Cp_c/delT * 1/HEX['\sigma']**2 * (K_f/K_h) * VR**2,
         F_HEX*flight_state['V_\infty'] == Q*.4/2 * flight_state['M_\infty']**2,
         D_N == rho_1 * flight_state['V_\infty']**2 * Cd * HEX['A_{face}'], # Basic drag relationship
         
         D_tot >= D_core + D_N + HEX['m_{HEX}']*g/18.4,
         ]
        
        with SignomialsEnabled():
            constraints += [
                K_f + 1.4314*z**2 >= 3.2323*z + 0.792,
                K_h + 0.6934*z**2 <= 0.815*z + 0.7258,
                ]
        

        return constraints
    
    
if __name__ == '__main__':

    HEX = HeatExchanger()
    HEX.substitutions.update({
     '\psi': 5.02*units('kg/ft^2'), # from Boeing (Chellappa)
     #'\sigma': .464,
     #'L_c': 1*units('m'),
     #'\rho_{HEX}': 2710*units('kg/m^3'), # density of aluminum
     
     
    })
     
    flight_state = FlightState(C_D_fuse, C_D_wing)
    
    HEXPerf = HeatExchangerPerformance(HEX, flight_state)
    HEXPerf.substitutions.update({
     
     'Q_{HEX}': 300*units('kW'), # notional heat; since we non-dimensionalize the thrust contribution, this value has no impact
     
     # Design inputs based on [Drela, 1996]
     #'K_f': 2, 
     #'K_h': 0.93,
     #'\phi': 20,
     'h': 25*units('W/m^2/K'),

     
     # Constants
     'Pr': 0.72,
     'C_p_\infty': 1006*units('J/kg/K'),
     'C_{p,h}': 2200*units('J/kg/K'),
     'C_d': 0.04,
     #'\epsilon': .85,
     '\dot{m}_h': 15*units('kg/s'),
     'T_{h,in}': (180+273.15)*units('K'),
     
     # Flight State @ Take-off
     #'T_{t\\infty}':       311.15 * (1 + 0.5*(gamma-1)*0.25**2) * units('K'),
     #'V_\\infty':          88.4*units('m/s'),
     #'M_\infty':           0.25,
     #'rho_\infty':         1.135*units('kg/m**3'),
     
     # Flight State @ Cruise (37.5 kft)
     'M_\infty':           0.773, 
     'T_{t\\infty}':       217.786 * (1 + 0.5*(gamma-1)*0.773**2) * units('K'),
     'V_\\infty':          443.369 * units('kts'),
     'rho_\infty':         .3509518*units('kg/m**3'),
     '\\mu_\infty':        0.00001447704*units('Pa*s'),
     
     # Flight State @ Cruise (20 kft)
     #'M_\infty':           0.7, 
     #'T_{t\\infty}':       248 * (1 + 0.5*(gamma-1)*0.7**2) * units('K'),
     #'V_\\infty':          221 * units('m/s'),
     #'rho_\infty':         .66*units('kg/m**3'),
     
     })
     
    init = {
            HEXPerf['L_c']: 3*units('m'),
            HEXPerf['K_f']: 2, 
            HEXPerf['K_h']: .9, 

        }
    
    model = Model(HEXPerf['D_{HEX}'], [HEX, HEXPerf])
    sol = model.localsolve(solver='mosek_conif', verbosity=2, skipsweepfailures=False)

print(sol.table())
                
            