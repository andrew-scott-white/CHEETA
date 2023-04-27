

from gpkit import Model, Variable, units

from flight_state import FlightState
from path import C_D_fuse, C_D_wing

import numpy as np
import matplotlib.pyplot as plt

gamma = 1.4
g = 9.81*units('m/s^2')

class HeatExchanger(Model):
    """A component sizing model of the heat exchanger.
    """
    def setup(self):
       # the HEX should be sized for take-off conditions on a hot day.
        A_r = Variable('A_{face}', 'm^2', 'HEX frontal area')
        m = Variable('m_{HEX}', 'kg', 'mass of HEX + plumbing + coolant + pumps')
        psi = Variable('\psi', 'kg/m^2', 'thermal system mass scaling')
        
        constraints = [ 
           m == psi*A_r,
           #A_r <= 1000*units('ft^2'),
           #A_r >= 1*units('in^2')
            ]

        return constraints
    
    def performance(self, flight_state_HEX):
        return HeatExchangerPerformance(self, flight_state_HEX)
        
    
class HeatExchangerPerformance(Model):
    """A component performance model of the heat exchanger.
    """
    def setup(self, HEX, flight_state, VR, T_r):
        
    # Constants
        Cd = Variable('C_d', '-', 'nacelle drag coefficient') # [Robinson et al., 2019]   
        Pr = Variable('Pr', '-', 'HEX Prandtl number')
        K_f = Variable('K_f', '-', 'non-dimensional pressure drop')
        K_h = Variable('K_h', '-', 'non-dimensional heat transfer')
        sigma = Variable('\sigma', '-', 'HEX blockage factor')
        eff = Variable('\epsilon', '-', 'HEX effectiveness')
        
        T_r = Variable('T_r', T_r, 'K', 'radiator temperature')
            
        Q = Variable('Q_{HEX}', 'kW', 'Primary heat dissipated')
                
        # drag parameters
        mdot = Variable('\dot{m}_{HEX}', 'kg/s', 'ram-air HEX flow rate')
        F_HEX = Variable('F_{HEX}', 'kN', 'HEX core heating-induced thrust')
        D_N = Variable('D_{Nacelle}', 'kN', 'nacelle drag')
        D_core = Variable('D_{core}', 'kN', 'HEX core friction drag')
        
        D_tot = Variable('D_{HEX}', 'kN', 'sum of HEX drag sources')
        
        # HEX drag power parameters
        
        VR = Variable('VR', VR, '-', 'HEX velocity ratio')
        delT = Variable('\DelT', 'K', 'HEX temperature differential')
        V_1 = Variable('v_1', 'm/s', 'radiator entrance velocity')
        rho_1 = Variable('\\rho_1', 'kg/m**3', 'air density at radiator face')
        
                
        constraints = [
         VR == V_1/flight_state['V_\infty'],
         
         T_r >= (delT + flight_state['T_{t\\infty}']),
         
         1 >= eff,
         0.000001 <= eff,
         #Q >= 1*units('W'),
         
         mdot >= Q/flight_state['C_p_\infty']/delT/eff, # **
         mdot == rho_1*HEX['A_{face}']*V_1,
         rho_1 == flight_state['rho_\infty'], # assume incompressible flow
         
         D_core*flight_state['V_\infty'] == Q*flight_state['V_\infty']**2 * Pr**(2/3) / flight_state['C_p_\infty']/delT * 1/sigma**2 * (K_f/K_h) * VR**2,
         F_HEX*flight_state['V_\infty'] == Q*.4/2 * flight_state['M_\infty']**2,
         D_N == flight_state['rho_\infty'] * flight_state['V_\infty']**2 * Cd * HEX['A_{face}'], # Basic drag relationship
         
         D_tot >= D_core + D_N
         ]
        

        return constraints
    
    
if __name__ == '__main__':

 # to run this script outside of the aircraft framework, add A_r and Q
 # as input arguments to HeatExchanger and HeatExchangerPerformance above    

    VR_high = 0.3 # cm^2
    VR_low = .05 # cm^2
    VR = np.linspace(VR_low, VR_high, num=30) # sq centimeters
    
    # Cruise Range
    T1 = 217.786 * (1 + 0.5*(gamma-1)*0.773**2) # Ambient Stagnation Temperature, K
    T_r_low = 0.5*T1 + T1 # K
    T_r_high = 1*T1 + T1 # K
    
    # Takeoff Range
    #T1 = 311.15 * (1 + 0.5*(gamma-1)*0.25**2) # Ambient Stagnation Temperature, K
    #T_r_low = 0.1*T1 + T1 # K
    #T_r_high = 1.2*T1 + T1 # K
    
    T_r = np.linspace(T_r_low, T_r_high, num=30) # kW
    
    F_HEX = np.zeros([len(VR),len(T_r)])
    D_nac = np.zeros([len(VR),len(T_r)])
    D_core = np.zeros([len(VR),len(T_r)])
    F_net = np.zeros([len(VR),len(T_r)])

    ii = 0 # initialize A_r counter

   # while ii < len(VR):

    jj = 0 # initialize Q counter

       # while jj < len(T_r):
    HEX = HeatExchanger()
    HEX.substitutions.update({
     '\psi': 5.02*units('kg/ft^2'), # from Boeing (Chellappa)
    })
     
    flight_state = FlightState(C_D_fuse, C_D_wing)
    
    HEXPerf = HeatExchangerPerformance(HEX, flight_state, VR, T_r)
    HEXPerf.substitutions.update({
     
     'Q_{HEX}': 500*units('kW'), # notional heat; since we non-dimensionalize the thrust contribution, this value has no impact
     
     # Design inputs based on [Drela, 1996]
     'K_f': 2, 
     'K_h': 0.93,
     '\sigma': .464,
     
     # Constants
     'Pr': 0.72,
     'C_p_\infty': 1006*units('J/kg/K'),
     'C_d': 0.04,
     '\epsilon': .85,
     
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
     
     # Flight State @ Cruise (20 kft)
     #'M_\infty':           0.7, 
     #'T_{t\\infty}':       248 * (1 + 0.5*(gamma-1)*0.7**2) * units('K'),
     #'V_\\infty':          221 * units('m/s'),
     #'rho_\infty':         .66*units('kg/m**3'),
     
     })
     
    model = Model(HEXPerf['D_{HEX}'], [HEX, HEXPerf])
    sol = model.solve(solver='mosek_conif', verbosity=2, skipsweepfailures=False)
     
          #   jj = jj + 1 # accumulator
         
        # ii = ii + 1 # accumulator
        
F_HEX = np.asarray(sol['variables']['F_{HEX}'])
D_nac = np.asarray(sol['variables']['D_{Nacelle}'])
D_core = np.asarray(sol['variables']['D_{core}'])



F_net = np.asarray(F_HEX - D_nac - D_core)

Q = np.asarray(sol['variables']['Q_{HEX}'])
V_infty = np.asarray(sol['variables']['V_\\infty'])

        
x = np.asarray(sol['variables']['VR'])
y = np.asarray(sol['variables']['\DelT'])/np.asarray(sol['variables']['T_{t\\infty}'])

X = np.reshape(x, (len(VR), len(T_r)))
Y = np.reshape(y, (len(VR), len(T_r)))

Z = np.reshape(F_net*V_infty/(Q), (len(VR), len(T_r)))

delT_low = (T_r_low - np.average(np.asarray(sol['variables']['T_{t\\infty}'])))/(np.average(np.asarray(sol['variables']['T_{t\\infty}'])))
delT_high = (T_r_high - np.average(np.asarray(sol['variables']['T_{t\\infty}'])))/(np.average(np.asarray(sol['variables']['T_{t\\infty}'])))

# formats font
# source: https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
small_size = 11
large_size = 14

cbar_ticks = [-.3,-.25,-.2,-.15, -.1, -.05, 0, .05]
plt.contourf(X, Y, Z, cbar_ticks)
plt.plot(VR, (180+273.15-np.average(np.asarray(sol['variables']['T_{t\\infty}'])))/(np.average(np.asarray(sol['variables']['T_{t\\infty}'])))*np.ones(len(VR)), 'k--') # line indicating CHEETA Fuel Cell Operating Temperature
plt.xticks(np.linspace(VR_low, VR_high, num=int(((VR_high - VR_low)/.05 + 1))))
plt.yticks(np.linspace(delT_low, delT_high,num=6))
plt.xlabel('$V_1$ / $V_\infty$', size=large_size)
plt.ylabel('$(T_r - T_\infty)$ / $T_\infty$', size=large_size)
cbar = plt.colorbar()
cbar.set_label('$\dfrac{\Sigma F_x V_\infty}{Q}$', size=large_size)

print(sol.table())
                
            