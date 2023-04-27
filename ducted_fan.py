from gpkit import Model, Variable, units, SignomialsEnabled
from numpy import pi

R = 287*units('J/kg/K')
gamma = 1.4

class DuctedFan(Model):
    """A component sizing model of the propulsor.
    """
    def setup(self):
        A_fan = Variable('A_{fan}', 'm^2', 'fan area')
        d_fan = Variable('d_{fan}', 'm', 'fan diameter')
        HTR = Variable('HTR', '-',
            'fan hub-to-tip radius ratio')
        K_fan = Variable('K_{fan}', 'kg/m**3', 'fan mass scaling constant') # m_ref / d_ref**2.4
        m_fan = Variable('m_{fan}', 'kg', 'mass of one propulsor')
        r_hub = Variable('r_{hub}', 'm', 'fan hub radius')
        r_tip = Variable('r_{tip}', 'm', 'fan tip radius')
    
        constraints = [
            #m_fan / (6130*.5) /units('lb') >= (2 * r_tip / 69.4/units('in'))**3, # sizing scaled to a fan weight equal to half of a LEAP-1B engine
            m_fan >= K_fan * d_fan**3, # mass scaling relationship
            d_fan == r_tip * 2, 
            r_tip == r_hub / HTR, 
            r_tip**2 >= r_hub**2 + A_fan / pi, # duct area relationship
            #14*units('ft') >= 3*d_fan, # fan diameter constrained by width of CHEETA tail
            ] 

        return constraints

    def performance(self, flight_state, f_BLI, f_wake, D_p, BLI):
        return DuctedFanPerformance(self, flight_state, f_BLI, f_wake, D_p, BLI)
    
    
class DuctedFanPerformance(Model):
    """A component performance model of the propulsor.
    """
    def setup(self, fan, flight_state, f_BLI, f_wake, D_p, BLI):
        
        # constants
        Cd = Variable('Cd', .04, '-', 'nacelle drag coefficient')
        eta_fan = Variable('\\eta_{fan}', .93, '-', 'fan efficiency')
        D = Variable('D_corr', 0.5137, '-', 'corrected flow per unit area')
        #D_p = Variable('D_p', 'kN', 'profile drag')
        K_inl = Variable('K_{inl}', 'kW', 'propulsor inlet kinetic energy defect')
        
        # cruise parameters
        F_t = Variable('F_{fan}', 'kN', 'thrust')
        F_d = Variable('F_{nac}', 'kN', 'nacelle drag')
        mdot = Variable('\\dot{m}', 'kg/s', 'fan mass flow')
        P_shaft = Variable('P_{shaft}', 'MW', 'fan shaft power')
        P_K = Variable('P_K', 'MW', 'fan mechanical flow power')
        u_jet = Variable('u_{jet}', 'm/s', 'jet perturbation velocity')
        
        A_fan = fan['A_{fan}']
        
        # cruise flight state parameters
        pt0 = flight_state['P_{t\\infty}']
        Tt0 = flight_state['T_{t\\infty}']
        V0 = flight_state['V_\\infty']
        rho0 = flight_state['rho_\\infty']
        
        # nacelle drag coefficient
        C = pi/8*Cd*rho0*V0**2
        
       
        constraints = [
            P_shaft == P_K / eta_fan,
            
            P_K >= 0.5 * mdot * u_jet * (2 * V0 + u_jet) + K_inl, # mechanical power to the flow = energy flux of jet + dissipation upstream of the inlet
            
            # fan sizing constraints for cruise
            A_fan >= mdot * (R * Tt0 / gamma)**0.5 / D / pt0, # corrected flow per unit area
            F_d == C*fan['d_{fan}']**2, # nacelle drag constraint Cd = d**2, D = C * Qinf (from flight_state)
            K_inl == f_BLI * f_wake * D_p * flight_state['V_\\infty'], # inlet kinetic energy defect
        ]        
        
        if BLI == 'on':
            # cruise constraints
            with SignomialsEnabled():
                constraints += [
                    mdot * u_jet + f_BLI * D_p >= F_t + F_d, # thrust force balance
                    ]            
        else:
            constraints += [
                    mdot * u_jet >= F_t + F_d, # thrust force balance
                    #K_inl == 1e-12*units('kW'),
            ]  
            
        return constraints