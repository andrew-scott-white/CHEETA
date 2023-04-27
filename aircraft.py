from gpkit import Model, Variable, units, SignomialsEnabled

class Aircraft(Model):
    """Aircraft sizing (static) model.
    """
    def setup(self, propulsion_system):
        m_MTO = Variable('MTOW', 'kg', 'aircraft max take-off mass')
        m_OE = Variable('m_{OE}', 'kg', 'operating empty mass')
        f_empty = Variable('f_{empty}', '-',
            'aircraft empty mass fraction (excl. propulsion)')
        b_fuse = Variable('b_{fuse}', 'm', 'fuselage span')
        b_wing = Variable('b_{wing}', 'm', 'wing span')
        S_ref = Variable('S_{ref}', 'm^2', 'reference area')

        self.propulsion_system = propulsion_system
        # expects propulsion_system a Model that has a variable with the string
        # name 'm_{\\rm prop.~sys.}' which is the total propulsion system mass
        # and performnace method with a variable names named '-\\dot{m}' which
        # is the total propulsion system fuel consumption rate and another
        # named 'F_{\\rm net}' which is the total propulsion system thrust
        
        constraints = [
            propulsion_system,
            m_OE >= propulsion_system['m_{propsys}'] \
                + f_empty * m_MTO,
            b_wing >= 0.5*propulsion_system['N_p']*2/3*propulsion_system['d_{fan}'],
            b_fuse >= propulsion_system['N_p']/3*propulsion_system['d_{fan}'],
                
                 m_OE >= 1*units('lb'),
                 m_OE <= 1e12*units('lb'),
        ]

        return constraints

    def performance(self, flight_state, h_fuel, BLI):
        return AircraftPerformance(self, flight_state, h_fuel, BLI)


class AircraftPerformance(Model):
    def setup(self, aircraft, flight_state, h_fuel, BLI):
        D = Variable('D', 'kN', 'drag')
        L = Variable('L', 'kN', 'lift')
        LDR = Variable('L/D', '-', 'lift-to-drag ratio')
        
        D_p_fuse = Variable('D_{p,fuse}', 'N', 'fuselage profile drag')
        D_p_wing = Variable('D_{p,wing}', 'N', 'wing profile drag')
       
        b_fuse = aircraft['b_{fuse}']
        b_wing = aircraft['b_{wing}']
        S_ref = aircraft['S_{ref}']
    
        prop_perf = self.prop = aircraft.propulsion_system.performance(flight_state, h_fuel, D_p_fuse, D_p_wing, BLI) # this is the important thing that makes it easy to swap in and out different propulsion systems...
        # Aircraft() expects an input propulsion_system that has a method
        # performance(self, flight_state) which has a variable with the string
        # name 'F_{\\rm net}' which is the total net propulsion system thrust
        

        constraints = [
            prop_perf,
            D_p_fuse == .5*flight_state['rho_\\infty']*flight_state['V_\infty']**2*S_ref*flight_state['C_{d,fuse}'],
            D_p_wing == .5*flight_state['rho_\\infty']*flight_state['V_\infty']**2*S_ref*flight_state['C_{d,wing}'],
            ]
        
        # Boundary-layer ingestion (BLI) fraction
        if BLI == 'on':
            constraints += [
            prop_perf['f_{fuse,BLI}'] == aircraft['d_{fan}']/b_fuse*prop_perf['\phi_{fuse}'], # BLI ingestion fraction per propulsor
            prop_perf['f_{wing,BLI}'] == aircraft['d_{fan}']/b_wing*prop_perf['\phi_{wing}'], # BLI ingestion fraction per propulsor
           # D == L / LDR,
           ]
        else:
            constraints += [
            prop_perf['f_{fuse,BLI}'] == 1e-32,
            prop_perf['f_{wing,BLI}'] == 1e-32,
         ]        

        return constraints
