from gpkit import Model, Variable, units, Vectorize

R1 = 287*units('J/kg/K')
gamma = 1.4

class FlightState(Model):
    """Set up free stream conditions used in component performance models.
    """
    def setup(self, C_D_fuse, C_D_wing):
    
        Pt0 = Variable('P_{t\\infty}', 'kPa',
            'free stream stagnation pressure')
        P0 = Variable('P_\\infty', 'kPa', 'free stream static pressure')
        Tt0 = Variable('T_{t\\infty}', 'K',
            'free stream stagnation temperature')
        T0 = Variable('T_\\infty', 'K', 'free stream static temperature')
        V0 = Variable('V_\infty', 'm/s', 'free stream velocity')
        M0 = Variable('M_\infty', '-', 'free stream Mach number')
        Cp_0 = Variable('C_p_\infty', 'J/kg/K', 'free stream air specific heat')
        rho0 = Variable('rho_\infty', 'kg/m**3', 'free stream air density')
        Pt2 = Variable('P_{t\\2}', 'kPa',
            'stagnation pressure')
        Tt2 = Variable('T_{t\\2}', 'K',
            'stagnation temperature')
        mu0 = Variable('\mu_\infty', 'Pa*s', 'free stream air viscosity')

        Cd_fuse = Variable('C_{d,fuse}', C_D_fuse, '-', 'fuselage profile drag coefficient')
        Cd_wing = Variable('C_{d,wing}', C_D_wing, '-', 'wing profile drag coefficient')
        dt = Variable('dt', 's', 'time')
        
        
        V0 = M0*(gamma*R1*Tt0)**0.5
        
        constraints = [
            ]

        return constraints
    

    
