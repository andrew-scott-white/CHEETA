from gpkit import Model, Variable, units

class Generator(Model):
    """A component sizing model of the generator
    """
    def setup(self):
        m_G = Variable('m_{gen}', 'kg', 'mass of generator') 
        T_m = Variable('(T/m)_{gen}', 'N*m/kg', 'generator specific torque')
        T_max = Variable('T_{max,gen}', 'N*m', 'generator maximum torque')
        
        constraints = [
            m_G >= T_max/T_m,
            ]
     
        return constraints
    
    
class GeneratorPerformance(Model):
    """A component performance model of the generator
    """
    def setup(self, generator, flight_state, core_perf):
        T = Variable('T', 'N*m', 'generator input shaft torque')
        omega_G = Variable('omega_{gen}', 'rpm', 'rotational speed')
        deta_G = Variable('1-eta_{gen}', '-', 'generator inefficiency')
        P_in = Variable('P_{in}', 'MW', 'input shaft power to generator')
        P_out = Variable('P_{out}', 'MW', 'output electrical power')
        Q = Variable('Q', 'MW', 'lost power')       

        
        constraints = [
            T == P_in / omega_G,
            P_in >= P_out + Q,
            Q == deta_G * P_in,
            generator['T_{max,gen}'] >= T,
            ]
        
        return constraints
    