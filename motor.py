from gpkit import Model, Variable

class Motor(Model):
    """A component sizing model of the motor
    """
    def setup(self):
        m_M = Variable('m_{motor}', 'kg', 'motor mass')
        T_m = Variable('(T/m)_{mot}', 'N*m/kg', 'motor specific torque')
        T_max = Variable('T_{max,mot}', 'N*m', 'motor maximum torque')
        
        constraints = [
            m_M >= T_max/T_m
            ]
        
        return constraints

    def performance(self, flight_state):
        return MotorPerformance(self, flight_state)
    
    
class MotorPerformance(Model):
    """A component performance model of the motor
    """
    def setup(self, motor, flight_state):
        T = Variable('T', 'N*m', 'motor output shaft torque')
        deta_M = Variable('1-eta_{motor}', '-', 'motor inefficiency')
        omega_M = Variable('omega_{motor}', 'rpm', 'rotational speed') 
        
        P_out = Variable('P_{out}', 'MW', 'output shaft power per motor')
        Q = Variable('Q_{mot}', 'kW', 'lost power')
        P_in = Variable('P_{in}', 'MW', 'input electrical power')
        
        constraints = [
            T == P_out/omega_M,
            P_in >= P_out + Q,
            Q == deta_M * P_in, # For fuel-cooled super conducting motors, this waste heat is dumped into the LH2. Otherwise, this cooling requirement must be book-kept.
            
            T <= motor['T_{max,mot}'], 
            ]

        return constraints