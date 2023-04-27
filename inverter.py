from gpkit import Model, Variable, units

class Inverter(Model):
    """A component sizing model of the battery.
    """
    def setup(self):
        m_inv = Variable('m_{inv}', 'lb', 'inverter mass')
        P_max = Variable('P_{max,inv}', 'kW', 'inverter max power')
        P_m = Variable('(P/m)_{inv}', 'kW/kg', 'inverter specific power')
        
        constraints = [
            m_inv >= P_max/P_m
            ]

        return constraints


class InverterPerformance(Model):
    """A component performance model of the battery.
    """
    def setup(self, battery):
        P_in = Variable('P_{in}', 'kW', 'inverter input power')
        P_out = Variable('P_{out}', 'kW', 'inverter output power')
        Q_inv = Variable('Q_{inv}', 'W', 'heat loss of inverter')
        eta_inv = Variable('eta_{inv}', '-', 'inverter efficiency')
    
        constraints = [
            
            P_in >= P_out + Q_inv,
            eta_inv == P_out/P_in,
            ]

        return constraints
    