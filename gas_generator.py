from gpkit import Model, Variable, units
from flight_state import FlightState
from ducted_fan import DuctedFan, DuctedFanPerformance

class GasGenerator(Model):
    """A component sizing model of the gas generator.
    """
    def setup(self, fan):
        m_gg = Variable('m_{gg}', 'kg', 'gas generator mass')
        P_m = Variable('(P/m)_{gg}', 'W/kg', 'gas generator specific power')
        P_max = Variable('P_{max,gg}', 'MW', 'max gas generator shaft power')
                           
        constraints = [
            m_gg >= P_max / P_m, # (1)
        ]

        return constraints
    
    
class GasGeneratorPerformance(Model):
    """A component performance model of the gas generator.
    """
    def setup(self, gas_generator, flight_state, h_fuel):
        eta_th = Variable('\\eta_{th}', '-','gas generator thermal efficiency')
        P_out = Variable('P_{out}', 'MW', 'total output shaft power')
        h_fuel = Variable('HV_{gg}', h_fuel, 'MJ/kg', 'fuel heating value')
        mdot_fuel = Variable('\\dot{m}_{f}', 'kg/s','gas generator fuel flow')
    
        constraints = [
            P_out == mdot_fuel*h_fuel*eta_th,
            P_out <= gas_generator['P_{max,gg}'],
            ]

        return constraints
    
    
# to test component, comment out gg sizing constraint (1)

if __name__ == '__main__':
    flight_state_gg = FlightState()
    fan = DuctedFan()
    fan_perf = DuctedFanPerformance(fan, flight_state_gg, 0.37)
    gg = GasGenerator(fan)
    ggPerf = GasGeneratorPerformance(gg, flight_state_gg, 43)
    ggPerf.substitutions.update({
        'HV_{gg}': 43*units('MJ/kg'),
        'P_{out}': 7.225*units('MW'),
        '\\eta_{th}':   .45
        })
    
    model = Model(ggPerf['\\dot{m}_{f}'], [ggPerf])
    sol = model.solve()
    print(sol.table())