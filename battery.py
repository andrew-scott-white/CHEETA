from gpkit import Model, Variable, units
from flight_state import FlightState
from path import C_D_fuse, C_D_wing

class Battery(Model):
    """A component sizing model of the battery.
    """
    def setup(self):
        m_b = Variable('m_{bat}', 'kg', 'battery mass')
        P_max = Variable('P_{max,bat}', 'kW', 'battery max power')
        P_m = Variable('(P/m)_{bat}', 'kW/kg', 'battery specific power')
        E_max = Variable('E_{bat,max}', 'kJ', 'maximum energy required of battery')
        E_m = Variable('(E/m)_{\\rm bat}', 'kJ/kg', 'battery specific energy')
        
        constraints = [
            m_b >= P_max/P_m,
            m_b == E_max/E_m,
            ]

        return constraints
    
    def performance(self, flight_state):
        return BatteryPerformance(self, flight_state)

class BatteryPerformance(Model):
    """A component performance model of the battery.
    """
    def setup(self, battery, flight_state):
        i = Variable('i', 'amp', 'battery current')
        i_max = Variable('i_{max}', 'amp', 'max battery current')
        e = Variable('e', 'V', 'battery electromotive force')
        P_out = Variable('P_{out}', 'kW', 'battery output power')
        Q_b = Variable('Q_b', 'kW', 'heat loss of battery')
        ohm = Variable('ohm_b', 'ohm', 'battery resistance')
        eta_b = Variable('eta_{b}', '-', 'battery efficiency')
        #dt = Variable('\Delta t', 's', 'time battery spends at a certain power setting')
        E = Variable('E_{bat}', 'kJ', 'energy required of battery')
    
        constraints = [
            battery,
            e == P_out/i,
            battery['P_{max,bat}'] >= P_out + Q_b,
            eta_b == P_out/battery['P_{max,bat}'], 
            Q_b == i**2*ohm,
            ohm == e/2/i_max,
            battery['P_{max,bat}'] >= i_max*e,
            i_max >= i,
            E == P_out*flight_state['dt'],
            battery['E_{bat,max}'] >= E,
            #dt >= .0001*units('s')
            ]

        return constraints
    
if __name__ == '__main__':
    battery = Battery()
    flight_state = FlightState(C_D_fuse, C_D_wing)
    batPerf = BatteryPerformance(battery, flight_state)
    batPerf.substitutions.update({
        '(P/m)_{bat}': 800*units('W/kg'), # conservative 2035 future estimate from LEARN report pg. 14
        '(E/m)_{\\rm bat}': 4*10**2*units('W*hr/kg'), # estimate of Li ion batteries based on [Viswanathan, Nature, 2022]
        'ohm_b': 10.2*units('ohm'), # estimate based on Tesla Model S Battery
        #'P_{max,bat}': 325*units('W'), # estimate based on Tesla Model S Battery
        'P_{out}': 100*units('kW'),       
        #'eta_{b}':.90,
        #'E_{bat,max}': 100000*units('kJ'),
        #'m_{bat}': 80*units('lb'),
        'dt': 1800*units('s'),
        
        })
    
    model = Model(batPerf['m_{bat}'], [batPerf])
    sol = model.solve()
    print(sol.table())
    