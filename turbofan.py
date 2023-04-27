from gpkit import Model, Variable, units
from gas_generator import GasGenerator, GasGeneratorPerformance
from ducted_fan import DuctedFan, DuctedFanPerformance

class Turbofan(Model):
    """Turbofan (gas generator core plus mechanically-coupled propulsor)
    propulsion system, consisting of N_eng turbofans.
    """
    def setup(self, N_eng):
        m_propsys = Variable('m_{propsys}', 'kg',
            'total propulsion system mass')
        N_eng = Variable('N_{eng}', N_eng, '-', 'number of engines')
        N_p = Variable('N_p', '-', 'number of propulsors')

        W_fan = Variable('W_{fan}', 'kg', 'total fan weight')
        W_core = Variable('W_{core}', 'kg', 'total core weight')
        
        # build from constituent components
        fan = DuctedFan()
        core = GasGenerator(fan)

        constraints = [
            core,
            fan,
            
            # propulsion system mass
            m_propsys >= N_eng*core['m_{gg}'] + N_p*fan['m_{fan}'], # one turbofan
            N_eng == N_p, # fan and core are mechanically connected
            
            # total component sizings
            W_fan == N_p*fan['m_{fan}'],
            W_core == N_eng*core['m_{gg}'],
        ]

        return constraints

    def performance(self, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI):
        return TurbofanPerformance(self, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI)


class TurbofanPerformance(Model):
    def setup(self, turbofan, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI):
        F_net = Variable('F_{net}', 'N', 'net propulsion system thrust')
        #F_net_toc = Variable('F_{net,toc}', 'N', 'net propulsion system thrust at top-of-climb')
        #F_net_to = Variable('F_{net,to}', 'N', 'net propulsion system thrust at take-off')
        mdot = Variable('-\\dot{m}', 'kg/s',
            'total propulsion system fuel consumption')
        SFC = Variable('{SFC}', 'kg/hr/N',
            'thrust-specific fuel consumption')
        P_m = Variable('(P_m)_{sys}', 'W/kg', 'propulsion system specific power')
        
        # BLI variables
        f_fuse_BLI = Variable('f_{fuse,BLI}', '-', 'fuselage BLI fraction')
        phi_fuse = Variable('\phi_{fuse}', '-', 'fuselage dissipation fraction')
        
        
        f_wing_BLI = Variable('f_{wing,BLI}', '-', 'wing BLI fraction')
        phi_wing = Variable('\phi_{wing}', '-', 'wing dissipation fraction')
        
        
        f_wake = Variable('1 - f_{wake}', '-', 'wake dissipation fraction')

        # fuselage only fans
        fan_perf = self.fan_perf = DuctedFanPerformance(turbofan, flight_state, f_fuse_BLI, f_wake, D_p_fuse, BLI)
        core_perf = GasGeneratorPerformance(turbofan, flight_state, h_fuel)
        
        

        constraints = [
            core_perf,
            fan_perf,
            
            # power propagation through propulsion system
            core_perf['P_{out}'] >= fan_perf['P_{shaft}'], # cruise
             
            # system performance
            fan_perf['F_{fan}'] >= F_net / turbofan['N_{eng}'],
            mdot >= turbofan['N_{eng}'] * core_perf['\\dot{m}_{f}'],
            
            # post-processing performance parameters
            P_m == turbofan['N_p'] * fan_perf['P_K'] / turbofan['m_{propsys}'],
            SFC == core_perf['\\dot{m}_{f}'] / fan_perf['F_{fan}'],
                   ]

        return constraints
    