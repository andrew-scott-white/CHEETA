from gpkit import Model, Variable, units, SignomialsEnabled
from gas_generator import GasGenerator, GasGeneratorPerformance
from elec_generator import Generator, GeneratorPerformance
from motor import Motor, MotorPerformance
from ducted_fan import DuctedFan, DuctedFanPerformance

class Turboelectric(Model):
    """Turboelectric (gas generator core and electrically-driven motors/propulsors) 
    propulsion system, consisting of N_eng cores and N_p propulsors.
    """
    def setup(self, N_eng, N_p):
        m_te = Variable('m_{te}', 'kg', 'turboelectric weight')
        m_propsys = Variable('m_{propsys}', 'kg', 'propulsion system weight')
        N_eng = Variable('N_{eng}', N_eng, '-', 'number of engines')
        N_p = Variable('N_p', N_p, '-', 'number of propulsors')
        
        W_fan = Variable('W_{fan}', 'kg', 'total fan weight')
        W_core = Variable('W_{core}', 'kg', 'total core weight')
        W_gen = Variable('W_{gen}', 'kg', 'total generator weight')
        W_mot = Variable('W_{mot}', 'kg', 'total motor weight')        
                 
        fan = DuctedFan()
        core = GasGenerator(fan)
        generator = Generator()
        motor = Motor()
        
    
        constraints = [
            core,
            generator,
            motor,
            fan,
            
            # propulsion system weight
            m_propsys >= N_eng*(core['m_{gg}'] + generator['m_{gen}'] + N_p/N_eng*(motor['m_{motor}'] + fan['m_{fan}'])),
            
            # total component sizings
            W_fan == N_p*fan['m_{fan}'],
            W_core == N_eng*core['m_{gg}'],
            W_mot == N_p*motor['m_{motor}'],
            W_gen == N_eng*generator['m_{gen}'],
            ]
    
        return constraints
    
    def performance(self, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI):
        return TurboelectricPerformance(self, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI)


class TurboelectricPerformance(Model):
    def setup(self, turboelectric, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI):
        F_net = Variable('F_{net}', 'N', 'net propulsion system thrust')
        mdot = Variable('-\\dot{m}', 'kg/s',
            'total propulsion system fuel consumption')
        SFC = Variable('SFC', 'kg/hr/N',
            'thrust-specific fuel consumption')
        P_m = Variable('(P_m)_{sys}', 'W/kg', 'propulsion system specific power')
        
        # BLI variables
        f_fuse_BLI = Variable('f_{fuse,BLI}', '-', 'fuselage BLI fraction')
        phi_fuse = Variable('\phi_{fuse}', '-', 'fuselage dissipation fraction')
        
        
        f_wing_BLI = Variable('f_{wing,BLI}', '-', 'wing BLI fraction')
        phi_wing = Variable('\phi_{wing}', '-', 'wing dissipation fraction')
        
        
        f_wake = Variable('1 - f_{wake}', '-', 'wake dissipation fraction')
        
        fan_fuse_perf = self.fan_fuse_perf = DuctedFanPerformance(turboelectric, flight_state, f_fuse_BLI, f_wake, D_p_fuse, BLI)
        fan_wing_perf = self.fan_wing_perf = DuctedFanPerformance(turboelectric, flight_state, f_wing_BLI, f_wake, D_p_wing, BLI)
        core_perf = GasGeneratorPerformance(turboelectric, flight_state, h_fuel)
        gen_perf = GeneratorPerformance(turboelectric, flight_state, core_perf)
        motor_perf = MotorPerformance(turboelectric, flight_state)
        
        constraints = [
            core_perf,
            gen_perf,
            motor_perf,
            fan_fuse_perf,
            fan_wing_perf,
            
            # power propagation through propulsion system
            core_perf['P_{out}'] >= gen_perf['P_{in}'],
            turboelectric['N_{eng}'] * gen_perf['P_{out}'] >= turboelectric['N_p'] * motor_perf['P_{in}'],
           turboelectric['N_p']* motor_perf['P_{out}'] >= turboelectric['N_p']*fan_fuse_perf['P_{shaft}']/3 + turboelectric['N_p']*fan_wing_perf['P_{shaft}']*2/3,
                      
            # system performance
            mdot >= turboelectric['N_{eng}'] * core_perf['\\dot{m}_{f}'],
            
            # post-processing performance parameters
          #  P_m == turboelectric['N_p'] * fan_perf['P_K'] / turboelectric['m_{propsys}'],
          #  SFC == turboelectric['N_{eng}'] * core_perf['\\dot{m}_{f}'] / (turboelectric['N_p']*fan_perf['F_{fan}']),
            ]
        
        with SignomialsEnabled():
            constraints += [
                fan_fuse_perf['F_{fan}']/3 + fan_wing_perf['F_{fan}']*2/3 >= F_net / turboelectric['N_p'],
                ]
    
        return constraints