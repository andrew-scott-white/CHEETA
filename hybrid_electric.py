from gpkit import Model, Variable, units, SignomialsEnabled

from ducted_fan import DuctedFan
from motor import Motor

from fuel_cell import FuelCell
from heat_exchanger import HeatExchanger

from gas_generator import GasGenerator, GasGeneratorPerformance
from elec_generator import Generator, GeneratorPerformance

g = 9.81*units('m/s/s')
gamma = 1.4

class HybridElectric(Model):
    """Fully electric propulsion powered by fuel cells"""
    
    # assumes...
        # one HEX per fuel cell stack
        # no mechanically-driven propulsors
        
    def setup(self, N_p):
        m_propsys = Variable('m_{propsys}', 'kg', 'propulsion system weight')
        
         
        N_s = Variable('N_s', '-', 'number of fuel cell stacks')
        N_eng = Variable('N_{eng}', '-', 'number of engines')
        N_p = Variable('N_p', N_p, '-', 'number of electrically-driven propulsors')
        P_m = Variable('(P/m)_{fc}', 'W/kg', 'fuel cell stack specific power')
        
        #nu_c = Variable('//nu_{conv}', '-', 'conversion electrification factor')
        
        W_fan = Variable('W_{fan}', 'kg', 'total fan weight')
        W_mot = Variable('W_{mot}', 'kg', 'total motor weight')
        W_stack = Variable('W_{stack}', 'kg', 'total stack weight')
        W_comp = Variable('W_{comp}', 'kg', 'total stack compression system weight')
        W_core = Variable('W_{core}', 'kg', 'total core weight')
        W_gen = Variable('W_{gen}', 'kg', 'total generator weight')
        W_HEX = Variable('W_{HEX}', 'kg', 'total HEX weight')
        A_HEX = Variable('A_{HEX}', 'm**2', 'total HEX area')
        
        fan = self.fan = DuctedFan()
        motor = self.motor = Motor()
        stack = self.stack = FuelCell()
        HEX = self.hx = HeatExchanger()
        core = GasGenerator(fan)
        generator = Generator()
        
        constraints = [
            fan,
            motor,
            stack,
            HEX,
            core,
            generator,
            
            P_m == stack['(P/m)_{fc}'], # stack specific power
            
            # propulsion system weight
            m_propsys >= N_s*(stack['m_{fc}'] + HEX['m_{HEX}']) + N_eng*(core['m_{gg}'] + generator['m_{gen}']) + N_p*(motor['m_{motor}'] + fan['m_{fan}']),
            
           
            
            # total component weights
            W_fan == N_p*fan['m_{fan}'], # total fan weight
            W_mot == N_p*motor['m_{motor}'], # total motor weight
            W_stack == N_s*stack['m_{fc}'], # total stack weight
            W_comp == N_s*stack['m_{comp}'], # total stack compressor weight
            W_HEX == N_s*HEX['m_{HEX}'], # total HEX weight
            A_HEX == N_s*HEX['A_{face}'], # total HEX area
            W_core == N_eng*core['m_{gg}'], # total core weight
            W_gen == N_eng*generator['m_{gen}'], # total generator weight
            ]

       # with SignomialsEnabled():
          #  constraints += [
                #N_s <= nu_c*(N_s + N_eng),
                #N_s + N_eng >= 2, # redundancy requirement in case of power system failure
               # ]

        return constraints
    
    def performance(self, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI):
        return HybridElectricPerformance(self, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI)

    
class HybridElectricPerformance(Model):
    def setup(self, hybrid, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI):
        F_net = Variable('F_{net}', 'N', 'net propulsion system thrust')
        mdot = Variable('-\\dot{m}', 'kg/s',
            'total propulsion system fuel consumption')
        #SFC = Variable('SFC', 'kg/hr/N',
        #    'thrust-specific fuel consumption')
       # P_m = Variable('(P_m)_{sys}', 'W/kg', 'propulsion system specific power')
        # BLI variables
        f_fuse_BLI = Variable('f_{fuse,BLI}', '-', 'fuselage BLI fraction')
        phi_fuse = Variable('\phi_{fuse}', '-', 'fuselage dissipation fraction')
        
        
        f_wing_BLI = Variable('f_{wing,BLI}', '-', 'wing BLI fraction')
        phi_wing = Variable('\phi_{wing}', '-', 'wing dissipation fraction')
        
        
        f_wake = Variable('1 - f_{wake}', '-', 'wake dissipation fraction')
        
        
        fan_fuse_perf = self.fan_fuse_perf = hybrid.fan.performance(flight_state, f_fuse_BLI, f_wake, D_p_fuse, BLI)
        fan_wing_perf = self.fan_wing_perf = hybrid.fan.performance(flight_state, f_wing_BLI, f_wake, D_p_wing, BLI)
        motor_perf = self.motor_perf = hybrid.motor.performance(flight_state)
        stack_perf = self.stack_perf = hybrid.stack.performance(flight_state, h_fuel)
        HEX_perf = self.HEX_perf = hybrid.hx.performance(flight_state)
        core_perf = self.core_perf = GasGeneratorPerformance(hybrid, flight_state, h_fuel)
        gen_perf = self.gen_perf = GeneratorPerformance(hybrid, flight_state, core_perf)
        
        constraints = [
            fan_fuse_perf,
            fan_wing_perf,
            motor_perf,
            stack_perf,
            HEX_perf,
            core_perf,
            gen_perf,
            
            stack_perf['HV_{fc}'] == core_perf['HV_{gg}'],
            
            # system fuel consumption
            mdot >= hybrid['N_s']*stack_perf['\\dot{m}_{f}'] + hybrid['N_{eng}']*core_perf['\\dot{m}_{f}'],
            
            core_perf['P_{out}'] >= gen_perf['P_{in}'],
            hybrid['N_p']*motor_perf['P_{out}'] >= hybrid['N_p']*fan_fuse_perf['P_{shaft}']/3 + hybrid['N_p']*fan_wing_perf['P_{shaft}']*2/3,
            ]
        
        with SignomialsEnabled():
            constraints += [
            # power propagation through propulsion system
            hybrid['N_s']*stack_perf['P_{out}'] + hybrid['N_{eng}']*gen_perf['P_{out}'] >= hybrid['N_p']*motor_perf['P_{in}'],
            
            hybrid['N_p']*fan_fuse_perf['F_{fan}']/3 + hybrid['N_p']*fan_wing_perf['F_{fan}']*2/3 + hybrid['N_s']*HEX_perf['F_{HEX}'] >= \
            hybrid['N_s']*HEX_perf['D_{core}'] + hybrid['N_s']*HEX_perf['D_{Nacelle}'] + F_net, # thrust distribution            
                     ]
            constraints += [
            # thermal system
            HEX_perf['T_{h,in}'] == flight_state['T_{t\\2}'],
            HEX_perf['Q_{HEX}'] == stack_perf['Q_{fc}'],
            
            # post-processing: performance metrics
            #SFC == mdot / (hybrid['N_p']*fan_perf['F_{fan}']),
            #P_m == hybrid['N_p']*fan_perf['P_K'] / hybrid['m_{propsys}'],
                 ]
        
        return constraints
  