from gpkit import Model, Variable, units, SignomialsEnabled
from ducted_fan import DuctedFan
from motor import Motor
from fuel_cell import FuelCell
from heat_exchanger import HeatExchanger

g = 9.81*units('m/s/s')
gamma = 1.4

class FullyElectric(Model):
    
      # assumes...
        # one HEX per fuel cell stack
        
    """Fully electric propulsion powered by fuel cells"""
    def setup(self, N_p):
        N_s = Variable('N_s', '-', 'number of fuel cell stacks')
        m_propsys = Variable('m_{propsys}', 'kg', 'propulsion system weight')
        N_p = Variable('N_p', N_p, '-', 'number of propulsors')
        P_m = Variable('(P/m)_{fc}', 'W/kg', 'fuel cell stack specific power')
        
        W_fan = Variable('W_{fan}', 'kg', 'total fan weight')
        W_mot = Variable('W_{mot}', 'kg', 'total motor weight')
        W_stack = Variable('W_{stack}', 'kg', 'total stack weight')
        W_comp = Variable('W_{comp}', 'kg', 'total stack compression system weight')
        W_HEX = Variable('W_{HEX}', 'kg', 'total HEX weight')
        A_HEX = Variable('A_{HEX}', 'm**2', 'total HEX area')
        A_fan = Variable('A_{fan}', 'm^2', 'total fan area')
        alpha = Variable('\\alpha_{HEX}', '-', 'HEX-fan area ratio')
        
        #span = Variable('l_{span}', 'ft', 'aircraft wing span')
        
        fan = self.fan = DuctedFan()
        motor = self.motor = Motor()
        stack = self.stack = FuelCell()
        HEX = self.hx = HeatExchanger()
        
        constraints = [
            fan,
            motor,
            stack,
            HEX,
            
            P_m == stack['(P/m)_{fc}'], # stack specific power
            
            m_propsys >= N_s*(stack['m_{fc}'] + HEX['m_{HEX}']) + N_p*(motor['m_{motor}'] + fan['m_{fan}']), # propulsion system weight
            N_s >= 2, # stack redundancy in case of power system failure
           
            # total component weights
            W_fan == N_p*fan['m_{fan}'],
            W_mot == N_p*motor['m_{motor}'],
            W_stack == N_s*stack['m_{fc}'],
            
            W_comp == N_s*stack['m_{comp}'],
            W_HEX == N_s*HEX['m_{HEX}'],
            A_HEX == N_s*HEX['A_{face}'],
            A_fan == N_p*3.1416/4*fan['d_{fan}']**2,
            alpha == A_HEX/A_fan,
            ]

        return constraints
    

    
    def performance(self, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI):
        return FullyElectricPerformance(self, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI)

    
class FullyElectricPerformance(Model):
    def setup(self, fully_electric, flight_state, h_fuel, D_p_fuse, D_p_wing, BLI):
        F_net = Variable('F_{net}', 'kN', 'net propulsion system thrust')
        mdot = Variable('-\\dot{m}', 'kg/s',
            'total propulsion system fuel consumption during cruise')
        SFC = Variable('SFC', 'kg/hr/N',
            'thrust-specific fuel consumption')
        P_m = Variable('(P_m)_{sys}', 'W/kg', 'propulsion system specific power')
        HV = Variable('HV_{gg}', 'MJ/kg', 'fuel heating value of fuel') # satisfies the PFEI constraint in mission
        
        # BLI variables
        f_fuse_BLI = Variable('f_{fuse,BLI}', '-', 'fuselage BLI fraction')
        phi_fuse = Variable('\phi_{fuse}', '-', 'fuselage dissipation fraction')
        
        
        f_wing_BLI = Variable('f_{wing,BLI}', '-', 'wing BLI fraction')
        phi_wing = Variable('\phi_{wing}', '-', 'wing dissipation fraction')
        
        
        f_wake = Variable('1 - f_{wake}', '-', 'wake dissipation fraction')
                
        fan_fuse_perf = self.fan_fuse_perf = fully_electric.fan.performance(flight_state, f_fuse_BLI, f_wake, D_p_fuse, BLI)
        fan_wing_perf = self.fan_wing_perf = fully_electric.fan.performance(flight_state, f_wing_BLI, f_wake, D_p_wing, BLI)
        motor_perf = fully_electric.motor.performance(flight_state)
        stack_perf = fully_electric.stack.performance(flight_state, h_fuel)
        HEX_perf = fully_electric.hx.performance(flight_state)
        
        constraints = [
            fan_fuse_perf,
            fan_wing_perf,
            motor_perf,
            stack_perf,
            HEX_perf,
            
            HV == stack_perf['HV_{fc}'],
        
            # power propagation through propulsion system
            fully_electric['N_s']*stack_perf['P_{out}'] >= fully_electric['N_p']*motor_perf['P_{in}'],
            fully_electric['N_p']*motor_perf['P_{out}'] >= fully_electric['N_p']*fan_fuse_perf['P_{shaft}']/3 + fully_electric['N_p']*fan_wing_perf['P_{shaft}']*2/3,
                            
            # Assume radiator temperature is equal to fuel cell operating temperature
            HEX_perf['T_{h,in}'] == flight_state['T_{t\\2}'],
            ]
        
        with SignomialsEnabled():
            constraints += [
            # system performance: force balance
            fully_electric['N_p']*fan_fuse_perf['F_{fan}']/3 + fully_electric['N_p']*fan_wing_perf['F_{fan}']*2/3 + fully_electric['N_s']*HEX_perf['F_{HEX}'] >= \
            fully_electric['N_s']*HEX_perf['D_{core}'] + fully_electric['N_s']*HEX_perf['D_{Nacelle}'] + F_net, # thrust distribution
                      ]
            
            constraints += [
            mdot == fully_electric['N_s']*stack_perf['\\dot{m}_{f}'], # total fuel consumption
            
            # thermal system (heat sent to radiator is equal to heat rejected by stack)
            HEX_perf['Q_{HEX}'] == stack_perf['Q_{fc}'],
            
            
            
            # post-processing performance parameters
            #SFC == fully_electric['N_s'] * stack_perf['\\dot{m}_{f}'] / (fully_electric['N_p']*fan_perf['F_{fan}']),
            #P_m == fully_electric['N_p']*fan_perf['P_K'] / fully_electric['m_{propsys}'],
            #fully_electric['N_s']*stack_perf['P_{max,fc}'] >= fully_electric['N_p']*fan_perf['P_{K,}'],
        ]
        
        return constraints
    

