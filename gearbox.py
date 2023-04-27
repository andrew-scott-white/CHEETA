from gpkit import Model, Variable, units

class Gearbox(Model):
    """A component model for a gearbox"""
    def setup(self, P_d, N_r, omega_eng, omega_rot):
        w_gbrs = Variable('w_{\\rm gbrs}', 'lb', 'gearbox + rotating shaft weight')
        P_d = Variable('P_{\\rm d}', 'hp', 'drive system power limit')
        f_q = Variable('f_{\\rm q}', '-', 'second rotor torque limit (fraction of total drive system torque limit)')
        N_r = Variable('N_{\\rm gb}', '-', 'number of rotors')
        omega_eng = Variable('omega_{\\rm eng}', 'rpm', 'engine output speed')
        omega_rot = Variable('omega_{\\rm rotor}', 'rpm', 'main rotor speed')
        f_rs = Variable('f_{\\rm rs}', 0.13, '-', 'rotor shaft weight (fraction of gearbox and rotor shaft)')
        f_gb = Variable('f_{\\rm gb}', '-', 'gearbox shaft weight (fraction of gearbox and rotor shaft)')
        W_gb = Variable('W_{\\rm gb}', 'lb', 'weight of gearbox')
        W_rs = Variable('W_{\\rm rs}', 'lb', 'weight of rotating shaft')
        
        constraints = [
            w_gbrs == 95.7634*N_r**.38553*P_d**.78137*omega_eng**.09899/omega_rot**.80686,
            W_gb == f_gb*w_gbrs,
            W_rs == f_rs*w_gbrs,
            1 >= f_gb + f_rs,
            ]
        
        return constraints