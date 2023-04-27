from gpkit import Model, Variable, units
from gpkit import SignomialsEnabled, SignomialEquality
from gpkit import Vectorize
from numpy import pi



mu_0 = 4*pi*1e-07 *units('H/m') # Vacuum permeability

class MotorSizing(Model):
    def setup(self):

        # radial build
        D_ag = Variable('D_{\\rm ag}', 'm', 'airgap diameter')
        D_out = Variable('D_{\\rm out}', 'm', 'machine outer diameter')
        d_c = Variable('d_c', 'm',
            'difference between outer dia. and airgap dia.')

        constraints = [
            D_out >= D_ag + d_c,
        ]

        # circumferential build
        N_a = Variable('N_a', '-', 'number of armature conductors per phase')
        N_P = Variable('N_P', '-', 'no. of poles')
        N_PP = Variable('N_{\\rm PP}', '-',
            'number of parallel paths in armature coils')
        N_phase = Variable('N_{\\rm phase}', '-', 'number of power phases')
        N_slot = Variable('N_{\\rm slot}', '-', 'number of slots')
        N_slotpole = Variable('N_{\\rm slot-pole}', '-',
            'No. of armature winding slots per pole')
        N_turns = Variable('N_{\\rm turns}', '-', 'number of turns per slot')

        constraints += [
            N_a == N_slot * N_turns / N_phase,
            N_slot == N_P * N_slotpole,
        ]
        
        # armature length
        A_s  = Variable('A_s', 'm*m', 'armature slot area')
        d_slot = Variable('d_{\\rm slot}', 'm',
            'radial depth of armature winding slots')
        d_iw = Variable('d_{\\rm iw}', 'm', 'insulation width')
        L = Variable('L', 'm', 'machine active length')
        O_D  = Variable('O_D', 'm', 'conductor diameter')

        with SignomialsEnabled():
            constraints += [
                SignomialEquality(A_s,
                    pi *((D_ag/2 + d_slot)**2 - (D_ag/2)**2
                        - N_slotpole * d_iw * d_slot * N_P)/(N_slotpole * N_P)),
            ]

        L_arm = Variable('L_{\\rm arm}', 'm', 'total armature length')
        L_e = Variable('L_e', 'm', 'end winding length')

        constraints += [
            L_arm >= (L + L_e) * (A_s * N_slotpole * N_P)/pi/(O_D/2)**2,
        ]

        # line and armature winding current
        I_aR = Variable('I_{a,R}', 'A', 'supply line rated current')
        J_R = Variable('J_R', 'A/m^2',
            'maximum operating current density in one conductor')
        lam_s = Variable('\\lambda_s', '-',
            'Fill factor of superconducting wire')

        constraints += [
            I_aR == A_s * J_R * lam_s / N_turns,
        ]

        # critical armature conductor current
        I_c = Variable('I_c', 'A', 'critical armature current')
        J_c = Variable('J_c', 'A/m^2',
            'critical current density in armature windings')
        lam_c = Variable('\\lambda_c', '-',
            'Fill factor of armature conductors')

        constraints += [
            I_c == lam_c * J_c * pi * (O_D/2)**2,
        ]

        # other constant parameters used in performance model
        B_g = Variable('B_g', 'T', 'airgap peak magnetic flux density')
        d_f = Variable('d_f', 'm', 'conductor filament diameter')
        k = Variable('k', '-', 'eddy current power dissipation constant')
        L_p = Variable('L_p', 'm', 'twist length of conducting wires')
        n = Variable('n', '-', 'coupling loss constant')
        rho_eff = Variable('\\rho_{\\rm eff}', 'ohm*m',
            'effective resistivity matrix of conductor material')

        return constraints


class MotorPerformance(Model):
    def setup(self, motor):

        # torque, speed, mechanical power, and physical dimensions
        A = Variable('A', 'A/m', 'electric loading')
        I_a = Variable('I_a', 'A', 'supply line current')
        I_W = Variable('I_W', 'A', 'winding bundle current')
        P_I = Variable('P_I', 'W', 'ideal machine power')
        P_0 = Variable ('P_0', 'W', 'shaft power')
        T_e = Variable('T_e', 'N*m', 'torque of electric origin')
        V = Variable('V', 'V', 'supply voltage')
        omega_m = Variable('\\omega_{\\rm mot}', 'rpm', 'motor angular speed')

        B_g = motor['B_g']
        D_ag = motor['D_{\\rm ag}']
        L = motor['L']
        N_a = motor['N_a']
        N_P = motor['N_P']

        constraints = [
            P_I == T_e * omega_m,
            P_I == V * I_a,
            T_e == 6 * pi * (D_ag**2 * L / 4) * B_g * A,
            A == N_a * I_W / pi / D_ag,
            I_a == motor['N_{\\rm PP}'] * I_W,
        ]

        # conductor volume
        L_arm = motor['L_{\\rm arm}']
        O_D = motor['O_D']
        vol = L_arm * pi * (O_D/2)**2

        # electrical frequency
        f = Variable('f', 'Hz', 'electrical frequency')
        constraints += [
            omega_m == 2 * pi * f / (N_P/2),
        ]

        # hysteresis loss
        J_c = motor['J_c']
        P_h = Variable('P_h', 'W', 'Hysteresis loss')

        lam_c = motor['\\lambda_c']
        d_f = motor['d_f']

        constraints += [
            omega_m == 2 * pi * f / (N_P/2),
            P_h >= (4/3) * B_g * J_c * d_f * f * lam_c * vol,
        ]

        # eddy current loss
        P_e = Variable('P_e', 'W', 'Eddy current loss')

        k = motor['k']
        rho_eff = motor['\\rho_{\\rm eff}']

        constraints += [
            P_e >= (pi**2/(k * rho_eff)) * (B_g * O_D * f)**2 * vol, 
        ]

        # coupling loss
        P_c = Variable('P_c', 'W', 'Coupling loss')

        n = motor['n']

        L_p = motor['L_p']

        constraints += [
            P_c >= (1/(n * rho_eff)) * (B_g * L_p * f)**2 * vol,
        ]

        # transport loss
        I_con = Variable('I_{\\rm con}', 'A', 'armature conductor current')
        J_O = Variable('J_O', 'A/m^2',
            'operating current density in one conductor')
        P_t = Variable('P_t', 'W', 'Transport current loss')

        A_s = motor['A_s']
        I_c = motor['I_c']
        J_R = motor['J_R']
        N_turns = motor['N_{\\rm turns}']
        lam_s = motor['\\lambda_s']

        constraints += [
            I_con == J_O * pi * (O_D/2)**2,
            J_O == I_W * N_turns / A_s / lam_s,
            J_O <= J_R,
        ]

        with SignomialsEnabled():
            constraints += [
                #SignomialEquality(
                    P_t >= 
                    vol * mu_0/pi * f * I_c**2 * ((1- I_con/I_c) 
                        * (-I_con/I_c - (I_con/I_c)**2/2 - (I_con/I_c)**3/3)
                        + I_con/I_c - 0.5*(I_con/I_c)**2) / (pi* (O_D/2)**2),
                #),
            ]

        # total power
        eta_mot = Variable('\\eta_{\\rm mot}', '-', 'motor efficiency')

        constraints += [
            P_I >= P_0 + P_h + P_e + P_c + P_t,
            eta_mot == P_0 / P_I,
        ]

        return constraints

if __name__ == '__main__':
    from x0 import x0_motor

    # multi-point CHEETA motor sizing
    motor = MotorSizing()
    motor.substitutions.update({
        # radial build
        'd_c': 0.19 * units('m'),
        # armature length calculation
        'd_{\\rm iw}': 0.001/100 * units('m'),
        'd_{\\rm slot}': 0.00252 * units('m'),
        #'L': 0.8569 * units('m'),
        'L_e': 0.094 * units('m'),
        'N_{\\rm PP}': 1,
        'N_{\\rm phase}': 3,
        'N_{\\rm slot-pole}': 6,
        'N_{\\rm turns}': 1,
        'O_D': 0.032/100 * units('m'),
        # line and armature winding current
        'J_R': 200.475e6 * units('A/m^2'),
        '\\lambda_s': 0.5,
        # other conductor parameters (used in loss models)
        'B_g': 0.35*units('T'),
        'd_f': 0.001/100 * units('m'),
        'J_c': 4.3831e09 * units('A/m^2'),
        'k': 4,
        'L_p': 0.005 * units('m'),
        'n': 2,
        '\\lambda_c': 0.15,
        '\\rho_{\\rm eff}': 1.25e-07 * units('ohm*m'),
    })

    D_ref = Variable('D_{\\rm ref}', 0.304, 'm', 'reference airgap diameter')
    N_ref = Variable('N_{\\rm ref}', 8, '-', 'reference number of poles')

    with Vectorize(5):
        motor_perf = MotorPerformance(motor)

    motor_perf.substitutions.update({
        motor_perf['P_0'][0]: 2.76*units('MW'), 
        motor_perf['\\omega_{\\rm mot}'][0]: 5430*units('rpm'),
        #
        motor_perf['P_0'][1]: 1.42*units('MW'), 
        motor_perf['\\omega_{\\rm mot}'][1]: 5530*units('rpm'),
        #
        motor_perf['P_0'][2]: 1.00*units('MW'), 
        motor_perf['\\omega_{\\rm mot}'][2]: 5320*units('rpm'),
        #
        motor_perf['P_0'][3]: 1.01*units('MW'), 
        motor_perf['\\omega_{\\rm mot}'][3]: 5480*units('rpm'),
        #
        motor_perf['P_0'][4]: 0.907*units('MW'), 
        motor_perf['\\omega_{\\rm mot}'][4]: 5260*units('rpm'),
    }) 

    model = Model(
        sum(motor_perf['P_I']) / 5
            * motor['D_{\\rm out}']**2*motor['L'],
        [
            motor,
            motor['D_{\\rm out}'] <= 5.606 * 2 * units('in'),
            motor_perf,
            motor['N_P'] == N_ref * (motor['D_{\\rm ag}'] / D_ref),
        ],
    )

    sol = model.localsolve(
        solver='mosek_conif',
        verbosity=2,
        x0 = x0_motor,
    )
    print(sol.table())

    ## single point 2.55 MW motor sizing
    #motor = MotorSizing()
    #motor.substitutions.update({
    #    # radial build
    #    'd_c': 0.19 * units('m'),
    #    # armature length calculation
    #    'd_{\\rm iw}': 0.001/100 * units('m'),
    #    'd_{\\rm slot}': 0.00252 * units('m'),
    #    'L': 0.8569 * units('m'),
    #    'L_e': 0.094 * units('m'),
    #    'N_{\\rm PP}': 1,
    #    'N_{\\rm phase}': 3,
    #    'N_{\\rm slot-pole}': 6,
    #    'N_{\\rm turns}': 1,
    #    'O_D': 0.032/100 * units('m'),
    #    # line and armature winding current
    #    'J_R': 200.475e6 * units('A/m^2'),
    #    '\\lambda_s': 0.5,
    #    # other conductor parameters (used in loss models)
    #    'B_g': 0.35*units('T'),
    #    'd_f': 0.001/100 * units('m'),
    #    'J_c': 4.3831e09 * units('A/m^2'),
    #    'k': 4,
    #    'L_p': 0.005 * units('m'),
    #    'n': 2,
    #    '\\lambda_c': 0.15,
    #    '\\rho_{\\rm eff}': 1.25e-07 * units('ohm*m'),
    #})

    #D_ref = Variable('D_{\\rm ref}', 0.304, 'm', 'reference airgap diameter')
    #N_ref = Variable('N_{\\rm ref}', 8, '-', 'reference number of poles')

    #motor_perf = MotorPerformance(motor)
    #motor_perf.substitutions.update({
    #    'P_0': 2.5452 * units('MW'),
    #    '\\omega_{\\rm mot}': 4500 * units('rpm'),
    #    'A': 40964 * units('A/m'),
    #})

    #model = Model(
    #    motor['D_{\\rm out}']**2*motor['L'] 
    #        * motor_perf['P_I'],
    #    [
    #        motor,
    #        motor_perf,
    #        motor['N_P'] == N_ref * (motor['D_{\\rm ag}'] / D_ref),
    #    ],
    #)

    #sol = model.localsolve(solver='mosek_conif',verbosity=2)
    #print(sol.table())