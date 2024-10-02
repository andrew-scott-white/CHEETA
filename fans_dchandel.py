from gpkit import Model, Variable, Vectorize, units
from gpkit import SignomialsEnabled, SignomialEquality
from numpy import array, exp, tan, pi



# constants
g_mps2 = 9.80665
R_JpkgK = 287.058
gamma = 1.4

# atmospheric conditions not part of GPkit model
def isa(h_m):
    """Return static pressure and temperature at a given altitude for the
    International Standard Atmosphere.
    """
    if h_m < 11000:
        # Troposphere
        a = 19 + 273 # Base temperature (K)
        b = -0.0065 # Temperature lapse rate (K/m)
        p0 = 108900 # Base pressure (Pa)
        h0 = -610 # Base altitude (m)
        p = p0 * ((a+b*h_m)/(a+b*h0))**(-g_mps2/b/R_JpkgK)
    elif h_m < 20000:
        # Tropopause
        a = -56.5 + 273 # Base temperature
        b = 0 # Temperature lapse rate
        p0 = 22632 # Base pressure
        h0 = 11000 # Base altitude
        p = p0 * exp(-g_mps2*(h_m-h0)/R_JpkgK/a)
    else:
        # Stratosphere (up to 32000 m)
        a = -56.5 + 273 # Base temperature
        b = 0.001 # Temperature lapse rate
        p0 = 5474.9 # Base pressure
        h0 = 20000 # Base altitude
        p = p0 * ((a+b*h_m)/(a+b*h0))**(-g_mps2/b/R_JpkgK)

    T = a + b * (h_m-h0)
        
    return p, T

# optimization parameters
g = g_mps2 * units('m/s^2')
R = R_JpkgK * units('J/kg/K')


# GPkit models
class FlightState(Model):
    def setup(self):
        M = Variable('M', '-', 'Mach number')
        p = Variable('p', 'Pa', 'static pressure')
        pt = Variable('p_t', 'Pa', 'stagnation pressure')
        q = Variable('q', 'Pa', 'dynamic pressure')
        T = Variable('T', 'K', 'static temperature')
        Tt = Variable('T_t', 'K', 'stagnation temperature')
        V = Variable('V', 'm/s', 'velocity')
        rho = Variable('\\rho', 'kg/m^3', 'density')

        return [
            V == M * (gamma * R * T)**0.5,
            p == rho * R * T,
            q == 0.5 * rho * V**2,
        ]


class Fan(Model):
    def setup(self):
        A2 = Variable('A_2', 'm^2', 'fan area')
        A8 = Variable('A_8', 'm^2', 'nozzle area')
        cd_nace = Variable('c_{d,{\\rm nace}}', '-',
            'nacelle section drag coefficient')
        d_fan = Variable('d_{\\rm fan}', 'in', 'fan diameter')
        l_nace = Variable('\\ell_{\\rm nace}', 'm',
            'nacelle section chord length')
        l_d = Variable('\\ell/d', 1.78, '-', 'nacelle length-to-diameter ratio')
        m_fan = Variable('m_{\\rm fan}', 'kg', 'fan mass')
        rhub = Variable('r_{\\rm hub}', 'm', 'fan hub radius')
        rtip = Variable('r_{\\rm tip}', 'm', 'fan tip radius')
        rhub_rtip = Variable('r_{\\rm hub}/r_{\\rm tip}', '-',
            'fan hub-to-tip ratio')
        tan_beta = Variable('\\tan\\beta', '-',
            'tangent of rotor-relative tip exit flow angle')

        constraints = [
            rtip**2 >= rhub**2 + A2 / pi,
            rhub_rtip <= rhub / rtip,
            m_fan / 2600/units('lb') >= (2 * rtip / 61/units('in'))**2.4,
            l_nace == l_d * d_fan,
            d_fan == 2 * rtip,
        ]

        return constraints

    def performance(self, flight_state):
        return FanPerformance(self, flight_state)


class FanPerformance(Model):
    def setup(self, fan, flight_state):
        # Internal flow ("thrust and power")
        D2 = Variable('D_2', '-', 'fan face corrected flow per unit area')
        D8 = Variable('D_8', '-', 'nozzle throat corrected flow per unit area')
        M8 = Variable('M_8', '-', 'nozzle throat Mach number')
        mdot = Variable('\\dot{m}', 'kg/s', 'mass flow')
        p8 = Variable('p_8', 'Pa', 'nozzle throat static pressure')
        pt1 = Variable('p_{t1}', 'Pa', 'inlet stagnation pressure')
        pt2 = Variable('p_{t2}', 'Pa', 'fan face stagnation pressure')
        pt8 = Variable('p_{t8}', 'Pa', 'nozzle throat stagnation pressure')
        T8 = Variable('T_8', 'K', 'nozzle throat static temperature')
        Tt2 = Variable('T_{t2}', 'K', 'fan face stagnation temperature')
        Tt8 = Variable('T_{t8}', 'K', 'nozzle throat stagnation temperature')
        u8 = Variable('u_8', 'm/s', 'nozzle throat velocity')
        pi_inl = Variable('\\pi_{\\rm inl}', 0.995, '-',
            'inlet pressure recovery')
        pi_fan = Variable('\\pi_{\\rm fan}', '-', 'fan pressure ratio')
        rho8 = Variable('\\rho_8', 'kg/m^3', 'nozzle throat static density')
        #rhot8 = Variable('\\rho_{t8}', 'kg/m^3',
        #    'nozzle throat stagnation density')

        p0 = flight_state['p']

        eta_fan = 0.93

        constraints = [
            # Fan state change
            pt2 == pt1 * pi_inl,
            Tt2 == flight_state['T_t'],
            pt8 == pt2 * pi_fan,
            pt8 / pt2 == (Tt8 / Tt2) ** (gamma * eta_fan / (gamma - 1)),
            # Corrected flow per unit area
            D2 == mdot * (R * Tt2 / gamma)**0.5 / fan['A_2'] / pt2,
            D2 <= 0.5, # M2 <= 0.628
            D8 == mdot * (R * Tt8 / gamma)**0.5 / fan['A_8'] / pt8,
            (1 + (gamma - 1) / 2 * M8**2) \
                * D8 ** (2 * (gamma - 1) / (gamma + 1)) \
                <= M8 ** (2 * (gamma - 1) / (gamma + 1)),
            # Nozzle exit conditions
            Tt8 / T8 == (M8 / D8)**(2 * (gamma - 1) / (gamma + 1)),
            #pt8 == rhot8 * R * Tt8,
            p8 == rho8 * R * T8,
            p8 == p0,
            #rhot8 / rho8 == (pt8 / p8)**(1 / gamma),
            mdot == rho8 * u8 * fan['A_8'],
            M8 == u8 / (gamma * R * T8)**0.5,
            # Nozzle choking conditions, requires SP F_net...
            #M8 <= 1, # If this constraint is tight...
            #p8 >= p0, # ...this one should be loose.
        ]

        # Fan shaft power (main metric, but SP...)
        P_shaft = Variable('P_{\\rm shaft}', 'MW', 'fan shaft power')
        dht = Variable('\\Delta h_t', 'J/kg', 'fan stagnation enthalpy rise')

        with SignomialsEnabled():
            constraints += [
                dht >= gamma / (gamma - 1) * R * (Tt8 - Tt2),
            ]

        constraints += [
            P_shaft >= mdot * dht,
        ]

        # fan characteristic
        M2 = Variable('M_2', '-', 'fan face Mach number')
        Mtip = Variable('M_{\\rm tip}', '-', 'fan tip Mach number')
        V2 = Variable('V_2', 'm/s', 'fan face velocity')
        T2 = Variable('T_2', 'K', 'fan face static temperature')
        p2 = Variable('p_2', 'Pa', 'fan face static pressure')
        rho2 = Variable('\\rho_2', 'kg/m^3', 'fan face density')
        tau = Variable('\\tau', 'N*m', 'fan shaft torque')
        Utip = Variable('U_{\\rm tip}', 'm/s', 'fan tip speed')
        phi = Variable('\\phi', '-', 'flow coefficient')
        psi = Variable('\\psi', '-',
            'stagnation enthalpy rise coefficient')
        Omega = Variable('\\Omega', 'rpm', 'fan angular speed')

        constraints += [
            Utip == Omega * fan['r_{\\rm tip}'],
            phi == V2 / Utip,
            phi == M2 / Mtip,
            psi == dht / Utip**2,
            P_shaft == tau * Omega,
            pt2 / p2 == (Tt2 / T2)**(gamma / (gamma-1)),
            p2 == rho2 * R * T2,
            mdot == rho2 * V2 * fan['A_2'],
            M2 == V2 / (gamma * R * T2)**0.5,
        ]

        with SignomialsEnabled():
            constraints += [
                SignomialEquality(
                    dht,
                    Utip**2 - Utip*V2*fan['\\tan\\beta']
                ),
                SignomialEquality(
                    Tt2 / T2, 1 + (gamma - 1) / 2 * M2**2
                )
            ]

        return constraints


class PropulsorSizing(Model):
    def setup(
        self,
        h_m, M0, hdot_fpm, L_D, CD_p_fuse, CD_p_wing, D_hx_lbf, m_lb,
        n_fuse_fans=-1, n_wing_fans=-1,
        fuse_bli='off', wing_bli='off',
        same_fan=True
    ):

        # fan sizing model(s)
        fuse_fan = self.fuse_fan = Fan()
        fuse_fan.substitutions.update({
            'c_{d,{\\rm nace}}': 0.005,
            'r_{\\rm hub}/r_{\\rm tip}': 0.38,
            '\\tan\\beta': tan(45 * pi / 180),
        })
        if n_fuse_fans >= 0:
            if n_fuse_fans == 0:
                n_fuse_fans = 1e-32
                same_fan = False
            N_fuse_fans = Variable('N_{\\rm fans,fuse}', n_fuse_fans, '-',
                'number of fuselage fans')
        else:
            N_fuse_fans = Variable('N_{\\rm fans,fuse}', '-',
                'number of fuselage fans')

        wing_fan = self.wing_fan = Fan()
        wing_fan.substitutions.update({
            'c_{d,{\\rm nace}}': 0.005,
            'r_{\\rm hub}/r_{\\rm tip}': 0.38,
            '\\tan\\beta': tan(45 * pi / 180),
        })
        if n_wing_fans >= 0:
            if n_wing_fans == 0:
                n_wing_fans = 1e-12
                same_fan = False
            N_wing_fans = Variable('N_{\\rm fans,wing}', n_wing_fans, '-',
                'number of wing fans')
        else:
            N_wing_fans = Variable('N_{\\rm fans,wing}', '-',
                'number of wing fans')

        # aircraft performance parameters (inputs), fan performance model(s)
        with Vectorize(len(hdot_fpm)):
            D_nace_fuse = Variable('D_{\\rm nace,fuse}', 'lbf',
                'fuselage propulsor nacelle drag')
            D_nace_wing = Variable('D_{\\rm nace,wing}', 'lbf',
                'wing propulsor nacelle drag')
            D = Variable('D', 'lbf', 'drag')
            D_mfan = Variable('D_{m_{\\rm fan}}', 'lbf',
                'drag increment due to fan mass')
            hdot = Variable('\\dot{h}', hdot_fpm, 'ft/min', 'climb rate')
            D_hx = Variable('D_{\\rm HX}', D_hx_lbf, 'lbf', 'heat exchanger drag')
            L_D = Variable('L/D', L_D, '-', 'airframe lift-to-drag ratio')
            m_aircraft = Variable('m_{\\rm aircraft}', m_lb, 'lb',
                'aircraft mass')

            if fuse_bli in ['on']:
                CD_p_fuse = Variable('C_{D_{p,{\\rm fuse}}}', CD_p_fuse, '-',
                    'fuselage profile drag coefficient')
                pi_BLI_fuse = Variable('\\pi_{\\rm BLI,fuse}', '-',
                    'fuselage BLI pressure ratio')

            if wing_bli in ['top', 'bottom', 'both']:
                CD_p_wing = Variable('C_{D_{p,{\\rm wing}}}', CD_p_wing, '-',
                    'wing profile drag coefficient')
                pi_BLI_wing = Variable('\\pi_{\\rm BLI,wing}', '-',
                    'wing BLI pressure ratio')

            F_tot = Variable('F_{\\rm tot}', 'lbf',
                'total propulsion system net thrust')
            f_wing = Variable('f_{\\rm wing}', '-',
                'wing propulsor thrust fraction')
            f_fuse = Variable('f_{\\rm fuse}', '-',
                'fuselage propulsor thrust fraction')

            flight_state = FlightState()
            fuse_fan_perf = self.fuse_fan_perf \
                = fuse_fan.performance(flight_state)
            wing_fan_perf = self.wing_fan_perf \
                = wing_fan.performance(flight_state)

        # set flight state parameters based on altitude and Mach number
        for ii in range(len(h_m)):
            p, T = isa(h_m[ii])
            Tt = T * (1 + (gamma - 1) / 2 * M0[ii]**2)
            pt = p * (Tt / T)**(gamma / (gamma - 1))
            flight_state.substitutions.update({
                flight_state['M'][ii]: M0[ii],
                flight_state['p'][ii]: p,
                flight_state['p_t'][ii]: pt,
                flight_state['T'][ii]: T,
                flight_state['T_t'][ii]: Tt
            })

        # component model constraints
        constraints = [
            wing_fan, fuse_fan,
            wing_fan_perf, fuse_fan_perf,
            flight_state
        ]

        # Boundary layers
        b_fuse = Variable('b_{\\rm fuse}', 'ft', 'fuselage TE span')
        b_wing = Variable('b_{\\rm wing}', 'ft', 'wing span')
        f_b_wing = Variable('f_{b_{\\rm wing}}', '-',
            'usable wing span fraction')
        S_ref = Variable('S_{\\rm ref}', 'ft^2', 'reference area')
        
        if fuse_bli in ['on']:
            f_BLI_fuse = Variable('f_{\\rm BLI,fuse}', '-', 'fuselage BLI fraction')
            f_surf_fuse = Variable('f_{\\rm surf,fuse}', 0.89, '-',
                'fuselage surface dissipation fraction')
            f_top_fuse = Variable('f_{\\rm top,fuse}', 0.5, '-',
                'fuselage top surface dissipation fraction')
            D_p_fuse = CD_p_fuse * flight_state['q'] * S_ref
            f_nace_fuse = 0.5

            constraints += [
                f_BLI_fuse == N_fuse_fans * fuse_fan['d_{\\rm fan}'] / b_fuse \
                    * f_top_fuse,
                fuse_fan_perf['p_{t1}'] \
                    + f_surf_fuse * f_BLI_fuse * D_p_fuse * flight_state['V'] \
                    * flight_state['\\rho'] / fuse_fan_perf['\\dot{m}'] / N_fuse_fans \
                    <= flight_state['p_t'],
                pi_BLI_fuse == fuse_fan_perf['p_{t1}'] / flight_state['p_t'],
                N_fuse_fans * fuse_fan['d_{\\rm fan}'] <= b_fuse, # TODO only for BLI?
            ]
        else:
            f_nace_fuse = 1
            constraints += [fuse_fan_perf['p_{t1}'] == flight_state['p_t']]

        if wing_bli in ['top', 'bottom', 'both']:
            f_BLI_wing = Variable('f_{\\rm BLI,wing}', '-', 'wing BLI fraction')
            f_surf_wing = Variable('f_{\\rm surf,wing}', 0.89, '-',
                'wing surface dissipation fraction')
            D_p_wing = CD_p_wing * flight_state['q'] * S_ref
            f_nace_wing = 0.5   

            constraints += [
                wing_fan_perf['p_{t1}'] \
                    + f_surf_wing * f_BLI_wing * D_p_wing * flight_state['V'] \
                    * flight_state['\\rho'] / wing_fan_perf['\\dot{m}'] \
                    / N_wing_fans <= flight_state['p_t'],
                pi_BLI_wing == wing_fan_perf['p_{t1}'] / flight_state['p_t'],
            ]

            if wing_bli == 'top':
                f_top_wing = Variable('f_{\\rm top,wing}', 0.64, '-',
                    'wing top surface dissipation fraction')
                constraints += [
                    f_BLI_wing == N_wing_fans * wing_fan['d_{\\rm fan}'] \
                        / b_wing * f_top_wing,
                ]
            elif wing_bli == 'bottom':
                f_bot_wing = Variable('f_{\\rm bot,wing}', 0.36, '-',
                    'wing bottom surface dissipation fraction')
                constraints += [
                    f_BLI_wing == N_wing_fans * wing_fan['d_{\\rm fan}'] \
                        / b_wing * f_bot_wing,
                ]
            elif wing_bli == 'both':
                constraints += [
                    f_BLI_wing == N_wing_fans * wing_fan['d_{\\rm fan}'] \
                        / b_wing,
                ]
        else:
            f_nace_wing = 1
            constraints += [wing_fan_perf['p_{t1}'] == flight_state['p_t']]

        constraints += [
            N_wing_fans * wing_fan['d_{\\rm fan}'] <= b_wing * f_b_wing, # TODO for BLI and non-BLI?
        ]

        # Nacelle drag
        with SignomialsEnabled():
            constraints += [
                D_nace_fuse / f_nace_fuse >= flight_state['q'] \
                    * (fuse_fan['c_{d,{\\rm nace}}'] \
                    * fuse_fan['\\ell_{\\rm nace}'] \
                    * fuse_fan['d_{\\rm fan}'] \
                    * (pi + 2 * (N_fuse_fans - 1))),
                D_nace_wing / f_nace_fuse >= flight_state['q'] \
                    * (wing_fan['c_{d,{\\rm nace}}'] \
                    * wing_fan['\\ell_{\\rm nace}'] * 2 \
                    * wing_fan['d_{\\rm fan}'] \
                    * (pi + 2 * (N_wing_fans / 2 - 1)))
            ]

        # airframe "drag"
        if fuse_bli == 'off' and wing_bli == 'off':
            constraints += [
                D >= m_aircraft * g / L_D + D_hx + D_mfan + D_nace_fuse + D_nace_wing
            ]
        elif fuse_bli == 'on' and wing_bli in ['top', 'bottom', 'both']:
            with SignomialsEnabled():
                constraints += [
                    D >= m_aircraft * g / L_D + D_hx + D_mfan + D_nace_fuse + D_nace_wing \
                        - f_BLI_fuse * D_p_fuse - f_BLI_wing * D_p_wing
                ]
        elif fuse_bli == 'on':
            with SignomialsEnabled():
                constraints += [
                    D >= m_aircraft * g / L_D + D_hx + D_mfan + D_nace_fuse + D_nace_wing \
                        - f_BLI_fuse * D_p_fuse
                ]
        elif wing_bli in ['top', 'bottom', 'both']:
            with SignomialsEnabled():
                constraints += [
                    D >= m_aircraft * g / L_D + D_hx + D_mfan + D_nace_fuse + D_nace_wing \
                        - f_BLI_wing * D_p_wing
                ]

        # force balance: lift >= weight, "thrust" >= "drag"
        constraints += [
            f_wing / wing_fan_perf['\\dot{m}'] / N_wing_fans * F_tot \
                + flight_state['V'] <= wing_fan_perf['u_8'],
            f_fuse / fuse_fan_perf['\\dot{m}'] / N_fuse_fans * F_tot \
                + flight_state['V'] <= fuse_fan_perf['u_8'],
            F_tot >= D + (m_aircraft \
                + N_wing_fans * wing_fan['m_{\\rm fan}'] \
                + N_fuse_fans * fuse_fan['m_{\\rm fan}']) \
                * g * hdot / flight_state['V'],
            D_mfan >= (wing_fan['m_{\\rm fan}'] * N_wing_fans \
                + fuse_fan['m_{\\rm fan}'] * N_fuse_fans) * g / L_D,
        ]

        # total thrust SP
        with SignomialsEnabled():
            constraints += [f_wing + f_fuse >= 1]

        if same_fan:
            constraints += [
                wing_fan['r_{\\rm tip}'] == fuse_fan['r_{\\rm tip}'],
            ]

        # TODO Motor model can go here
        #   set up motor sizing model
        #   set up motor performance model (vectorize like fan performance)
        #   sizing constraint(s)
        #   performance constraints (speed and power)

        # minimize weighted average of total propulsive shaft power
        # TODO change to motor electrical power
        P_wing = wing_fan_perf['P_{\\rm shaft}']
        P_fuse = fuse_fan_perf['P_{\\rm shaft}']
        self.cost = \
            N_wing_fans * (0.9 * P_wing[-1] \
            + 0.1 / (len(P_wing) - 1) * sum(P_wing[:-1])) \
            + N_fuse_fans * (0.9 * P_fuse[-1] \
            + 0.1 / (len(P_fuse) - 1) * sum(P_fuse[:-1]))

        # SP initialization
        self.x0 = {
            # fan stagnation enthalpy rise SP
            wing_fan_perf['M_2']: 0.5,
            wing_fan_perf['T_{t8}']: 288 * units('K'),
            wing_fan_perf['\\Delta h_t']: 1e12 * units('J/kg'),
            fuse_fan_perf['M_2']: 0.5,
            fuse_fan_perf['T_{t8}']: 288 * units('K'),
            fuse_fan_perf['\\Delta h_t']: 1e12 * units('J/kg'),
        }

        return constraints


if __name__ == '__main__':
    from numpy import array
    from cheeta_performance import h_m, M0, hdot_fpm, L_D, CD_p_fuse, CD_p_wing, D_hx_lbf, m_lb

    model = PropulsorSizing(
        h_m, M0, hdot_fpm, L_D, CD_p_fuse, CD_p_wing, D_hx_lbf, m_lb,
        n_fuse_fans=3, n_wing_fans=6, fuse_bli='on', wing_bli='top',
        same_fan=True
    )

    sol = model.localsolve(x0=model.x0, verbosity=2, iteration_limit=1000)
    print(sol.table())