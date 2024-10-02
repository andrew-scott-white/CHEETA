# fan BLI optimization script by Durgesh Chandel for NASA CHEETA
# D. Chandel et al., “Conceptual Design of Distributed Electrified Boundary Layer Ingesting Propulsors for the CHEETA Aircraft Concept,” 2021, doi: 10.2514/6.2021-3287.
# D. Chandel et al., “Fan and Motor Co-optimization for a Distributed Electric Aircraft Propulsion System,” IEEE Trans. Transp. Electrif., pp. 1–11, 2022, doi: 10.1109/TTE.2022.3204202.

# Adapted by Andrew Scott White on 4 October 2022

import numpy as np
from numpy import array

# aircraft scalar inputs from Elias
span_fuse = 14 # ft
span_wing = 118 # ft
S_ref = 2000 # sq ft

# aircraft vector inputs from Elias
h_m = array([0, 14934.72218, 29623.40467, 36937.52169, 37000]) * 0.3048
M0 = array([0.25, 0.398315526, 0.546062782, 0.771144278, 0.773])
hdot_fpm = array([988.7660847, 1100.035258, 951.6387309, 141.9973407, 0])
CL = array([1.05135995, 0.728831109, 0.723232455, 0.508639028, 0.507563935])
CD = array([0.063728775, 0.037345509, 0.037288927, 0.027149456, 0.027137671])
CD_p_fuse = array([0.007739158, 0.007622935, 0.007670017, 0.00739524, 0.007394066])
CD_p_wing = array([0.004306984, 0.004608784, 0.00493297, 0.005237413, 0.005240396])
D_hx_lbf = array([3373.925121, 2790.819298, 3217.626928, 3848.654638, 962.1636596])
m_lb = array([194681.2805, 193840.7004, 193001.0401, 192162.2906, 192162.2906]) - 14380
L_D = CL / CD


# -------------------------------------------------------------


if __name__ == '__main__':
    from gpkit import units
    from fans_dchandel import PropulsorSizing
    from numpy import zeros

    n_fuse_fans = array([0.0, 3])
    n_wing_fans = array([2.0, 6])
    fuse_bli = array(['off', 'on'])
    wing_bli = array(['off', 'top'])

    p_cruise_mw = zeros(len(n_fuse_fans))
    psc = zeros(len(n_fuse_fans))
    p_max_fuse_mw = zeros(len(n_fuse_fans))
    p_max_wing_mw = zeros(len(n_fuse_fans))
    d_fan_fuse_in = zeros(len(n_fuse_fans))
    d_fan_wing_in = zeros(len(n_fuse_fans))
    fpr_fuse_min = zeros(len(n_fuse_fans))
    fpr_fuse_max = zeros(len(n_fuse_fans))
    fpr_wing_min = zeros(len(n_fuse_fans))
    fpr_wing_max = zeros(len(n_fuse_fans))
    f_fuse_BLI   = zeros(len(n_fuse_fans))
    f_wing_BLI   = zeros(len(n_fuse_fans))

    for ii in range(len(n_fuse_fans)):
        if n_fuse_fans[ii] == 0 or n_wing_fans[ii] == 0: same_fan = False
        else: same_fan = True
        model = PropulsorSizing(
            h_m, M0, hdot_fpm, L_D, CD_p_fuse, CD_p_wing, D_hx_lbf, m_lb,
            n_fuse_fans=n_fuse_fans[ii], n_wing_fans=n_wing_fans[ii],
            fuse_bli=fuse_bli[ii], wing_bli=wing_bli[ii], same_fan=same_fan
        )
        
        model.substitutions.update({
            'b_{\\rm fuse}': span_fuse*units('ft'), # fuselage TE span
            'b_{\\rm wing}': span_wing*units('ft'), # wing span
            'f_{b_{\\rm wing}}': 0.5, # usable wing span fraction (a.s.w. - assumption of wing fraction that can support propulsors?)
            'S_{\\rm ref}': S_ref*units('ft**2'), # reference area (a.s.w. - aircraft surface area?)
            })
        sol = model.localsolve(x0=model.x0, verbosity=1, iteration_limit=1000, rel_tol=1e-9)

        n_fuse_fans[ii] = sol['variables']['N_{\\rm fans,fuse}']
        n_wing_fans[ii] = sol['variables']['N_{\\rm fans,wing}']
        p_cruise_mw[ii] \
            = n_fuse_fans[ii] \
            * sol['variables'][model.fuse_fan_perf['P_{\\rm shaft}']][-1] \
            + n_wing_fans[ii] \
            * sol['variables'][model.wing_fan_perf['P_{\\rm shaft}']][-1]
        psc[ii] = (1 - p_cruise_mw[ii] / p_cruise_mw[0])*100
        p_max_fuse_mw[ii] \
            = max(sol['variables'][model.fuse_fan_perf['P_{\\rm shaft}']])
        p_max_wing_mw[ii] \
            = max(sol['variables'][model.wing_fan_perf['P_{\\rm shaft}']])
        d_fan_fuse_in[ii] = sol['variables'][model.fuse_fan['d_{\\rm fan}']]
        d_fan_wing_in[ii] = sol['variables'][model.wing_fan['d_{\\rm fan}']]
        fpr_fuse_max[ii] \
            = max(sol['variables'][model.fuse_fan_perf['\\pi_{\\rm fan}']])
        fpr_fuse_min[ii] \
            = min(sol['variables'][model.fuse_fan_perf['\\pi_{\\rm fan}']])
        fpr_wing_max[ii] \
            = max(sol['variables'][model.wing_fan_perf['\\pi_{\\rm fan}']])
        fpr_wing_min[ii] \
            = min(sol['variables'][model.wing_fan_perf['\\pi_{\\rm fan}']])
        if ii > 0:
            f_wing_BLI[ii] = sol['variables'][model['f_{\\rm BLI,wing}']]
            f_fuse_BLI[ii] = sol['variables']['f_{\\rm BLI,fuse}']
            
        #D_p_wing = CD_p_wing * flight_state['q'] * S_ref
        
        


    print('                                 & Non-BLI Twin          & CHEETA')
    print('Cruise Power [MW]                & {0:5.2f}              & {1:5.2f}                                  \\\\'.format(
                                            p_cruise_mw[0],         p_cruise_mw[1]                              ))
    print('Power Savings Coefficient        & ---                   & {1:4.1f}\\%                               \\\\'.format(
                                            psc[0],                 psc[1]                                      ))
    print('Fan diameter [in]                & {0:5.1f}              & {1:4.1f}                                  \\\\'.format(
                                            d_fan_wing_in[0],       d_fan_wing_in[1]                            ))
    print('No. of wing fans                 &  2                    & {1:4.1f}                                  \\\\'.format(
                                            n_wing_fans[0],         n_wing_fans[1]                              ))
    print('No. of fuselage fans             & ---                   & {1:4.1f}                                  \\\\'.format(
                                            n_fuse_fans[0],         n_fuse_fans[1]                              ))
    print('Max motor power [MW]             & {0:5.2f}              & {1:5.2f}                                  \\\\'.format(
                                            p_max_wing_mw[0],       max(p_max_wing_mw[1], p_max_fuse_mw[1])     ))
    print('FPR range                        & {0:4.2f} -- {1:4.2f}  & {2:4.2f} -- {3:4.2f}                      \\\\'.format(
                                            fpr_wing_min[0], 
                                            fpr_wing_max[0],
                                                                    min(fpr_wing_min[1], fpr_fuse_min[1]), 
                                                                    max(fpr_wing_max[1], fpr_fuse_max[1]),      ))
    print('Wing BLI Fraction                     & {0:4.3f}              & {1:4.3f}                                  \\\\'.format(
                                            f_wing_BLI[0],          f_wing_BLI[1]                               ))
    print('Fuselage BLI Fraction                 & {0:4.3f}              & {1:4.3f}                                  \\\\'.format(
                                            f_fuse_BLI[0],          f_fuse_BLI[1]                               ))
    
 
    
    