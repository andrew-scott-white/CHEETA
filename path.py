import numpy as np

# Design mission for CHEETA
        
h_alt = np.asarray(     [0,      5000,   10000,  15000,  20000,  25000,  30000,  35000,     36000,  36867.5,  37000])                     # altitude (ft)
hdot_in = np.asarray(   [1142,   1186,   1162,   1080,   969,    849,    685,    484,       433,    370,    0])*.00508                  # climb rate (m/s)
M = np.asarray(         [0.25,   0.3,    0.35,   0.40,   0.45,   0.50,   0.55,   0.65,      .68,    .693,   0.773])                     # flight Mach number 
P = np.asarray(         [101325, 84309,  69692,  57206,  46602,  37652,  30151,  23911,     22753,  21876,  21876])                     # ambient static pressure (Pa)
T = np.asarray(         [311,    278,    268,    258,    249,    239,    229,    219,       218,    218,    218])                       # ambient static temperature (K)
rho = np.asarray(       [23.77,  20.48,  17.56,  14.96,  12.67,  10.66,  8.91,   7.38,      7.078,  6.776,  6.625])*10**(-4)*515.379    # ambient air density (kg/m^3)
mu = np.asarray(        [3.737,  3.637,  3.534,  3.430,  3.324,  3.217,  3.107,  2.995,     2.99,   2.9846,   2.98])*10**(-7)*47.8803    # ambient air dynamic viscosity (kg/m/s)
Cp = np.asarray(        [1006,   1006,   1006,   1006,   1006,   1006,   1006,   1006,      1006,   1006,   1006])                      # air specific heat (J/kg/K)
C_D = np.asarray(       [.0624,  .0473,  .0406,  .0376,  .0362,  .0366,  .0374,  .0334,     .033,   .033,  .0283])                     # aircraft drag coefficient
C_D_wing = np.asarray(  [.00431, .00431, .00431, .00461, .00461, .00461, .00493, .00493,    .00493, .00524, .00524])
C_D_fuse = np.asarray(  [.00774, .00770, .00766, .00762, .00764, .00766, .00767, .00756,    .00748, .00740, .00740])

C_L = np.asarray(       [1.060,  .883,   .787,   .735,   .709,   .714,   .726,   .650,      .637,   .635,   .514])                      # aircraft lift coefficient
LDR = C_L/C_D

        
