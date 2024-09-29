import numpy as np

def SensSort(SENS):
    """Function for converting array sensitivities to sclar sensitivities by assigning maximum value as the sens
    """

    for kk in range(len(SENS)):
        if hasattr(SENS[kk], "__len__"):  
            
            if len(SENS[kk]) > 1:
                sens_max = np.max(SENS[kk])
                sens_min = np.min(SENS[kk])
                
                if abs(sens_max) > abs(sens_min):
                    SENS[kk] = sens_max
                else:
                    SENS[kk] = sens_min
            
    return SENS