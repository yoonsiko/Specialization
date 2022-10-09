# Author: Yoonsik Oh
# Project: SUBPRO Blue hydrogen
# https://www.sciencedirect.com/science/article/pii/S0360319913014456

import numpy as np
import autograd.numpy as np


def K_smr(T):
    return (np.exp(-22790 / T + 8.156 * np.log(T) - 4.421 / 10 ** 3 * T
                   - 4.330 * 10 ** 3 / (T ** 2) - 26.030))

"""
def K_wgsr(T):
    return (np.exp(-5087 / T + 1.560 * np.log(T) - 1.509 / 10 ** 4 * T
                   - 4.762 * 10 ** 4 / (T ** 2) - 13.9333))
"""

def K_wgsr(T):
    return np.exp(5693.5/T + 1.077*np.log(T) + 5.44e-4*T - 1.125e-7*T**2 - 49170/(T**2)-13.148)


#print(K_wgsr(700))
print(K_smr(1073))
print(K_wgsr(1073))
# WGSR: https://www.researchgate.net/figure/Equilibrium-constants-for-WGSR-Twigg-1989_tbl5_228506095