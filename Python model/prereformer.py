from enthalpy import const, enthalpy, enthalpy_heavy, heavy_const
from Keq import K_smr, K_wgsr
from cp import cp_const, cp, single_cp, H2Ocoeff
import scipy.optimize as opt
import autograd.numpy as np
from autograd import grad, jacobian

def PR_Reactor(inlet, outlet):
    n4CH4, n4H2O, n4H2, n4CO, n4CO2, n1C2H6, n1C3H8, n1n_C4H10, n1i_C4H10, n1_C5,T4 = inlet
    n5CH4, n5H2O, n5H2, n5CO, n5CO2, T5 = outlet

    n4stream = inlet[:10]
    n5 = n5CH4 + n5H2O + n5H2 + n5CO + n5CO2

    pr_ksi1 = n1C2H6
    pr_ksi2 = n1C3H8
    pr_ksi3 = n1n_C4H10
    pr_ksi4 = n1i_C4H10
    pr_ksi5 = n1_C5
    pr_ksi6 = n4CH4 - n5CH4
    pr_ksi7 = n5CO2 - n4CO2

    eq1 = n5H2O - n4H2O + 2 * pr_ksi1 + 3 * pr_ksi2 + 4 * pr_ksi3 + 4 * pr_ksi4 + 5 * pr_ksi5 + pr_ksi6 + pr_ksi7
    eq2 = n5H2 - n4H2 - 5 * pr_ksi1 - 7 * pr_ksi2 - 9 * pr_ksi3 - 9 * pr_ksi4 - 11 * pr_ksi5 - 3 * pr_ksi6 - pr_ksi7
    eq3 = n5CO - n4CO - 2 * pr_ksi1 - 3 * pr_ksi2 - 4 * pr_ksi3 - 4 * pr_ksi4 - 5 * pr_ksi5 - pr_ksi6 + pr_ksi7
    eq4 = K_smr(T5) * ((n5CH4 / n5) * (n5H2O / n5)) - (((n5CO / n5) * (n5H2 / n5) ** 3))
    eq5 = K_wgsr(T5) * ((n5CO / n5) * (n5H2O / n5)) - (((n5CO2 / n5) * (n5H2 / n5)))
    eq6 = np.dot(n4stream, enthalpy_heavy(heavy_const, T4)) - np.dot(outlet[:5], enthalpy(const, T5))

    return np.array([eq1,eq2,eq3,eq4,eq5,eq6])


inletstream = np.array([2685.39, 8444.18, 0.00, 0.00, 46.77, 208.86, 229.41, 84.92, 48.28, 126.69,693.00])
guess = np.array([4026.22, 7561.41, 1150.02, 975.24, 0.52, 700.0])

def function(outlet):
    return PR_Reactor(inletstream,outlet)

def solve(f,guess):
    sol = opt.fsolve(f,guess)
    print(sol)

#solve(function,guess)