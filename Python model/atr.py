# Author: Yoonsik Oh
# Project: SUBPRO Blue hydrogen

#import numpy as np
from enthalpy import const, enthalpy
from Keq import K_smr, K_wgsr
import scipy.optimize as opt
import autograd.numpy as np
from autograd import grad, jacobian


# inlet/outlet = [molar flow, T, P]
# molar flow = CH4, H2O, H2, CO, CO2
Mm = np.array([16, 18, 2, 28, 44])

def ATR_Reactor(inlet, outlet):
    n0CH4, n0H2O, n0H2, n0CO, n0CO2, P0, T0= inlet
    nO2, nCH4, nH2O, nH2, nCO, nCO2, P = outlet

    ksi1 = n0CH4 - nCH4
    ksi2 = nCO2 - n0CO2
    ksi3 = -nO2
    nH2O += 2*ksi3
    nH2 += -2*ksi3
    ntot = nCH4 + nH2O + nH2 + nCO + nCO2

    T = 1600

    eq1 = K_smr(T)*((nCH4/ntot) * (nH2O/ntot )) - (((nCO/ntot) * (nH2/ntot) ** 3))

    eq2 = K_wgsr(T)*((nCO/ntot) * (nH2O/ntot)) - (((nCO2/ntot) * (nH2/ntot)))

    eq3 = nH2O - n0H2O + ksi1 + ksi2

    eq4 = nH2 - n0H2 - 3 * ksi1 - ksi2

    eq5 = nCO - n0CO - ksi1 + ksi2

    eq6 = np.dot(inlet[:5],enthalpy(const, T0)) - np.dot(outlet[1:6],enthalpy(const, T)) + nO2*44.2

    eq7 = P-P0

    balances = np.array([eq1,eq2,eq3,eq4,eq5,eq6,eq7])

    return balances

# https://www.ohio.edu/mechanical/thermo/property_tables/combustion/oxygen_enth.html
inletstream = np.array([5.834e+01, 1.384e+03, 2.907e+03,
 9.897e+02, 8.457e-02, 30, 973])
outletstream = np.array([5.83315576e+01, 1.38399198e+03, 2.90740371e+03,
 9.89702263e+02, 8.45788845e-02, 3.00000000e+01,30])


def function(outlet):

    return ATR_Reactor(inletstream,outlet)

def solve(f,guess,grad_function):
    sol = opt.fsolve(f,guess,fprime=grad_function)
    print(sol)
