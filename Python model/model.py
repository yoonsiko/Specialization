# Author: Yoonsik Oh
# Project: SUBPRO Blue hydrogen

from enthalpy import const, enthalpy, enthalpy_heavy, heavy_const
from Keq import K_smr, K_wgsr
import scipy.optimize as opt
import autograd.numpy as np
from autograd import grad, jacobian

# inlet/outlet = [molar flow, P, T/Q/nO2]
C = 12
H = 1.
O = 16
# molar flow component order = CH4, H2O, H2, CO, CO2
Mm = np.array([C * 1 + H * 4, H * 2 + O, H * 2, C + O, C + O * 2])  # This is to check the conservation of mass
# molar flow component order = CH4, H2O, H2, CO, CO2, C2H6, C3H8, n-C4H10, i-C4H10, C5H12

Mm_heavy = np.array(
    [C * 1 + H * 4, H * 2 + O, H * 2, C + O, C + O * 2, C * 2 + H * 6, C * 3 + H * 8, C * 4 + H * 10, C * 4 + H * 10,
     C * 5 + H * 12])

inletcomposition = np.array([78.24, 0.0, 0.0, 0.0, 1.34, 6.10, 6.7, 2.48, 1.41, 3.7]) / 100 # [Sm^3/h]
inletstream_kg = 4000 * 0.829
avgMm = np.dot(Mm_heavy,inletcomposition)
inletstream_mol = inletstream_kg/avgMm

inletstream = np.multiply(inletstream_mol,inletcomposition)
other = np.array([30.0, 311.0])

inletstream = np.append(inletstream, other)  # This is the given values


guess = np.array([400, 30.0, 10000,  # n2
                  100, 400, 0.0, 0.0, 1, 30.0, 400.0,  # n3
                  100, 400, 0.0, 0.0, 1, 30.0, 4000.0,  # n4
                  200, 400, 100, 1, 10, 30.0, 700.0,  # n5
                  200, 300, 100.0, 1, 10, 30.0, 20000.0,  # n6
                  10, 300, 500.0, 100.0, 1.0, 30.0, -5000,  # n7
                  1, 200, 500.0, 100.0, 1.0, 30.0, 30,  # n8
                  1, 200, 500.0, 100.0, 1.0, 30.0, 200,  # n9
                  1, 10, 500.0, 200, 1.0, 30.0, -6000.0,  # n10
                  1, 10, 500.0, 10.0, 1.0, 30.0,  # n11
                  1, 1.0, 500.0, 10.0, 100, 30.0, 313.0,  # n12
                  30.0, 313.0,  # n13
                  1.0, 1.0, 500.0, 1.0, 1.0, 30.0, 313.0,  # n14
                  30.0, 313.0])  # n15

splitratio = np.array([[1.0, 0.001, 1.0, 1.0, 1.0], [0.01, 0.01, 0.999, 0.01, 0.01]])


# Input is given, Q and output is the initial guess values
def system(input, output, f):
    # Defining and unpacking all the variables
    n1CH4, n1H2O, n1H2, n1CO, n1CO2, \
    n1C2H6, n1C3H8, n1i_C4H10, n1n_C4H10, n1_C5, P1, T1 = input
    f1, f2 = f
    # Note: on stream 8, we set the reactor temperature depending on O2 stream
    n2H2O, P2, Q1, \
    n3CH4, n3H2O, n3H2, n3CO, n3CO2, P3, T3, \
    n4CH4, n4H2O, n4H2, n4CO, n4CO2, P4, Q2, \
    n5CH4, n5H2O, n5H2, n5CO, n5CO2, P5, T5, \
    n6CH4, n6H2O, n6H2, n6CO, n6CO2, P6, Q3, \
    n7CH4, n7H2O, n7H2, n7CO, n7CO2, P7, Q4, \
    n8CH4, n8H2O, n8H2, n8CO, n8CO2, P8, nO2, \
    n9CH4, n9H2O, n9H2, n9CO, n9CO2, P9, T9, \
    n10CH4, n10H2O, n10H2, n10CO, n10CO2, P10, Q5, \
    n11CH4, n11H2O, n11H2, n11CO, n11CO2, P11, \
    n12CH4, n12H2O, n12H2, n12CO, n12CO2, P12, T12, \
    P13, T13, \
    n14CH4, n14H2O, n14H2, n14CO, n14CO2, P14, T14, \
    P15, T15 = output

    # Defining variables
    n2 = n2H2O
    n5 = n5CH4 + n5H2O + n5H2 + n5CO + n5CO2
    n7 = n7CH4 + n7H2O + n7H2 + n7CO + n7CO2
    n8 = n8CH4 + n8H2O + n8H2 + n8CO + n8CO2
    n10 = n10CH4 + n10H2O + n10H2 + n10CO + n10CO2

    nCarbon = n1CH4 + n1C2H6 + n1C3H8 + n1n_C4H10 + n1i_C4H10 + n1_C5

    # BLOCK 1: mixing purified natural gas and steam and setting S/C to 2.5
    stcr = 2.5  # steam to carbon ratio (molar ratio)
    eq0 = n3H2O - stcr * nCarbon
    eq1 = n3CH4 - n1CH4
    eq2 = n3H2O - n2H2O - n1H2O
    eq3 = n3H2 - n1H2
    eq4 = n3CO - n1CO
    eq5 = n3CO2 - n1CO2
    eq6 = P2 - P1
    eq7 = P3 - P2
    n3stream = output[3:8]
    n3stream = np.append(n3stream, input[5:10])
    T2 = 423
    # eq8 = np.dot(input[:10], cp(cp_const, T1)) * 8.314 * (T1 - 298) + n2 * 8.314 * single_cp(H2Ocoeff, T2) * (T2 - 298) \
    #      - np.dot(n3stream, cp(cp_const, T3)) * 8.314 * (T3 - 298)

    eq8 = np.dot(input[:10], enthalpy_heavy(heavy_const, T1)) + \
          np.dot([0, n2H2O, 0, 0, 0, 0, 0, 0, 0, 0], enthalpy_heavy(heavy_const, T2)) - \
          np.dot(n3stream, enthalpy_heavy(heavy_const, T3))

    # BLOCK 2: Pre-prereformer heat exchanger
    T4 = 693
    eq9 = n4CH4 - n3CH4
    eq10 = n4H2O - n3H2O
    eq11 = n4H2 - n3H2
    eq12 = n4CO - n3CO
    eq13 = n4CO2 - n3CO2
    n4stream = output[10:15]
    n4stream = np.append(n4stream, input[5:10])
    eq14 = np.dot(n3stream, enthalpy_heavy(heavy_const, T3)) - \
           np.dot(n4stream, enthalpy_heavy(heavy_const, T4)) + Q1
    eq15 = P4 - P3

    # BLOCK 3: Prereformer (Removing all heavier carbons to methane, (full conversion))
    pr_ksi1 = n1C2H6
    pr_ksi2 = n1C3H8
    pr_ksi3 = n1n_C4H10
    pr_ksi4 = n1i_C4H10
    pr_ksi5 = n1_C5
    pr_ksi6 = n4CH4 - n5CH4
    pr_ksi7 = n5CO2 - n4CO2

    eq16 = n5H2O - n4H2O + 2 * pr_ksi1 + 3 * pr_ksi2 + 4 * pr_ksi3 + 4 * pr_ksi4 + 5 * pr_ksi5 + pr_ksi6 + pr_ksi7
    eq17 = n5H2 - n4H2 - 5 * pr_ksi1 - 7 * pr_ksi2 - 9 * pr_ksi3 - 9 * pr_ksi4 - 11 * pr_ksi5 - 3 * pr_ksi6 - pr_ksi7
    eq18 = n5CO - n4CO - 2 * pr_ksi1 - 3 * pr_ksi2 - 4 * pr_ksi3 - 4 * pr_ksi4 - 5 * pr_ksi5 - pr_ksi6 + pr_ksi7
    eq19 = K_smr(T5) * ((n5CH4 / n5) * (n5H2O / n5)) - (((n5CO / n5) * (n5H2 / n5) ** 3))
    eq20 = K_wgsr(T5) * ((n5CO / n5) * (n5H2O / n5)) - (((n5CO2 / n5) * (n5H2 / n5)))
    eq21 = np.dot(n4stream, enthalpy_heavy(heavy_const, T4)) - np.dot(output[17:22], enthalpy(const, T5))
    eq22 = P5 - P4

    # BLOCK 4: Pre-GHR Heat Exchanger
    T6 = 973
    eq23 = n6CH4 - n5CH4
    eq24 = n6H2O - n5H2O
    eq25 = n6H2 - n5H2
    eq26 = n6CO - n5CO
    eq27 = n6CO2 - n5CO2
    eq28 = np.dot(output[17:22], enthalpy(const, T5)) - \
           np.dot(output[24:29], enthalpy(const, T6)) + Q2
    eq29 = P6 - P5

    # BLOCK 5: GHR
    T7 = 1073
    ghr_ksi1 = n6CH4 - n7CH4
    ghr_ksi2 = n7CO2 - n6CO2
    eq30 = K_smr(T7) * ((n7CH4 / n7) * (n7H2O / n7)) - (((n7CO / n7) * (n7H2 / n7) ** 3))
    eq31 = K_wgsr(T7) * ((n7CO / n7) * (n7H2O / n7)) - (((n7CO2 / n7) * (n7H2 / n7)))
    eq32 = n7H2O - n6H2O + ghr_ksi1 + ghr_ksi2
    eq33 = n7H2 - n6H2 - 3 * ghr_ksi1 - ghr_ksi2
    eq34 = n7CO - n6CO - ghr_ksi1 + ghr_ksi2
    eq35 = np.dot(output[24:29], enthalpy(const, T6)) - np.dot(output[31:36], enthalpy(const, T7)) + Q3
    eq36 = P7 - P6

    # BLOCK 6: ATR

    atr_ksi1 = n7CH4 - n8CH4 - nO2 / 2
    atr_ksi2 = n8CO2 - n7CO2 - nO2 / 2

    T8 = 1323

    eq37 = K_smr(T8) * ((n8CH4 / n8) * (n8H2O / n8)) - ((n8CO / n8) * (n8H2 / n8) ** 3)
    eq38 = K_wgsr(T8) * ((n8CO / n8) * (n8H2O / n8)) - ((n8CO2 / n8) * (n8H2 / n8))
    eq39 = n8H2O - n7H2O + atr_ksi1 + atr_ksi2 - nO2
    eq40 = n8H2 - n7H2 - 3 * atr_ksi1 - atr_ksi2
    eq41 = n8CO - n7CO - atr_ksi1 + atr_ksi2
    eq42 = np.dot(output[31:36], enthalpy(const, T7)) - np.dot(output[38:43], enthalpy(const, T8)) + nO2 * 23.7

    eq43 = P8 - P7

    #eq37, eq38, eq39, eq40, eq41, eq42, eq43 = ATR_Reactor(np.array([n7CH4,n7H2O,n7H2, n7CO, n7CO2, P7, T7]),
    #np.array([nO2,n8CH4, n8H2O, n8H2, n8CO, n8CO2, P8]))

    # BLOCK 7: Post-ATR heat exchanger for making the model make sense by setting the
    # post temperature after ATR as the duty required for heating up GHR
    eq44 = n9CH4 - n8CH4
    eq45 = n9H2O - n8H2O
    eq46 = n9H2 - n8H2
    eq47 = n9CO - n8CO
    eq48 = n9CO2 - n8CO2
    eq49 = np.dot(output[38:43], enthalpy(const, T8)) - \
           np.dot(output[45:50], enthalpy(const, T9)) + Q3
    eq50 = P9 - P8

    # BLOCK 8: Isothermal temperature shift reactor
    T10 = 523  # Setting the reactor temperature and solving for the heat required instead
    itsr_ksi1 = n9H2O - n10H2O

    eq51 = n10CH4 - n9CH4
    eq52 = n10H2 - n9H2 - itsr_ksi1
    eq53 = n10CO - n9CO + itsr_ksi1
    eq54 = n10CO2 - n9CO2 - itsr_ksi1
    eq55 = K_wgsr(T10) * ((n10CO / n10) * (n10H2O / n10)) - ((n10CO2 / n10) * (n10H2 / n10))
    eq56 = np.dot(output[45:50], enthalpy(const, T9)) - \
           np.dot(output[52:57], enthalpy(const, T10)) + Q4
    eq57 = P10 - P9

    # BLOCK 9: Precondensate heat exchanger, same as ITSR, the temperature is set and we solve for heat gained
    T11 = 313
    eq58 = n11CH4 - n10CH4
    eq59 = n11H2O - n10H2O
    eq60 = n11H2 - n10H2
    eq61 = n11CO - n10CO
    eq62 = n11CO2 - n10CO2
    eq63 = np.dot(output[52:57], enthalpy(const, T10)) - \
           np.dot(output[59:64], enthalpy(const, T11)) + Q5
    eq64 = P11 - P10

    # BLOCK 10: Condensate
    eq65 = n12CH4 - f1[0] * n11CH4
    eq66 = n12H2O - f1[1] * n11H2O
    eq67 = n12H2 - f1[2] * n11H2
    eq68 = n12CO - f1[3] * n11CO
    eq69 = n12CO2 - f1[4] * n11CO2
    eq70 = P12 - P11
    eq71 = T12 - T11

    eq72 = P13 - P11
    eq73 = T13 - T11

    # BLOCK 11 PSA
    eq74 = n14CH4 - f2[0] * n12CH4
    eq75 = n14H2O - f2[1] * n12H2O
    eq76 = n14H2 - f2[2] * n12H2
    eq77 = n14CO - f2[3] * n12CO
    eq78 = n14CO2 - f2[4] * n12CO2
    eq79 = P14 - P12
    eq80 = T14 - T12

    eq81 = P15 - P12
    eq82 = T15 - T12

    system = np.array([eq0, eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14, eq15,
                       eq16, eq17, eq18, eq19, eq20, eq21, eq22, eq23, eq24, eq25, eq26, eq27, eq28, eq29,
                       eq30, eq31, eq32, eq33, eq34, eq35, eq36, eq37, eq38, eq39, eq40, eq41, eq42, eq43,
                       eq44, eq45, eq46, eq47, eq48, eq49, eq50, eq51, eq52, eq53, eq54, eq55, eq56, eq57,
                       eq58, eq59, eq60, eq61, eq62, eq63, eq64, eq65, eq66, eq67, eq68, eq69, eq70, eq71,
                       eq72, eq73, eq74, eq75, eq76, eq77, eq78, eq79, eq80, eq81, eq82])
    return system


def function(guess):
    return system(inletstream, guess, splitratio)


def solve(f, guess, grad_function):
    sol = opt.fsolve(f, guess)
    # print(sol)
    return sol


grad_function = jacobian(function)

# sol = opt.fsolve(function,guess,fprime=grad_function)


# print("error:",function(sol))
sol = solve(function, guess, grad_function)


def error_output():
    sol = opt.fsolve(function, guess, fprime=grad_function, xtol=1e-8)
    # print(sol)
    errors = function(sol)
    print(sum(abs(errors)))
    for i in range(len(errors)):
        print("eq", i, ":", errors[i])


def summary():
    print('{:<12}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}'
          .format("Stream", "CH4", "H2O", "H2", "CO", "CO2", "P", "T", "C2H6", "C3H8", "n-C4H10", "i-C4H10", "C5H12"))
    print('{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}'
          .format("", "[mol/s]", "[mol/s]", "[mol/s]", "[mol/s]", "[mol/s]", "[bar]", "[K]", "[mol/s]", "[mol/s]",
                  "[mol/s]", "[mol/s]", "[mol/s]"))
    print("------------------------------------------------------------------------------------------------"
          "--------------------------------------")

    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 1:", inletstream[0], inletstream[1], inletstream[2], inletstream[3], inletstream[4],
                    inletstream[10], inletstream[11], inletstream[5], inletstream[6], inletstream[7], inletstream[8],
                    inletstream[9]))

    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 2:", 0.00, sol[0], 0.00, 0.00, 0.00, sol[1], 423.00, 0.00, 0.00, 0.00, 0.00, 0.00))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 3:", sol[3], sol[4], sol[5], sol[6], sol[7], sol[8], sol[9], inletstream[5], inletstream[6],
                    inletstream[7], inletstream[8], inletstream[9]))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 4:", sol[10], sol[11], sol[12], sol[13], sol[14], sol[15], 693.00, inletstream[5],
                    inletstream[6], inletstream[7], inletstream[8], inletstream[9]))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 5:", sol[17], sol[18], sol[19], sol[20], sol[21], sol[22], sol[23], 0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 6:", sol[24], sol[25], sol[26], sol[27], sol[28], sol[29], 753, 0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 7:", sol[31], sol[32], sol[33], sol[34], sol[35], sol[36], 973, 0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 8:", sol[38], sol[39], sol[40], sol[41], sol[42], sol[43], 1323, 0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 9:", sol[45], sol[46], sol[47], sol[48], sol[49], sol[50], sol[51], 0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 10:", sol[52], sol[53], sol[54], sol[55], sol[56], sol[57], 523, 0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 11:", sol[59], sol[60], sol[61], sol[62], sol[63], sol[64], 313, 0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 12:", sol[65], sol[66], sol[67], sol[68], sol[69], sol[70], sol[71], 0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 13:", sol[59] - sol[65], sol[60] - sol[66], sol[61] - sol[67], sol[62] - sol[68],
                    sol[63] - sol[69], sol[72], sol[73], 0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 14:", sol[74], sol[75], sol[76], sol[77], sol[78], sol[79], sol[80], 0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 15:", sol[65] - sol[74], sol[66] - sol[75], sol[67] - sol[76], sol[68] - sol[77],
                    sol[69] - sol[78], sol[79], sol[80], 0, 0, 0, 0, 0))

    print('{:<10}{:.2f}'.format("Q1 [kJ]: ", sol[2]))
    print('{:<10}{:.2f}'.format("Q2 [kJ]: ", sol[16]))
    print('{:<10}{:.2f}'.format("Q3 [kJ]: ", sol[30]))
    print('{:<10}{:.2f}'.format("Q4 [kJ]: ", sol[37]))
    print('{:<10}{:.2f}'.format("Q5 [kJ]: ", sol[58]))
    print('{:<10}{:.2f}'.format("nO2 [mol/s]: ", sol[44]))


def mass():
    print('{:<15}{:.4f}'.format("Stream 1 [g/s]: ", Mm_heavy @ inletstream[:10]))
    print('{:<15}{:.4f}'.format("Stream 2 [g/s]: ", Mm @ [0, sol[0], 0, 0, 0]))
    print('{:<15}{:.4f}'.format("Stream 3 [g/s]: ", Mm @ sol[3:8] + (Mm_heavy[5:] @ inletstream[5:10])))
    print('{:<15}{:.4f}'.format("Stream 4 [g/s]: ", Mm @ sol[10:15] + (Mm_heavy[5:] @ inletstream[5:10])))
    print('{:<15}{:.4f}'.format("Stream 5 [g/s]: ", Mm @ sol[17:22]))
    print('{:<15}{:.4f}'.format("Stream 6 [g/s]: ", Mm @ sol[24:29]))
    print('{:<15}{:.4f}'.format("Stream 7 [g/s]: ", Mm @ sol[31:36]))
    print('{:<15}{:.4f}'.format("Stream O2 [g/s]: ", sol[44] * 32))
    print('{:<15}{:.4f}'.format("Stream 8 [g/s]: ", Mm @ sol[38:43]))
    print('{:<15}{:.4f}'.format("Stream 9 [g/s]: ", Mm @ sol[45:50]))
    print('{:<15}{:.4f}'.format("Stream 10 [g/s]: ", Mm @ sol[52:57]))
    print('{:<15}{:.4f}'.format("Stream 11 [g/s]: ", Mm @ sol[59:64]))
    print('{:<15}{:.4f}'.format("Stream 12 [g/s]: ", Mm @ sol[65:70]))
    print('{:<15}{:.4f}'.format("Stream 13 [g/s]: ",
                                Mm @ [sol[59] - sol[65], sol[60] - sol[66], sol[61] - sol[67], sol[62] - sol[68],
                                      sol[63] - sol[69]]))
    print('{:<15}{:.4f}'.format("Stream 14 [g/s]: ", Mm @ sol[74:79]))
    print('{:<15}{:.4f}'.format("Stream 15 [g/s]: ",
                                Mm @ [sol[65] - sol[74], sol[66] - sol[75], sol[67] - sol[76], sol[68] - sol[77],
                                      sol[69] - sol[78]]))


summary()
error_output()
mass()
