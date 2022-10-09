# Author: Yoonsik Oh
# Project: SUBPRO Blue hydrogen

from enthalpy import const, enthalpy, enthalpy_heavy, heavy_const
from Keq import K_smr, K_wgsr
import scipy.optimize as opt
import autograd.numpy as np
from autograd import grad, jacobian
from atr import ATR_Reactor
from mixing import mixingTemp
from prereformer import PR_Reactor
from ghr import GHR_Reactor
from postATR import postATR_temp
from itsr import itsr


# inlet/outlet = [molar flow, T, P]
C = 12
H = 1.
O = 16
# molar flow component order = CH4, H2O, H2, CO, CO2
Mm = np.array([C * 1 + H * 4, H * 2 + O, H * 2, C + O, C + O * 2])  # This is to check the conservation of mass
# molar flow component order = CH4, H2O, H2, CO, CO2  C2H6      C3H8      n-C4H10   i-C4H10   C5H12


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

splitratio = np.array([[1.0, 0.001, 1.0, 1.0, 1.0], [0.01, 0.01, 0.999, 0.01, 0.01]])


# Input is given, Q and output is the initial guess values
def system(input, f):
    # ----------------------------------------------------------------------------------
    # Defining and unpacking all the variables
    # ----------------------------------------------------------------------------------
    n1CH4, n1H2O, n1H2, n1CO, n1CO2, \
    n1C2H6, n1C3H8, n1i_C4H10, n1n_C4H10, n1_C5, P1, T1 = input
    f1, f2 = f

    # ----------------------------------------------------------------------------------
    # Defining variables
    # ----------------------------------------------------------------------------------

    nCarbon = n1CH4 + n1C2H6 + n1C3H8 + n1n_C4H10 + n1i_C4H10 + n1_C5
    stcr = 2.5  # steam to carbon ratio (molar ratio)
    # ----------------------------------------------------------------------------------
    # BLOCK 1: mixing purified natural gas and steam and setting S/C to 2.5
    # ----------------------------------------------------------------------------------

    n2H2O = stcr * nCarbon
    n3CH4 = n1CH4
    n3H2O = n2H2O + n1H2O
    n3H2 = n1H2
    n3CO = n1CO
    n3CO2 = n1CO2
    P2 = P1
    P3 = P2

    T2 = 423
    inlet1 = input
    inlet1 = np.append(inlet1, np.array([n2H2O, T2]))

    def T3function(guess):
        return mixingTemp(inlet1, guess)

    sol1 = opt.fsolve(T3function, np.array([T2]))
    print("T3", T3function([381.41]))

    T3 = sol1[0]

    # ----------------------------------------------------------------------------------
    # BLOCK 2: Pre-prereformer heat exchanger
    # ----------------------------------------------------------------------------------
    T4 = 693
    n4CH4 = n3CH4
    n4H2O = n3H2O
    n4H2 = n3H2
    n4CO = n3CO
    n4CO2 = n3CO2

    n4stream = np.array([n4CH4, n4H2O, n4H2, n4CO, n4CO2])
    n4stream = np.append(n4stream, input[5:10])
    n3stream = n4stream
    Q1 = np.dot(n4stream, enthalpy_heavy(heavy_const, T4)) - np.dot(n3stream, enthalpy_heavy(heavy_const, T3))
    P4 = P3

    # ----------------------------------------------------------------------------------
    # BLOCK 3: Pre-reformer (Removing all heavier carbons to methane, (full conversion))
    # ----------------------------------------------------------------------------------
    stream4 = np.array([n4CH4, n4H2O, n4H2, n4CO, n4CO2])
    stream4 = np.append(stream4, input[5:10])
    stream4 = np.append(stream4, T4)
    guess3 = np.array([100, 500, 100, 1, 100, 700])

    def PR_function(guess):
        return PR_Reactor(stream4, guess)

    sol3 = opt.fsolve(PR_function, guess3)

    n5CH4, n5H2O, n5H2, n5CO, n5CO2, T5 = sol3
    stream5 = np.array([n5CH4, n5H2O, n5H2, n5CO, n5CO2])
    P5 = P4

    # ----------------------------------------------------------------------------------
    # BLOCK 4: Pre-GHR Heat Exchanger
    # ----------------------------------------------------------------------------------
    T6 = 973
    n6CH4 = n5CH4
    n6H2O = n5H2O
    n6H2 = n5H2
    n6CO = n5CO
    n6CO2 = n5CO2

    stream6 = np.array([n6CH4, n6H2O, n6H2, n6CO, n6CO2])
    Q2 = np.dot(stream6, enthalpy(const, T6)) - np.dot(stream5, enthalpy(const, T5))
    P6 = P5

    # ----------------------------------------------------------------------------------
    # BLOCK 5: GHR
    # ----------------------------------------------------------------------------------

    guess5 = np.array([100, 500, 500, 100, 100, 28000.0])
    T7 = 1073
    input5 = np.append(stream6, np.array([T6]))

    def GHR_function(guess):
        return GHR_Reactor(input5, guess)


    sol5 = opt.fsolve(GHR_function, guess5)


    output_for_error = [22.543406840790453, 122.15875765109533, 544.9262854431126, 139.04829250631568, 50.60368611824971, 34472.61];
    print("GHR function error:", GHR_function(output_for_error))
    n7CH4, n7H2O, n7H2, n7CO, n7CO2, Q3 = sol5
    P7 = P6

    # ----------------------------------------------------------------------------------
    # BLOCK 6: ATR
    # ----------------------------------------------------------------------------------
    stream7 = np.array([n7CH4, n7H2O, n7H2, n7CO, n7CO2])
    guess6 = np.array([100,100, 500, 500, 100, 100, 30])
    inlet6 = np.append(stream7, np.array([30, T7]))
    T8 = 1600

    def ATR_function(guess):
        return ATR_Reactor(inlet6, guess)

    sol6 = opt.fsolve(ATR_function, guess6)

    nO2, n8CH4, n8H2O, n8H2, n8CO, n8CO2, P8 = sol6
    stream8 = np.array([n8CH4, n8H2O, n8H2, n8CO, n8CO2])

    # ----------------------------------------------------------------------------------
    # BLOCK 7: Post-ATR heat exchanger for making the model make sense by setting the
    # post temperature after ATR as the duty required for heating up GHR
    # ----------------------------------------------------------------------------------
    n9CH4 = n8CH4
    n9H2O = n8H2O
    n9H2 = n8H2
    n9CO = n8CO
    n9CO2 = n8CO2
    stream9 = np.array([n9CH4, n9H2O, n9H2, n9CO, n9CO2])
    inlet7 = np.append(stream8, stream9)
    inlet7 = np.append(inlet7, np.array([T8, -Q3]))


    def postAtr_function(guess):
        return postATR_temp(inlet7, guess)


    guess7 = np.array([400])

    sol7 = opt.fsolve(postAtr_function, guess7)


    T9 = sol7[0]

    P9 = P8

    # ----------------------------------------------------------------------------------
    # BLOCK 8: Isothermal temperature shift reactor
    # ----------------------------------------------------------------------------------
    T10 = 523  # Setting the reactor temperature and solving for the heat required instead

    input8 = np.append(stream9, T9)
    guess8 = np.array([1, 500, 500, 100, 100, -100000])

    def itsr_function(guess):
        return itsr(input8, guess)


    sol8 = opt.fsolve(itsr_function, guess8)

    n10CH4, n10H2O, n10H2, n10CO, n10CO2, Q4 = sol8

    P10 = P9

    # ----------------------------------------------------------------------------------
    # BLOCK 9: Precondensate heat exchanger, same as ITSR, the temperature is set
    # and we solve for heat gained
    # ----------------------------------------------------------------------------------
    T11 = 313
    n11CH4 = n10CH4
    n11H2O = n10H2O
    n11H2 = n10H2
    n11CO = n10CO
    n11CO2 = n10CO2
    stream10 = np.array([n10CH4, n10H2O, n10H2, n10CO, n10CO2])
    stream11 = np.array([n11CH4, n11H2O, n11H2, n11CO, n11CO2])
    Q5 = np.dot(stream11, enthalpy(const, T11)) - \
         np.dot(stream10, enthalpy(const, T10))

    P11 = P10

    # ----------------------------------------------------------------------------------
    # BLOCK 10: Condensate
    # ----------------------------------------------------------------------------------
    n12CH4 = f1[0] * n11CH4
    n12H2O = f1[1] * n11H2O
    n12H2 = f1[2] * n11H2
    n12CO = f1[3] * n11CO
    n12CO2 = f1[4] * n11CO2
    P12 = P11
    T12 = T11

    n13CH4 = (1 - f1[0]) * n11CH4
    n13H2O = (1 - f1[1]) * n11H2O
    n13H2 = (1 - f1[2]) * n11H2
    n13CO = (1 - f1[3]) * n11CO
    n13CO2 = (1 - f1[4]) * n11CO2

    P13 = P11
    T13 = T11

    # ----------------------------------------------------------------------------------
    # BLOCK 11 PSA
    # ----------------------------------------------------------------------------------
    n14CH4 = f2[0] * n12CH4
    n14H2O = f2[1] * n12H2O
    n14H2 = f2[2] * n12H2
    n14CO = f2[3] * n12CO
    n14CO2 = f2[4] * n12CO2
    P14 = P12
    T14 = T12

    n15CH4 = (1 - f2[0]) * n12CH4
    n15H2O = (1 - f2[1]) * n12H2O
    n15H2 = (1 - f2[2]) * n12H2
    n15CO = (1 - f2[3]) * n12CO
    n15CO2 = (1 - f2[4]) * n12CO2
    P15 = P12
    T15 = T12

    # ----------------------------------------------------------------------------------
    # Return the variables
    # ----------------------------------------------------------------------------------
    system = np.array([n1CH4, n1H2O, n1H2, n1CO, n1CO2, T1, P1,
                       0.0, n2H2O, 0.0, 0.0, 0.0, T2, P2,
                       n3CH4, n3H2O, n3H2, n3CO, n3CO2, T3, P3,
                       n4CH4, n4H2O, n4H2, n4CO, n4CO2, T4, P4,
                       n5CH4, n5H2O, n5H2, n5CO, n5CO2, T5, P5,
                       n6CH4, n6H2O, n6H2, n6CO, n6CO2, T6, P6,
                       n7CH4, n7H2O, n7H2, n7CO, n7CO2, T7, P7,
                       n8CH4, n8H2O, n8H2, n8CO, n8CO2, T8, P8,
                       n9CH4, n9H2O, n9H2, n9CO, n9CO2, T9, P9,
                       n10CH4, n10H2O, n10H2, n10CO, n10CO2, T10, P10,
                       n11CH4, n11H2O, n11H2, n11CO, n11CO2, T11, P11,
                       n12CH4, n12H2O, n12H2, n12CO, n12CO2, T12, P12,
                       n13CH4, n13H2O, n13H2, n13CO, n13CO2, T13, P13,
                       n14CH4, n14H2O, n14H2, n14CO, n14CO2, T14, P14,
                       n15CH4, n15H2O, n15H2, n15CO, n15CO2, T15, P15,
                       nO2, Q1, Q2, Q3, Q4, Q5])
    return system

data = system(inletstream,splitratio)

def summary(data):
    print("------------------------------------------------------------------------------------------------"
          "--------------------------------------")
    print('{:<12}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}'
          .format("Stream", "CH4", "H2O", "H2", "CO", "CO2", "T", "P", "C2H6", "C3H8", "n-C4H10", "i-C4H10", "C5H12"))
    print('{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}'
          .format("", "[mol/s]", "[mol/s]", "[mol/s]", "[mol/s]", "[mol/s]", "[K]", "[bar]", "[mol/s]", "[mol/s]",
                  "[mol/s]", "[mol/s]", "[mol/s]"))
    print("------------------------------------------------------------------------------------------------"
          "--------------------------------------")

    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 1:", data[0], data[1], data[2], data[3], data[4], data[5], data[6],
                    inletstream[5], inletstream[6], inletstream[7], inletstream[8], inletstream[9]))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 2:", data[7], data[8], data[9], data[10], data[11], data[12], data[13],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 3:", data[14], data[15], data[16], data[17], data[18], data[19], data[20],
                    inletstream[5], inletstream[6], inletstream[7], inletstream[8], inletstream[9]))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 4:", data[21], data[22], data[23], data[24], data[25], data[26], data[27],
                    inletstream[5], inletstream[6], inletstream[7], inletstream[8], inletstream[9]))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 5:", data[28], data[29], data[30], data[31], data[32], data[33], data[34],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 6:", data[35], data[36], data[37], data[38], data[39], data[40], data[41],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 7:", data[42], data[43], data[44], data[45], data[46], data[47], data[48],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 8:", data[49], data[50], data[51], data[52], data[53], data[54], data[55],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 9:", data[56], data[57], data[58], data[59], data[60], data[61], data[62],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 10:", data[63], data[64], data[65], data[66], data[67], data[68], data[69],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 11:", data[70], data[71], data[72], data[73], data[74], data[75], data[76],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 12:", data[77], data[78], data[79], data[80], data[81], data[82], data[83],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 13:", data[84], data[85], data[86], data[87], data[88], data[89], data[90],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 14:", data[91], data[92], data[93], data[94], data[95], data[96], data[97],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}'
            .format("Stream 15:", data[98], data[99], data[100], data[101], data[102], data[103], data[104],
                    0, 0, 0, 0, 0))
    print('{:<10}{:.2f}'.format("nO2 [mol/s]: ", data[105]))
    print('{:<10}{:.2f}'.format("Q1 [kJ]: ", data[106]))
    print('{:<10}{:.2f}'.format("Q2 [kJ]: ", data[107]))
    print('{:<10}{:.2f}'.format("Q3 [kJ]: ", data[108]))
    print('{:<10}{:.2f}'.format("Q4 [kJ]: ", data[109]))
    print('{:<10}{:.2f}'.format("Q5 [kJ]: ", data[110]))

#summary(data)

def mass(data):
    print('{:<15}{:.4f}'.format("Stream 1 [g/s]: ", Mm_heavy @ inletstream[:10]))
    print('{:<15}{:.4f}'.format("Stream 2 [g/s]: ", Mm @ data[7:12]))
    print('{:<15}{:.4f}'.format("Stream 3 [g/s]: ", Mm @ data[14:19] + (Mm_heavy[5:] @ inletstream[5:10])))
    print('{:<15}{:.4f}'.format("Stream 4 [g/s]: ", Mm @ data[21:26] + (Mm_heavy[5:] @ inletstream[5:10])))
    print('{:<15}{:.4f}'.format("Stream 5 [g/s]: ", Mm @ data[28:33]))
    print('{:<15}{:.4f}'.format("Stream 6 [g/s]: ", Mm @ data[35:40]))
    print('{:<15}{:.4f}'.format("Stream 7 [g/s]: ", Mm @ data[42:47]))
    print('{:<15}{:.4f}'.format("Stream O2 [g/s]: ", data[105] * 32))
    print('{:<15}{:.4f}'.format("Stream 8 [g/s]: ", Mm @ data[49:54]))
    print('{:<15}{:.4f}'.format("Stream 9 [g/s]: ", Mm @ data[56:61]))
    print('{:<15}{:.4f}'.format("Stream 10 [g/s]: ", Mm @ data[63:68]))
    print('{:<15}{:.4f}'.format("Stream 11 [g/s]: ", Mm @ data[70:75]))
    print('{:<15}{:.4f}'.format("Stream 12 [g/s]: ", Mm @ data[77:82]))
    print('{:<15}{:.4f}'.format("Stream 13 [g/s]: ", Mm @ data[84:89]))
    print('{:<15}{:.4f}'.format("Stream 14 [g/s]: ", Mm @ data[91:96]))
    print('{:<15}{:.4f}'.format("Stream 15 [g/s]: ", Mm @ data[98:103]))

def composition(data):
    n1 = np.sum(data[:5]+inletstream[5:10])
    n2 = np.sum(data[7:12])
    n3 = np.sum(data[14:19]+inletstream[5:10])
    n4 = np.sum(data[21:26]+inletstream[5:10])
    n5 = np.sum(data[28:33])
    n6 = np.sum(data[35:40])
    n7 = np.sum(data[42:47])
    n8 = np.sum(data[49:54])
    n9 = np.sum(data[56:61])
    n10 = np.sum(data[63:68])
    n11 = np.sum(data[70:75])
    n12 = np.sum(data[77:82])
    n13 = np.sum(data[84:89])
    n14 = np.sum(data[91:96])
    n15 = np.sum(data[98:103])

    print("------------------------------------------------------------------------------------------------"
          "--------------------------------------")
    print('{:<12}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}'
          .format("Stream", "xCH4", "xH2O", "xH2", "xCO", "xCO2", "T", "P", "xC2H6", "xC3H8", "xn-C4H10", "xni-C4H10", "xnC5H12"))
    print("------------------------------------------------------------------------------------------------"
          "--------------------------------------")

    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 1:", data[0]/n1, data[1]/n1, data[2]/n1, data[3]/n1, data[4]/n1, data[5], data[6],
                    inletstream[5]/n1, inletstream[6]/n1, inletstream[7]/n1, inletstream[8]/n1, inletstream[9]/n1))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 2:", data[7]/n2, data[8]/n2, data[9]/n2, data[10]/n2, data[11]/n2, data[12], data[13],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 3:", data[14]/n3, data[15]/n3, data[16]/n3, data[17]/n3, data[18]/n3, data[19], data[20],
                    inletstream[5]/n3, inletstream[6]/n3, inletstream[7]/n3, inletstream[8]/n3, inletstream[9]/n3))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 4:", data[21]/n4, data[22]/n4, data[23]/n4, data[24]/n4, data[25]/n4, data[26], data[27],
                    inletstream[5]/n4, inletstream[6]/n4, inletstream[7]/n4, inletstream[8]/n4, inletstream[9]/n4))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 5:", data[28]/n5, data[29]/n5, data[30]/n5, data[31]/n5, data[32]/n5, data[33], data[34],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 6:", data[35]/n6, data[36]/n6, data[37]/n6, data[38]/n6, data[39]/n6, data[40], data[41],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 7:", data[42]/n7, data[43]/n7, data[44]/n7, data[45]/n7, data[46]/n7, data[47], data[48],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 8:", data[49]/n8, data[50]/n8, data[51]/n8, data[52]/n8, data[53]/n8, data[54], data[55],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 9:", data[56]/n9, data[57]/n9, data[58]/n9, data[59]/n9, data[60]/n9, data[61], data[62],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 10:", data[63]/n10, data[64]/n10, data[65]/n10, data[66]/n10, data[67]/n10, data[68], data[69],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 11:", data[70]/n11, data[71]/n11, data[72]/n11, data[73]/n11, data[74]/n11, data[75], data[76],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 12:", data[77]/n12, data[78]/n12, data[79]/n12, data[80]/n12, data[81]/n12, data[82], data[83],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 13:", data[84]/n13, data[85]/n13, data[86]/n13, data[87]/n13, data[88]/n13, data[89], data[90],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 14:", data[91]/n14, data[92]/n14, data[93]/n14, data[94]/n14, data[95]/n14, data[96], data[97],
                    0, 0, 0, 0, 0))
    print(
        '{:<12}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.2f}{:<10.2f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}{:<10.4f}'
            .format("Stream 15:", data[98]/n15, data[99]/n15, data[100]/n15, data[101]/n15, data[102]/n15, data[103], data[104],
                    0, 0, 0, 0, 0))


summary(data)
mass(data)
composition(data)