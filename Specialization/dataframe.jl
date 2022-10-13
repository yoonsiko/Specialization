#using Pkg
#Pkg.add("DataFrames")
#Pkg.add("PrettyTables")

function printTable(model)
    # Stream table
    CH₄ = [value(model[:mix_in_mol][1]), 0.0, value(model[:prePR_in_mol][1]), value(model[:pr_in_mol][1]),
     value(model[:preGHR_in_mol][1]), value(model[:ghr_in_mol][1]), value(model[:atr_in_mol][1]),
     value(model[:postATR_in_mol][1]), value(model[:itsr_in_mol][1]), value(model[:preCond_in_mol][1]),
     value(model[:cond_in_mol][1]), value(model[:cond_outPurge_mol][1]), value(model[:psa_in_mol][1]),
     value(model[:psa_outProduct_mol][1]), value(model[:psa_outPurge_mol][1])];

    H₂O = [value(model[:mix_in_mol][2]), value(model[:H2Ostream]), value(model[:prePR_in_mol][2]),
    value(model[:pr_in_mol][2]), value(model[:preGHR_in_mol][2]), value(model[:ghr_in_mol][2]),
    value(model[:atr_in_mol][2]),value(model[:postATR_in_mol][2]), value(model[:itsr_in_mol][2]),
    value(model[:preCond_in_mol][2]), value(model[:cond_in_mol][2]), value(model[:cond_outPurge_mol][2]),
    value(model[:psa_in_mol][2]), value(model[:psa_outProduct_mol][2]), value(model[:psa_outPurge_mol][2])];

    H₂ = [value(model[:mix_in_mol][3]), 0.0, value(model[:prePR_in_mol][3]), value(model[:pr_in_mol][3]),
    value(model[:preGHR_in_mol][3]), value(model[:ghr_in_mol][3]), value(model[:atr_in_mol][3]),
    value(model[:postATR_in_mol][3]), value(model[:itsr_in_mol][3]), value(model[:preCond_in_mol][3]),
    value(model[:cond_in_mol][3]), value(model[:cond_outPurge_mol][3]), value(model[:psa_in_mol][3]),
    value(model[:psa_outProduct_mol][3]), value(model[:psa_outPurge_mol][3])];

    CO = [value(model[:mix_in_mol][4]), 0.0, value(model[:prePR_in_mol][4]), value(model[:pr_in_mol][4]),
    value(model[:preGHR_in_mol][4]), value(model[:ghr_in_mol][4]), value(model[:atr_in_mol][4]),
    value(model[:postATR_in_mol][4]), value(model[:itsr_in_mol][4]), value(model[:preCond_in_mol][4]),
     value(model[:cond_in_mol][4]), value(model[:cond_outPurge_mol][4]), value(model[:psa_in_mol][4]),
    value(model[:psa_outProduct_mol][4]), value(model[:psa_outPurge_mol][4])];

    CO₂ = [value(model[:mix_in_mol][5]), 0.0, value(model[:prePR_in_mol][5]), value(model[:pr_in_mol][5]),
    value(model[:preGHR_in_mol][5]), value(model[:ghr_in_mol][5]), value(model[:atr_in_mol][5]),
    value(model[:postATR_in_mol][5]), value(model[:itsr_in_mol][5]), value(model[:preCond_in_mol][5]),
     value(model[:cond_in_mol][5]), value(model[:cond_outPurge_mol][5]), value(model[:psa_in_mol][5]),
    value(model[:psa_outProduct_mol][5]), value(model[:psa_outPurge_mol][5])];

    C₂H₆ = [value(model[:mix_in_mol][6]), 0.0, value(model[:prePR_in_mol][6]), value(model[:pr_in_mol][6]),
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    C₃H₈ = [value(model[:mix_in_mol][7]), 0.0, value(model[:prePR_in_mol][7]), value(model[:pr_in_mol][7]),
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    nC₄H₁₀ = [value(model[:mix_in_mol][8]), 0.0, value(model[:prePR_in_mol][8]), value(model[:pr_in_mol][8]),
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    iC₄H₁₀ = [value(model[:mix_in_mol][9]), 0.0, value(model[:prePR_in_mol][9]), value(model[:pr_in_mol][9]),
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    C₅₊ = [value(model[:mix_in_mol][10]), 0.0, value(model[:prePR_in_mol][10]), value(model[:pr_in_mol][10]),
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    T = [value(model[:mix_in_T]), value(model[:H2O_T]), value(model[:prePR_in_T]), value(model[:pr_in_T]),
    value(model[:preGHR_in_T]), value(model[:ghr_in_T]), value(model[:atr_in_T]), value(model[:postATR_in_T]),
    value(model[:itsr_in_T]), value(model[:preCond_in_T]), value(model[:cond_in_T]), value(model[:cond_outPurge_T]),
    value(model[:psa_in_T]), value(model[:psa_outProduct_T]), value(model[:psa_outPurge_T])];

    streamdf = DataFrame(T = T, CH₄ = CH₄, H₂O = H₂O, H₂ = H₂, CO = CO, CO₂ = CO₂,
     C₂H₆ = C₂H₆, C₃H₈ = C₃H₈, nC₄H₁₀ = nC₄H₁₀, iC₄H₁₀ = iC₄H₁₀, C₅₊ = C₅₊);

    # Other variable table
    variable = ["prePR_Q","preGHR_Q","ghr_Q", "postATR_Q", "itsr_Q","preCond_Q","nO2","additional_Q"]
    values = [value(model[:prePR_Q]),value(model[:preGHR_Q]),value(model[:ghr_Q]), value(model[:postATR_Q]),
    value(model[:itsr_Q]),value(model[:preCond_Q]),value(model[:nO2]),value(model[:additional_Q])];
    otherdf = DataFrame(Variable = variable, Value = values);

    # Mass table
    C = 12.01;
    H = 1.008;
    O = 16;
    # CH4, H2O, H2, CO, CO2, C2H6, C3H8, n-C4H10, i-C4H10, C5+
    Mm = [C+H*4, H*2+O, H*2, C+O, C+O*2, C*2+H*6, C*3+H*8, C*4+H*10, C*4+H*10, C*5+H*12];
    mass = zeros(15);
    for i=1:10
        mass[1] += value(model[:mix_in_mol][i]) * Mm[i];
        mass[3] += value(model[:prePR_in_mol][i]) * Mm[i];
        mass[4] += value(model[:pr_in_mol][i]) * Mm[i];
    end

    mass[2] = value(model[:H2Ostream])*Mm[2];

    for j=1:5
        mass[5] += value(model[:preGHR_in_mol][j])*Mm[j]
        mass[6] += value(model[:ghr_in_mol][j])*Mm[j]
        mass[7] += value(model[:atr_in_mol][j])*Mm[j]
        mass[8] += value(model[:postATR_in_mol][j])*Mm[j]
        mass[9] += value(model[:itsr_in_mol][j])*Mm[j]
        mass[10] += value(model[:preCond_in_mol][j])*Mm[j]
        mass[11] += value(model[:cond_in_mol][j])*Mm[j]
        mass[12] += value(model[:psa_in_mol][j])*Mm[j]
        mass[13] += value(model[:cond_outPurge_mol][j])*Mm[j]
        mass[14] += value(model[:psa_outProduct_mol][j])*Mm[j]
        mass[15] += value(model[:psa_outPurge_mol][j])*Mm[j]
    end
    oxygenstream = zeros(15);
    oxygenstream[7] = value(model[:nO2])*32;
    massdf = DataFrame(Mass= mass, O₂ = oxygenstream);

    return streamdf, otherdf, massdf
end