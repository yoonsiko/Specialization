#using Pkg
#Pkg.add("Plots")
#Pkg.add("PyPlot")
using JuMP, Ipopt, MathOptInterface, DataFrames, PrettyTables, Plots
include("enthalpy.jl")
include("0par.jl")
include("1MIX.jl")
include("2PrePR.jl")
include("3PR.jl")
include("4PreGHR.jl")
include("5GHR.jl")
include("6ATR.jl")
include("7PostATR.jl")
include("8ITSR.jl")
include("9PreCondensate.jl")
include("10Condensate.jl")
include("11PSA.jl")
include("dataframe.jl")
include("equilibrium.jl")
include("compFunc.jl")

const MOI = MathOptInterface

Base.@kwdef mutable struct _par
    mix::mix_par=mix_par();
    prePR::prePR_par=prePR_par();
    pr::pr_par = pr_par();
    preGHR::preGHR_par = preGHR_par();
    ghr::ghr_par = ghr_par();
    atr::atr_par = atr_par();
    postATR::postATR_par = postATR_par();
    itsr::itsr_par = itsr_par();
    preCond::preCond_par = preCond_par();
    cond::cond_par = cond_par();
    psa::psa_par = psa_par();
    hconst = heavy_const;
    smr_const = smr_const;
    wgsr_const = wgsr_const;
end


global_x_CH4_init = [];
global_x_H2_final = [];
global_n_CH4_init = [];
global_n_H2_final = [];

initial_value = [113.76022309722791, 0.0, 0.0, 0.0, 1.9483473792214394, 8.869342547202073,
9.741736896107197, 3.605896642141171, 2.0501267199270368, 5.379765151581587];
total_mole = sum(initial_value);

init = [0.6726347904371311, 0.0, 0.0, 0.0, 0.020187270170927026, 0.09189727465869764,
     0.10093635085463512, 0.03736151494320822, 0.02124182906045306, 0.05574096987494777]*total_mole;

for i = 1:30
    initcomp = comp(init, i*0.01);

    init_val = initcomp*total_mole

    par = _par();
    par.mix.ini_mix_in = init_val;
    optimizer = optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "constr_viol_tol" => 1e-8)
    m = Model(optimizer);

    MIX_model(m, par);
    prePR_model(m, par);
    PR_model(m, par);
    preGHR_model(m, par);
    GHR_model(m, par);
    ATR_model(m, par);
    postATR_model(m, par);
    ITSR_model(m, par);
    preCond_model(m, par);
    Cond_model(m, par);
    PSA_model(m, par);

    for i = 1:10
        @NLconstraint(m, m[:mix_out_mol][i] - m[:prePR_in_mol][i] == 0)
        @NLconstraint(m, m[:prePR_out_mol][i] - m[:pr_in_mol][i] == 0)
    end
    #@NLconstraint(m, m[:mix_out_mol][i] == 0)
    for i = 1:5 # After all heavier carbons are removed
        @NLconstraint(m, m[:pr_out_mol][i] - m[:preGHR_in_mol][i] == 0)
        @NLconstraint(m, m[:preGHR_out_mol][i] - m[:ghr_in_mol][i] == 0)
        @NLconstraint(m, m[:ghr_out_mol][i] - m[:atr_in_mol][i] == 0)
        @NLconstraint(m, m[:atr_out_mol][i] - m[:postATR_in_mol][i] == 0)
        @NLconstraint(m, m[:postATR_out_mol][i] - m[:itsr_in_mol][i] == 0)
        @NLconstraint(m, m[:itsr_out_mol][i] - m[:preCond_in_mol][i] == 0)
        @NLconstraint(m, m[:preCond_out_mol][i] - m[:cond_in_mol][i] == 0)
        @NLconstraint(m, m[:cond_outProduct_mol][i] - m[:psa_in_mol][i] == 0)
    end
    # Same for the temperature
    @NLconstraint(m, m[:mix_out_T] - m[:prePR_in_T] == 0);
    @NLconstraint(m, m[:prePR_out_T] - m[:pr_in_T] == 0);
    @NLconstraint(m, m[:pr_out_T] - m[:preGHR_in_T] == 0);
    @NLconstraint(m, m[:preGHR_out_T] - m[:ghr_in_T] == 0);
    @NLconstraint(m, m[:ghr_out_T] - m[:atr_in_T] == 0);
    @NLconstraint(m, m[:atr_out_T] - m[:postATR_in_T] == 0);
    @NLconstraint(m, m[:postATR_out_T] - m[:itsr_in_T] == 0);
    @NLconstraint(m, m[:itsr_out_T] - m[:preCond_in_T] == 0);
    @NLconstraint(m, m[:preCond_out_T] - m[:cond_in_T] == 0);
    @NLconstraint(m, m[:cond_outProduct_T] - m[:psa_in_T] == 0);


    for i = 1:10
        @NLconstraint(m, par.mix.ini_mix_in[i] - m[:mix_in_mol][i] == 0);
    end
    @NLconstraint(m, par.mix.in_T - m[:mix_in_T] == 0);
    @NLconstraint(m, par.mix.H2O_T - m[:H2O_T] == 0);


    @NLconstraint(m, m[:ghr_Q] + m[:postATR_Q] == 0); #- additional_Q == 0);
    @NLconstraint(m, m[:atr_out_T] - m[:ghr_out_T] >= 25);
    @NLconstraint(m, m[:postATR_out_T] - m[:ghr_in_T] >= 25);
    @NLobjective(m, Max, m[:psa_outProduct_mol][3]);
    optimize!(m)
    
    #total_final_mole_stream = sum(value(m[:psa_outProduct_mol][i]) for i = 1:5);
    append!(global_x_CH4_init, initcomp[1]); 
    #append!(global_x_H2_final, m[:psa_outProduct_mol][3]/total_final_mole_stream); 
    append!(global_n_CH4_init, initcomp[1]*total_mole); 
    append!(global_n_H2_final, value(m[:psa_outProduct_mol][3])); 
end 

#Plots.pyplot()
display(plot(global_n_CH4_init, global_n_H2_final,
            title = "Initial CH4 vs final H2",
            xlabel = "CH4 [kmol/h]",
            ylabel = "H2 [kmol/h]",
            label = "Final H2"));