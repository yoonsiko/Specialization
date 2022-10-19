using JuMP, Ipopt, MathOptInterface, DataFrames, PrettyTables
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

const MOI = MathOptInterface

Base.@kwdef mutable struct _par
    mix::mix_par=mix_par();
    prePR::prePR_par=prePR_par();
    pr::pr_par = pr_par();
    preGHR::preGHR_par = preGHR_par();
    ghr::ghr_par = ghr_par();
    atr::atr_par = atr_par(out_T = 1450);
    postATR::postATR_par = postATR_par();
    itsr::itsr_par = itsr_par();
    preCond::preCond_par = preCond_par();
    cond::cond_par = cond_par();
    psa::psa_par = psa_par();
    hconst = heavy_const;
    const_u = const_u;
    const_o = const_o;
end

par = _par();
m = Model(Ipopt.Optimizer);
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
#all_variables(m)
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


# Initial value conditions
for i = 1:10
    @NLconstraint(m, par.mix.ini_mix_in[i] - m[:mix_in_mol][i] == 0);
end
@NLconstraint(m, par.mix.in_T - m[:mix_in_T] == 0);
@NLconstraint(m, par.mix.H2O_T - m[:H2O_T] == 0);
# Flowsheet


# constraint for postATR_Q and ghr_Q
@variable(m, 0 <= additional_Q, start = 0);
@NLconstraint(m, m[:ghr_Q] + m[:postATR_Q] - additional_Q == 0);
@NLobjective(m, Min, additional_Q);
optimize!(m)
#@show m
streamdf, otherdf, massdf = printTable(m);
println("Stream table")
show(streamdf, allrows=true);
println("\n\nOther variables")
show(otherdf, allrows=true);
println("\n\nMass table")
show(massdf, allrows=true);