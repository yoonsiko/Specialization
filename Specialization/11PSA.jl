#import Pkg
#Pkg.add("JuMP")
#Pkg.add("Ipopt")
#Pkg.add("MathOptInterface")
"""
using JuMP, Ipopt, MathOptInterface
include("enthalpy.jl")
const MOI = MathOptInterface


Base.@kwdef mutable struct psa_par
    outPurge_T::Float64 = 313.0;
    outProduct_T::Float64 = 313.0;
    splitratio::Vector = [0.01, 0.01, 0.999, 0.01, 0.01]; # Almost ideal
end

Base.@kwdef mutable struct _par
  psa::psa_par=psa_par();
  hconst = heavy_const;
end

par = _par()
#"""

function PSA_model(model, par)
  @variable(model, 0 <= psa_in_mol[1:5]); # Stream 12
  @variable(model, 0 <= psa_outPurge_mol[1:5]); # Stream 15 (Purge)
  @variable(model, 0 <= psa_outProduct_mol[1:5]); # Stream 14 (H2)

  ini_psa_in = [0.00032515179570229987, 0.05589429423256371, 656.2769122396338, 23.343132770546816, 188.85192754301332];  # CH4 H2O H2 CO CO2
  ini_psa_outPurge = [0.0003219002777452769, 0.05533535129023807, 0.6562769122396344, 23.10970144284135, 186.9634082675832];
  ini_psa_outProduct = [3.2515179570229986e-06, 0.0005589429423256371, 655.6206353273942, 0.23343132770546818, 1.8885192754301332];

  for i=1:5
    set_start_value(psa_in_mol[i] , ini_psa_in[i]);
    set_start_value(psa_outPurge_mol[i], ini_psa_outPurge[i]);
    set_start_value(psa_outProduct_mol[i], ini_psa_outProduct[i]);
  end

  @variable(model, 273 <= psa_in_T, start = 313.00);
  @variable(model, 273 <= psa_outPurge_T, start = 313.00);
  @variable(model, 273 <= psa_outProduct_T, start = 313.00)

  psa_H_outPurge = build_enthalpy(model, psa_outPurge_T, par)
  psa_H_outProduct = build_enthalpy(model, psa_outProduct_T, par)
  psa_H_in = build_enthalpy(model, psa_in_T, par)
  
  # Constraints
  # Mass balance
  @NLconstraint(model, par.psa.splitratio[1]*psa_in_mol[1] - psa_outProduct_mol[1] == 0);
  @NLconstraint(model, par.psa.splitratio[2]*psa_in_mol[2] - psa_outProduct_mol[2] == 0);
  @NLconstraint(model, par.psa.splitratio[3]*psa_in_mol[3] - psa_outProduct_mol[3] == 0);
  @NLconstraint(model, par.psa.splitratio[4]*psa_in_mol[4] - psa_outProduct_mol[4] == 0);
  @NLconstraint(model, par.psa.splitratio[5]*psa_in_mol[5] - psa_outProduct_mol[5] == 0);
  @NLconstraint(model, (1-par.psa.splitratio[1])*psa_in_mol[1] - psa_outPurge_mol[1] == 0);
  @NLconstraint(model, (1-par.psa.splitratio[2])*psa_in_mol[2] - psa_outPurge_mol[2] == 0);
  @NLconstraint(model, (1-par.psa.splitratio[3])*psa_in_mol[3] - psa_outPurge_mol[3] == 0);
  @NLconstraint(model, (1-par.psa.splitratio[4])*psa_in_mol[4] - psa_outPurge_mol[4] == 0);
  @NLconstraint(model, (1-par.psa.splitratio[5])*psa_in_mol[5] - psa_outPurge_mol[5] == 0);

  #@NLconstraint(model, sum(psa_H_outPurge[i]*psa_outPurge_mol[i] - psa_H_in[i]*(1-par.psa.splitratio[i])*psa_in_mol[i] for i=1:5) == 0);
  #@NLconstraint(model, sum(psa_H_outProduct[i]*psa_outProduct_mol[i] - psa_H_in[i]*par.psa.splitratio[i]*psa_in_mol[i] for i=1:5) == 0);

  @NLconstraint(model, psa_outPurge_T - par.psa.outPurge_T == 0);
  @NLconstraint(model, psa_outProduct_T - par.psa.outProduct_T == 0);


  return model;

end

"""
m = Model(Ipopt.Optimizer);
m = PSA_model(m, par);
#optimize!(m);
d = JuMP.NLPEvaluator(m);
MOI.initialize(d, [:Jac]);
constraint_values = zeros(1,14);
inp = [0.00032515179570229987, 0.05589429423256371, 656.2769122396338, 23.343132770546816, 188.85192754301332,
      0.0003219002777452769, 0.05533535129023807, 0.6562769122396344, 23.10970144284135, 186.9634082675832,
      3.2515179570229986e-06, 0.0005589429423256371, 655.6206353273942, 0.23343132770546818, 1.8885192754301332, 
      313.00, 313.00, 313.00];
MOI.eval_constraint(d, constraint_values, inp[:]);
@show constraint_values
#"""