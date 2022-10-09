#import Pkg
#Pkg.add("JuMP")
#Pkg.add("Ipopt")
#Pkg.add("MathOptInterface")
"""
using JuMP, Ipopt, MathOptInterface
include("enthalpy.jl")
const MOI = MathOptInterface


Base.@kwdef mutable struct cond_par
    outPurge_T::Float64 = 313.0;
    outProduct_T::Float64 = 313.0;
    splitratio::Vector = [1.0, 0.001, 1.0, 1.0, 1.0]; # Ideal split ratio, switch later
end

Base.@kwdef mutable struct _par
  cond::cond_par=cond_par();
  hconst = heavy_const;
end

par = _par()
#"""

function Cond_model(model, par)

  @variable(model, 0 <= cond_in_mol[1:5]); # Stream 11
  @variable(model, 0 <= cond_outPurge_mol[1:5]); # Stream 13 (H2O stream)
  @variable(model, 0 <= cond_outProduct_mol[1:5]); # Stream 12 (PSA )

  ini_cond_in = [0.00032515179570229987, 55.89429423256371, 656.2769122396338, 23.343132770546816, 188.85192754301332];  # CH4 H2O H2 CO CO2
  ini_cond_outPurge = [0.0, 55.83839993833114, 0.0, 0.0, 0.0];
  ini_cond_outProduct = [0.00032515179570229987, 0.05589429423256371, 656.2769122396338, 23.343132770546816, 188.85192754301332];

  for i=1:5
    set_start_value(cond_in_mol[i] , ini_cond_in[i]);
    set_start_value(cond_outPurge_mol[i], ini_cond_outPurge[i]);
    set_start_value(cond_outProduct_mol[i], ini_cond_outProduct[i]);
  end

  @variable(model, 273 <= cond_in_T, start = 313.00);
  @variable(model, 273 <= cond_outPurge_T, start = 313.00);
  @variable(model, 273 <= cond_outProduct_T, start = 313.00)

  cond_H_outPurge = build_enthalpy(model, cond_outPurge_T, par)
  cond_H_outProduct = build_enthalpy(model, cond_outProduct_T, par)
  cond_H_in = build_enthalpy(model, cond_in_T, par)


  # Constraints
  # Mass balance
  @NLconstraint(model, par.cond.splitratio[1]*cond_in_mol[1] - cond_outProduct_mol[1] == 0)
  @NLconstraint(model, par.cond.splitratio[2]*cond_in_mol[2] - cond_outProduct_mol[2] == 0)
  @NLconstraint(model, par.cond.splitratio[3]*cond_in_mol[3] - cond_outProduct_mol[3] == 0)
  @NLconstraint(model, par.cond.splitratio[4]*cond_in_mol[4] - cond_outProduct_mol[4] == 0)
  @NLconstraint(model, par.cond.splitratio[5]*cond_in_mol[5] - cond_outProduct_mol[5] == 0)
  @NLconstraint(model, (1-par.cond.splitratio[1])*cond_in_mol[1] - cond_outPurge_mol[1] == 0)
  @NLconstraint(model, (1-par.cond.splitratio[2])*cond_in_mol[2] - cond_outPurge_mol[2] == 0)
  @NLconstraint(model, (1-par.cond.splitratio[3])*cond_in_mol[3] - cond_outPurge_mol[3] == 0)
  @NLconstraint(model, (1-par.cond.splitratio[4])*cond_in_mol[4] - cond_outPurge_mol[4] == 0)
  @NLconstraint(model, (1-par.cond.splitratio[5])*cond_in_mol[5] - cond_outPurge_mol[5] == 0)

  # Energy balance
  #@NLconstraint(model, sum(cond_H_outPurge[i]*cond_outPurge_mol[i] - cond_H_in[i]*(1-par.cond.splitratio[i])*cond_in_mol[i] for i=1:5) == 0);
  #@NLconstraint(model, sum(cond_H_outProduct[i]*cond_outProduct_mol[i] - cond_H_in[i]*par.cond.splitratio[i]*cond_in_mol[i] for i=1:5) == 0);

  @NLconstraint(model, cond_outPurge_T - par.cond.outPurge_T == 0);
  @NLconstraint(model, cond_outProduct_T - par.cond.outProduct_T == 0);

  return model;
end

"""
m = Model(Ipopt.Optimizer);
m = Cond_model(m, par);
#optimize!(m);
d = JuMP.NLPEvaluator(m);
MOI.initialize(d, [:Jac]);
constraint_values = zeros(1,14);
inp = [0.00032515179570229987, 55.89429423256371, 656.2769122396338, 23.343132770546816, 188.85192754301332, 0.0, 55.83839993833114, 0.0, 0.0, 0.0, 0.00032515179570229987, 0.05589429423256371, 656.2769122396338, 23.343132770546816, 188.85192754301332, 313.00, 313.00, 313.00];
MOI.eval_constraint(d, constraint_values, inp[:]);
@show constraint_values
#"""