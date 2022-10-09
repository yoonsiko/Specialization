#import Pkg
#Pkg.add("JuMP")
#Pkg.add("Ipopt")
#Pkg.add("MathOptInterface")
"""
using JuMP, Ipopt, MathOptInterface
include("enthalpy.jl")
const MOI = MathOptInterface


Base.@kwdef mutable struct ghr_par
    out_T::Float64 = 973.0;
end

Base.@kwdef mutable struct _par
  ghr::ghr_par=ghr_par();
  hconst = heavy_const;
end

par = _par()
"""
function GHR_model(model, par)

  # Variables
  # CH4, H2O, H2, CO, CO2
  @variable(model, 0 <= ghr_in_mol[1:5]); # Stream 6
  @variable(model, 0 <= ghr_out_mol[1:5]); # Stream 7

  #ini_ghr_in = [178.58 295.52 59.49 0.33 33.28];  # CH4 H2O H2 CO CO2
  ini_ghr_in = [178.58124131020057, 295.5198912988238, 59.489482856563995, 0.33375721522392365, 33.28038693993136];
  #ini_ghr_out = [22.54 122.16 544.93 139.05 50.60];
  ini_ghr_out = [22.543406840790453, 122.15875765109533, 544.9262854431126, 139.04829250631568, 50.60368611824971];

  for i=1:5
      set_start_value(ghr_in_mol[i] , ini_ghr_in[i]);
      set_start_value(ghr_out_mol[i] , ini_ghr_out[i]);
  end

  @variable(model, 273 <= ghr_in_T, start = 973);
  @variable(model, 273 <= ghr_out_T, start = 973);

  @variable(model, 0 <= ghr_Q, start = 34472.61482403047);

  # Expressions

  ghr_Ksmr_model = @NLexpression(model, exp(-22790 / ghr_out_T + 8.156 * log(ghr_out_T) - 4.421 / 10^3 * ghr_out_T

  - 4.330 * 10^3 / (ghr_out_T^2) - 26.030));

  ghr_Kwgsr_model = @NLexpression(model, exp(5693.5/ghr_out_T + 1.077*log(ghr_out_T) + 5.44e-4*ghr_out_T - 1.125e-7*ghr_out_T^2 - 49170/(ghr_out_T^2)-13.148));
  
  ghr_ksi_smr = @NLexpression(model, ghr_in_mol[1] - ghr_out_mol[1]);

  ghr_ksi_wgsr = @NLexpression(model, ghr_out_mol[5] - ghr_in_mol[5]);

  ghr_ntot = @NLexpression(model, sum(ghr_out_mol[i] for i=1:5));

  ghr_H_out = build_enthalpy(model, ghr_out_T, par);
  ghr_H_in = build_enthalpy(model, ghr_in_T, par);

  # Constraints
  # Mass balance

  @NLconstraint(model, ghr_Ksmr_model*((ghr_out_mol[1]/ghr_ntot) * (ghr_out_mol[2]/ghr_ntot )) -
   (((ghr_out_mol[4]/ghr_ntot) * (ghr_out_mol[3]/ghr_ntot)^3)) == 0);

  @NLconstraint(model, ghr_Kwgsr_model*((ghr_out_mol[4]/ghr_ntot) * (ghr_out_mol[2]/ghr_ntot)) -
   (((ghr_out_mol[5]/ghr_ntot) * (ghr_out_mol[3]/ghr_ntot))) == 0);

  @NLconstraint(model, ghr_out_mol[2] - ghr_in_mol[2] + ghr_ksi_smr + ghr_ksi_wgsr == 0);

  @NLconstraint(model, ghr_out_mol[3] - ghr_in_mol[3] - 3 * ghr_ksi_smr - ghr_ksi_wgsr == 0);

  @NLconstraint(model, ghr_out_mol[4] - ghr_in_mol[4] - ghr_ksi_smr + ghr_ksi_wgsr == 0);

  # Energy balance

  @NLconstraint(model, sum(ghr_H_out[i]*ghr_out_mol[i] - ghr_H_in[i]*ghr_in_mol[i] for i=1:5) - ghr_Q==0);

  # energy balance equipment specification

  @NLconstraint(model, ghr_out_T - par.ghr.out_T == 0);

  return model;

end
                
"""
m = Model(Ipopt.Optimizer);
m = GHR_model(m, par);
#optimize!(m);
d = JuMP.NLPEvaluator(m);
MOI.initialize(d, [:Jac]);
constraint_values = [0.0 0.0 0.0 0.0 0.0 0.0 0.0];
inp = [178.58124131020057; 295.5198912988238; 59.489482856563995; 0.33375721522392365;
 33.28038693993136; 22.543406840790453; 122.15875765109533; 544.9262854431126; 139.04829250631568; 50.60368611824971; 973; 973; 34472.61482403047];
MOI.eval_constraint(d, constraint_values, inp);
"""