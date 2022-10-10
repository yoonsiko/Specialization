# Order of components: CH4, H2O, H2, CO, CO2, C2H6, C3H8, n-C4H10, i-C4H10, C5+

# Order of coefficients: A, B, C, H298, Cp298/R <- remember to multiply with 8.314
using JuMP, Ipopt
heavy_const = [1.702 9.081 -2.164 -74.520 4.217; 
                3.470 1.450 0.000 -241.818 4.038;
                3.249 0.422 0.000 0.000 3.468;
                3.376 0.557 0.000 -110.525 3.507;
                5.457 0.557 0.000 -393.509 4.467;
                1.131 19.225 -5.561 -83.820 6.369;
                1.213 28.785 -8.824 -104.680 9.011;
                1.935 36.915 -11.402 -125.790 11.298;
                1.935 36.915 -11.402 -125.790 11.298;
                2.464 45.351 -14.111 -146.760 14.731]
"""
function heavy_enthalpy(T, v=heavy_const)
    h_vec = v[:,4] + (((v[:,1]*T + v[:,2]/2*T^2*10^(-3) + v[:,3]/3*T^3*10^(-6)) 
    - (v[:,1]*298 + v[:,2]/2*298^2*10^(-3) + v[:,3]/3*298^3*10^(-6)))*8.314 /1000)

    return h_vec
end
"""

#print(heavy_enthalpy(298))
using JuMP

using Ipopt 
"""
Structure for parameters specific to GHR
"""
Base.@kwdef mutable struct ghr_par
    out_T::Float64 = 973.0
end

"""
Structure for parameters of flow sheet
"""
Base.@kwdef mutable struct _par
    ghr::ghr_par=ghr_par();
    hconst = heavy_const
end


par = _par()
"""
Build nonlinear expression for enthalpy calculation.
Pass the symbolic temperature.
"""
function build_enthalpy(model, T, par)
    return @NLexpression(model, [i = 1:5], par.hconst[i,4] + (((par.hconst[i,1]*T + par.hconst[i,2]/2*T^2*10^(-3) + par.hconst[i,3]/3*T^3*10^(-6)) 
    - (par.hconst[i,1]*298 + par.hconst[i,2]/2*298^2*10^(-3) + par.hconst[i,3]/3*298^3*10^(-6)))*8.314 /1000))
end

model = Model(Ipopt.Optimizer)


function GHR_model(model, par)

    # Variables
    # Ch4, H20, H2, CO, CO2
    @variable(model, 0 <= ghr_in_mol[1:5]);
    @variable(model, 0 <= ghr_out_mol[1:5], start = 1);

    const ini_ghr_in = [1. 2. 3. 4. 5.];
    const ini_ghr_out = [1. 2. 3. 4. 5.];

    for i=1:5
        set_start_value(ghr_in_mol[i] , ini_ghr_in[i]);
        set_start_value(ghr_out_model[i], ini_ghr_out[i])
    end

    @variable(model, 273 <= ghr_in_T, start = 600);
    @variable(model, 273 <= ghr_out_T, start = 700);

    @variable(model, 0 <= ghr_Q, start = 1)

    # Expressions

    ghr_Ksmr_model = @NLexpression(model, exp(-22790 / ghr_out_T + 8.156 * log(ghr_out_T) - 4.421 / 10^3 * ghr_out_T

    - 4.330 * 10^3 / (ghr_out_T^2) - 26.030))

    ghr_Kwgsr_model = @NLexpression(model, exp(5693.5/ghr_out_T + 1.077*log(ghr_out_T) + 5.44e-4*ghr_out_T - 1.125e-7*ghr_out_T^2 - 49170/(ghr_out_T^2)-13.148))

    ghr_ksi_smr = @NLexpression(model, ghr_in_mol[1] - ghr_out_mol[1])

    ghr_ksi_wgsr = @NLexpression(model, ghr_out_mol[5] - ghr_in_mol[5])

    ghr_ntot = @NLexpression(model, sum(ghr_out_mol[i] for i=1:5))

    ghr_H_out = build_enthalpy(model, ghr_out_T, par)
    ghr_H_in = build_enthalpy(model, ghr_in_T, par)

    

    # Constraints
    
    # mass balance

    @NLconstraint(model, ghr_Ksmr_model*((ghr_out_mol[1]/ghr_ntot) * (ghr_out_mol[2]/ghr_ntot )) - (((ghr_out_mol[4]/ghr_ntot) * (ghr_out_mol[3]/ghr_ntot)^3)) == 0)

    @NLconstraint(model, ghr_Kwgsr_model*((ghr_out_mol[4]/ghr_ntot) * (ghr_out_mol[2]/ghr_ntot)) - (((ghr_out_mol[5]/ghr_ntot) * (ghr_out_mol[3]/ghr_ntot))) == 0)

    @NLconstraint(model, ghr_out_mol[2] - ghr_in_mol[2] + ghr_ksi_smr + ghr_ksi_wgsr == 0)

    @NLconstraint(model, ghr_out_mol[3] - ghr_in_mol[3] - 3 * ghr_ksi_smr - ghr_ksi_wgsr == 0)

    @NLconstraint(model, ghr_out_mol[4] - ghr_in_mol[4] - ghr_ksi_smr + ghr_ksi_wgsr == 0)

    # energy balance

    @NLconstraint(model, sum(ghr_H_out[i]*ghr_out_mol[i] - ghr_H_in[i]*ghr_in_mol[i] for i=1:5) - ghr_Q==0)

    # energy balance equipment specification

    @NLconstraint(model, ghr_out_T - par.ghr.out_T == 0)

    return model

end

model = GHR_model(model, par)
