function Build_separator(SSEl, par)

    #Variables
    @variable(SSEl, I, start = 500);
    @variable(SSEl, U_rev, start = 1.4);
    @variable(SSEl, U);
    @variable(SSEl, 0 <= n_H2_i, start = 1);
    @variable(SSEl, 0 <= n_O2_i, start = 1);
    @variable(SSEl, 0 <= n_c_H2_o, start = 1);
    @variable(SSEl, 0 <= n_a_H2_o, start = 0.01);
    @variable(SSEl, 0 <= n_c_O2_o, start = 0.01);
    @variable(SSEl, 0 <= n_a_O2_o, start = 1);
    @variable(SSEl, 0 <= n_sepa_H2_o, start = 0.01);
    @variable(SSEl, 0 <= n_sepa_O2_o, start = 1);
    @variable(SSEl, 0 <= n_sepc_H2_o, start = 0.01);
    @variable(SSEl, 0 <= n_sepc_O2_o, start = 1);
    @variable(SSEl, TEl_i, start = 340);
    @variable(SSEl, Tcw_o, start = 305);
    @variable(SSEl, par[:minlye] <= mel_lye_i <= par[:maxlye], start = 10);
    @variable(SSEl, 0 <= m_cw <= 26.6, start = 5);
    @variable(SSEl, par[:minlye] <= mel_lye_o <= par[:maxlye], start = 9.9);
    @variable(SSEl, par[:minlye] <= mbu_lye_i <= par[:maxlye] , start = 9.9);
    @variable(SSEl, 0 <= mmu_H2O <= 5, start = 0.1); #Water makeup
    @variable(SSEl, par[:minP] <= P <= par[:maxP], start = 4);
    @variable(SSEl, θ);

    #Expressions
    @expression(SSEl, n_r_H2, I*(1-θ)*par[:ael]/(2*par[:F]));
    @expression(SSEl, n_r_O2, I*(1-θ)*par[:ael]/(4*par[:F]));
    @expression(SSEl,S_H2, 1000/(0.018*7.1698*10^4)/(10^(3.14*0.312)));
    @expression(SSEl, S_O2, S_H2*0.5);
    @NLexpression(SSEl, n_diff_H2, P*S_H2*par[:DH2]/(6.18*500e-6)*par[:ael]);
    @NLexpression(SSEl, n_diff_O2, P*S_O2*par[:DH2]/(6.18*500e-6)*par[:ael]);
    @NLexpression(SSEl, n_buff_O2_i, 0.5*mbu_lye_i*S_O2*P/(par[:sepeff]*par[:rhol]));
    @NLexpression(SSEl, n_buff_H2_i, 0.5*mbu_lye_i*S_H2*P/(par[:sepeff]*par[:rhol]));
    @NLexpression(SSEl,Pcomp, (n_sepc_H2_o+n_sepa_O2_o)/par[:effcomp]*(1.4/0.4)*8.314*TEl*((25/P)^(0.4/1.4)-1));
    @expression(SSEl, logT, (Tbuff + TEl_i)/2 - (Tcw_o + par[:Tcw])/2)  #arithmethic dT!
    @NLexpression(SSEl, HTO, n_a_H2_o/n_sepa_O2_o);
    @NLexpression(SSEl, etasys, n_c_H2_o*par[:LHV]*1e6/(par[:Uk]*I*par[:ael]+Pcomp));
    @NLexpression(SSEl, v, mel_lye_i/(par[:rhol]*3.63*0.1));

  
    #Separators
    @NLconstraint(SSEl, n_a_H2_o - n_sepa_H2_o == 0);
    @NLconstraint(SSEl, n_a_O2_o - n_sepa_O2_o - n_buff_O2_i == 0);
    @NLconstraint(SSEl, n_c_H2_o - n_sepc_H2_o - n_buff_H2_i == 0); 
    @NLconstraint(SSEl, n_c_O2_o - n_sepc_O2_o == 0);
    @NLconstraint(SSEl, par[:MH2]*(n_a_H2_o - n_sepa_H2_o) + par[:MO2]*(n_a_O2_o - n_sepa_O2_o - n_buff_O2_i) + 0.5*mel_lye_o - 0.5*mbu_lye_i == 0);
    @NLconstraint(SSEl, par[:MH2]*(n_c_H2_o - n_sepc_H2_o - n_buff_H2_i) + par[:MO2]*(n_c_O2_o - n_sepc_O2_o) + 0.5*mel_lye_o - 0.5*mbu_lye_i == 0);


    return SSEl
end 
