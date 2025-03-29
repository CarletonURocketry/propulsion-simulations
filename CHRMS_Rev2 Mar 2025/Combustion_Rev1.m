%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combustion_Rev1
% 2024/08/06
% Adrian Comisso
%
% Desccription: 
% This function calcualtes C* based on OF and P_chamber 
%****this method interpolates see about CEA, ProPrep or RPA or combustion
%toolbox??******************
%
% Inputs:
% param - stores the static parameters of the scenario
% state - stores the current state of each varible
%
% Outputs:
% new_state - returns the updated state of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state] = Combustion(param,state,fuel)

%Interepolate k, molar mass, temperature from tables in propellant config
%first two are axes, 3rd is desired varible fill of table, last 2 are input
%specific heat ratio (aka gamma)***********
state.gamma = interp2x(fuel.prop_OF,fuel.prop_Pc,fuel.prop_k,state.OF,state.P_cmbr);
%molar mas of mixture
%interepret is a custom matlab fucntion that this guy wrote*****
MM = interp2x(fuel.prop_OF,fuel.prop_Pc,fuel.prop_M,state.OF,state.P_cmbr);
%temparture of mixture
T = interp2x(fuel.prop_OF,fuel.prop_Pc,fuel.prop_T,state.OF,state.P_cmbr);

%calcualte the specific gas constant of the mixture
R_mixture = 8314.5/MM;
%find density using ideal gas law???
rho_mixture = state.P_cmbr/(R_mixture*T);
%find c star
state.c_star = param.c_star_eff*sqrt((R_mixture*T)/(state.gamma*(2/(state.gamma+1))^((state.gamma+1)/(state.gamma-1))));

end