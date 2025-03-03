%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nozzle_Rev1
% 2024/08/06
% Adrian Comisso
%
% Desccription: 
% This function calcualte the nozzle exit velocity
%****check against tatmo book*****
% *****do divergence stuff ******
% ****do unchoked regimes *****
%
% Inputs:
% param - stores the static parameters of the scenario
% state - stores the current state of each varible
%
% Outputs:
% new_state - returns the updated state of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state] = nozzle(param,state)

%define that unsolvable agibraic mach number function
%for any expansion ratio greater than one two solutions exit to this
%equation the higher one is the desired mach number.
%the lower one is subsonic and proably represents a bernouli case which is
%invalid for choked flow so its ignored.
%check validity of equation
exp_ratio_relation = @(M) ((state.gamma+1)/2)^(-(state.gamma+1)/(2*(state.gamma-1)))*...
                     (1+(state.gamma-1)/2*M^2)^((state.gamma+1)/(2*(state.gamma-1)))/...
                     M-param.exp_ratio_noz;
%solve using fzeros rather than iteration****
%will always find second one cuz exp ratio is always above mach number(prove using desmos)
M_exit = fzero(exp_ratio_relation,param.exp_ratio_noz);
%find pressure at exit
P_exit = state.P_cmbr*(1+0.5*(state.gamma-1)*M_exit^2)^(-state.gamma/(state.gamma-1));
%forget what this is called not cstar *****
Cf              = sqrt(((2*state.gamma^2)/(state.gamma-1))*(2/(state.gamma+1))^((state.gamma+1)/...
                  (state.gamma-1))*(1-(P_exit/state.P_cmbr)^((state.gamma-1)/state.gamma)))+...
                  ((P_exit-param.P_atm)*(0.25*pi*param.d_noz_throat^2*param.exp_ratio_noz))/...
                  (state.P_cmbr*0.25*pi*param.d_noz_throat^2);
%find thrust
state.F_rocket = param.eff_noz*Cf*0.25*pi*param.d_noz_throat^2*state.P_cmbr*param.Cd_noz;

if state.F_rocket < 0
    state.F_rocket = 0;
end

%what is this???***************
if isreal(M_exit) == 0
    disp("m exit")
end


end