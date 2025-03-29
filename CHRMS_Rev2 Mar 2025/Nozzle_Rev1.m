%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nozzle_Rev1
% 2024/08/06
% Adrian Comisso
%
% Desccription: 
% This function calcualte the nozzle exit velocity
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


%define that agibraic mach number function
%****for any expansion ratio greater than one two solutions exit to this
%****equation the higher one is the desired mach number.
%****the lower one is subsonic and proably represents a bernouli case which is
%****invalid for choked flow so its ignored.
%*****check validity of equation
exp_ratio_relation = @(M) ((state.gamma+1)/2)^(-(state.gamma+1)/(2*(state.gamma-1)))*...
                     (1+(state.gamma-1)/2*M^2)^((state.gamma+1)/(2*(state.gamma-1)))/...
                     M-param.exp_ratio_noz;
%solve for Mach at exit
%****solve using fzeros rather than iteration****
%****will always find second one cuz exp ratio is always above mach number(prove using desmos)
M_exit = fzero(exp_ratio_relation,param.exp_ratio_noz);
%find pressure at exit (Pa)
P_exit = state.P_cmbr*(1+0.5*(state.gamma-1)*M_exit^2)^(-state.gamma/(state.gamma-1));
%****forget what this is called not cstar *****
Cf              = sqrt(((2*state.gamma^2)/(state.gamma-1))*(2/(state.gamma+1))^((state.gamma+1)/...
                  (state.gamma-1))*(1-(P_exit/state.P_cmbr)^((state.gamma-1)/state.gamma)))+...
                  ((P_exit-param.P_atm)*(0.25*pi*param.d_noz_throat^2*param.exp_ratio_noz))/...
                  (state.P_cmbr*0.25*pi*param.d_noz_throat^2);
%find motor thrust
if strcmp(param.model_noz,'1D')
    %compute F_rocket (N)
    state.F_rocket = param.eff_noz*Cf*0.25*pi*param.d_noz_throat^2*state.P_cmbr*param.Cd_noz;
else %param.half_angle_noz = 1
    F_rocket_ideal = Cf*0.25*pi*param.d_noz_throat^2*state.P_cmbr;
    %compute devergence factor
    lambda = (1/2)*(1+cosd(param.alpha_noz));
    %find velocity component of F_rocket (N) ****reverse F_rocket to apply devergence to only velocity component ******
    F_rocket_vel = F_rocket_ideal - (0.25*pi*param.d_noz_throat^2*param.exp_ratio_noz)*(P_exit-param.P_atm);
    %compute F_rocket (N)
    state.F_rocket = param.eff_noz*param.Cd_noz*(F_rocket_vel*lambda + (0.25*pi*param.d_noz_throat^2*param.exp_ratio_noz)*(P_exit-param.P_atm));
end

%correct end of burn thrust
if state.F_rocket < 0
    state.F_rocket = 0;
end

%what is this???***************
if isreal(M_exit) == 0
    disp("m exit")
end

end