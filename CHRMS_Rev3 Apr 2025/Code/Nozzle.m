%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nozzle
% 2024/08/06
% Adrian Comisso
%
% Desccription: 
% This function calcualte the nozzle exit velocity
%
% Inputs:
% state - a struct that stores the current state of each varible
%
% Outputs:
% state - returns the updated state of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state] = Nozzle(state)


%define that agibraic mach number function
exp_ratio_relation = @(M) ((state.gamma+1)/2)^(-(state.gamma+1)/(2*(state.gamma-1)))*...
                     (1+(state.gamma-1)/2*M^2)^((state.gamma+1)/(2*(state.gamma-1)))/...
                     M-state.exp_ratio_noz;
%solve for Mach at exit
M_exit = fzero(exp_ratio_relation,state.exp_ratio_noz);
%find pressure at exit (Pa)
P_exit = state.P_cmbr*(1+0.5*(state.gamma-1)*M_exit^2)^(-state.gamma/(state.gamma-1));
%thrust coefficient
Cf              = sqrt(((2*state.gamma^2)/(state.gamma-1))*(2/(state.gamma+1))^((state.gamma+1)/...
                  (state.gamma-1))*(1-(P_exit/state.P_cmbr)^((state.gamma-1)/state.gamma)))+...
                  ((P_exit-state.P_atm)*(0.25*pi*state.d_noz_throat^2*state.exp_ratio_noz))/...
                  (state.P_cmbr*0.25*pi*state.d_noz_throat^2);
%find motor thrust
if strcmp(state.model_noz,'1D')
    %compute F_rocket (N)
    state.F_rocket = state.eff_noz*Cf*0.25*pi*state.d_noz_throat^2*state.P_cmbr*state.Cd_noz;
else %state.half_angle_noz = 1
    F_rocket_ideal = Cf*0.25*pi*state.d_noz_throat^2*state.P_cmbr;
    %compute devergence factor
    lambda = (1/2)*(1+cosd(state.alpha_noz));
    %find velocity component of F_rocket (N)
    F_rocket_vel = F_rocket_ideal - (0.25*pi*state.d_noz_throat^2*state.exp_ratio_noz)*(P_exit-state.P_atm);
    %compute F_rocket (N)
    state.F_rocket = state.eff_noz*state.Cd_noz*(F_rocket_vel*lambda + (0.25*pi*state.d_noz_throat^2*state.exp_ratio_noz)*(P_exit-state.P_atm));
end

%correct end of burn thrust
if state.F_rocket < 0
    state.F_rocket = 0;
end

end