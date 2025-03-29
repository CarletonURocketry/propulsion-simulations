%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COM
% 2024/08/06
% Adrian Comisso
%
% Desccription: 
%
% Inputs:
% state - a struct that stores the current state of each varible
%
% Outputs:
% state - returns the updated state of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state] = Mass(state)

%find CG of oxidizer liquid (m)
state.cg_ox_l = ((state.V_ox_l/(0.25*pi*state.d_tank^2))/2) + state.l_tank_bottom;
%find CG of oxidizer vapor (m)
state.cg_ox_v = ((state.V_ox_v/(0.25*pi*state.d_tank^2))/2) + ((state.V_ox_l/(0.25*pi*state.d_tank^2))/2) + state.l_tank_bottom;
%find CG of fuel (m)
state.cg_fuel = state.l_grn_bottom + (state.l_grn/2);
%find CG of gas in chamber (m)
state.cg_gas = state.cg_fuel;
%find total mass (kg)
state.m_motor = state.m_dry_motor + state.m_ox_l + state.m_ox_v + state.m_gas + state.m_fuel;
%find total CG from nozzle (m)
state.cg_motor = ((state.m_dry_motor*state.cg_dry_motor) + (state.m_ox_l*state.cg_ox_l) + (state.m_ox_v*state.cg_ox_v) ...
           + (state.m_gas*state.cg_gas) + (state.m_fuel*state.cg_fuel))/state.m_motor;

end