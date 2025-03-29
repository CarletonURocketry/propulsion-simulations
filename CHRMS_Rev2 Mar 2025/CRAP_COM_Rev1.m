%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COM_Rev1
% 2024/08/06
% Adrian Comisso
%
% Desccription: 
%
% Inputs:
% param - stores the static parameters of the scenario
% state - stores the current state of each varible
%
% Outputs:
% new_state - returns the updated state of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state] = COM(param,state)

%find CG of oxidizer liquid (m)
state.cg_ox_l = ((state.V_ox_l/(0.25*pi*param.d_tank^2))/2) + param.l_tank_bottom;
%find CG of oxidizer vapor (m)
state.cg_ox_v = ((state.V_ox_v/(0.25*pi*param.d_tank^2))/2) + ((state.V_ox_l/(0.25*pi*param.d_tank^2))/2) + param.l_tank_bottom;
%find CG of fuel (m)
state.cg_fuel = param.l_grn_bottom + (param.l_grn/2);
%find CG of gas in chamber (m)
state.cg_gas = state.cg_fuel;
%find total mass (kg)
state.m_motor = param.m_dry_motor + state.m_ox_l + state.m_ox_v + state.m_gas + state.m_fuel;
%find total CG from nozzle (m)
state.cg_motor = ((param.m_dry_motor*param.cg_dry_motor) + (state.m_ox_l*state.cg_ox_l) + (state.m_ox_v*state.cg_ox_v) ...
           + (state.m_gas*state.cg_gas) + (state.m_fuel*state.cg_fuel))/state.m_motor;

end