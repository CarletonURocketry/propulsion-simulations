%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chamber
% 2024/08/06
% Adrian Comisso
%
% Desccription: 
% This function determines the chamber pressure.
%
% Inputs:
% state - a struct that stores the current state of each varible
%
% Outputs:
% state - returns the updated state of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state] = Chamber(state)

%calculate chamber volume (m^3)
state.V_cmbr_grn = state.V_cmbr_empty - 0.25*pi*(state.od_grn^2 - state.id_grn^2)*state.l_grn;

%calculate change in chamber volume (m^3/s)***HRAP uses form last time step to
V_dot_cmbr_grn = 0.25*pi*(state.id_grn^2 - (state.id_grn - 2*state.r_dot_grn*state.dt)^2)*state.l_grn/state.dt;

%calcualte mss flow through nozzle (kg/s)
state.m_dot_noz = state.P_cmbr*state.Cd_noz*(0.25*pi*state.d_noz_throat^2)/state.c_star;

%calcualte the rate of increase of gas in the chamber(kg/s)
state.m_dot_gas = state.m_dot_fuel + state.m_dot_inj - state.m_dot_noz;

%calcualte the mass of the gas in the chamber (kg)
state.m_gas = state.m_gas + state.m_dot_gas*state.dt;

%calcualte the change in pressure (Pa/s) 
P_dot_cmbr = state.P_cmbr*(state.m_dot_gas/state.m_gas - V_dot_cmbr_grn/state.V_cmbr_grn);

%calcualte the new chamber pressure (Pa)
state.P_cmbr = state.P_cmbr + P_dot_cmbr*state.dt;

%make sure it doesnt go below atm
if state.P_cmbr <= state.P_atm
    state.P_cmbr = state.P_atm;
    state.mdot_n = 0;
end

end