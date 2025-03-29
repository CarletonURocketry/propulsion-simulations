%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chamber_Rev1
% 2024/08/06
% Adrian Comisso
%
% Desccription: 
% This function determines the chamber pressure.
%****check against tatmo book*****
%
% Inputs:
% param - stores the static parameters of the scenario
% state - stores the current state of each varible
%
% Outputs:
% new_state - returns the updated state of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state] = chamber(param,state)

%seems like there might be some redundant paths to certain values in the
%chamber make sure they all agree. eg mass flow by sum of mass flows or by
%regression volume and density****************************

%calculate chamber volume (m^3)
state.V_cmbr_grn = param.V_cmbr_empty - 0.25*pi*(param.od_grn^2 - state.id_grn^2)*param.l_grn;

%calculate change in chamber volume (m^3/s)***HRAP uses form last time step to
%current rather than current to next?*******************
V_dot_cmbr_grn = 0.25*pi*(state.id_grn^2 - (state.id_grn - 2*state.r_dot_grn*param.dt)^2)*param.l_grn/param.dt;

%calcualte mss flow through nozzle (kg/s)
state.m_dot_noz = state.P_cmbr*param.Cd_noz*(0.25*pi*param.d_noz_throat^2)/state.c_star;

%calcualte the rate of increase of gas in the chamber(kg/s)
state.m_dot_gas = state.m_dot_fuel + state.m_dot_inj - state.m_dot_noz;

%calcualte the mass of the gas in the chamber (kg)
state.m_gas = state.m_gas + state.m_dot_gas*param.dt;

%calcualte the change in pressure (Pa/s) ******what is happening here?
P_dot_cmbr = state.P_cmbr*(state.m_dot_gas/state.m_gas - V_dot_cmbr_grn/state.V_cmbr_grn);

%calcualte the new chamber pressure (Pa)
state.P_cmbr = state.P_cmbr + P_dot_cmbr*param.dt;

%make sure it doesnt go below atm************ is this needed?
if state.P_cmbr <= param.P_atm
    state.P_cmbr = param.P_atm;
    state.mdot_n = 0;
end

end