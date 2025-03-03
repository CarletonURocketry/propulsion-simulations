%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regression_Rev1
% 2024/08/05
% Adrian Comisso
%
% Desccription: 
% This function calcualtes the regression rate the id of the grain and the
% mass flow of the fuel away from the grain.  
%
% Inputs:
% param - stores the static parameters of the scenario
% state - stores the current state of each varible
%
% Outputs:
% new_state - returns the updated state of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state] = Regression(param,state)

%calculate curent surface area of grain (m^2)
A_port_grn = pi*(state.id_grn/2)^2;
%calculate oxidizer mass flux on surface of grain (kg/s*m^2)
m_ox_flux = state.m_dot_inj/A_port_grn;
%calucalte regression rate (m/s)
%rocket bible pg 262****
state.r_dot_grn = 0.001*param.a*(m_ox_flux^param.n)*(param.l_grn^param.m);
%calcualte mass flow of fuel away from grain (kg/s)
state.m_dot_fuel = param.rho_fuel*state.r_dot_grn*pi*state.id_grn*param.l_grn;
%calcualte oxidizer to fuel ratio
state.OF = state.m_dot_inj/state.m_dot_fuel;
%calculate grain inner diameter
state.id_grn = state.id_grn + 2*state.r_dot_grn*param.dt;
%calucalte remaining mass of fuel
state.m_fuel = state.m_fuel - state.m_dot_fuel*param.dt;

end