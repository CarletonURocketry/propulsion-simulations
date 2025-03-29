%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regression
% 2024/08/05
% Adrian Comisso
%
% Desccription: 
% This function calcualtes the regression rate the id of the grain and the
% mass flow of the fuel away from the grain.  
%
% Inputs:
% state - a struct that stores the current state of each varible
%
% Outputs:
% state - returns the updated state of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state] = Regression(state)

%calculate curent surface area of grain (m^2)
A_port_grn = pi*(state.id_grn/2)^2;

if strcmp(state.model_reg,'shift_OF')
    %calculate oxidizer mass flux on surface of grain (kg/s*m^2)
    m_ox_flux = state.m_dot_inj/A_port_grn;
    %calucalte regression rate (m/s)****rocket bible pg 262****
    state.r_dot_grn = 0.001*state.a*(m_ox_flux^state.n)*(state.l_grn^state.m);
    %calcualte mass flow of fuel away from grain (kg/s)
    state.m_dot_fuel = state.rho_fuel*state.r_dot_grn*pi*state.id_grn*state.l_grn;
    %calcualte oxidizer to fuel ratio
    state.OF = state.m_dot_inj/state.m_dot_fuel;
else %state.model_reg = const_OF
    %calcualte mass flow of fuel away from grain (kg/s)
    state.m_dot_fuel = state.m_dot_inj/state.const_OF;
    %calucalte regression rate (m/s)
    state.r_dot_grn = state.m_dot_fuel/(state.rho_fuel*pi*state.id_grn*state.l_grn);
end

%calculate grain inner diameter
state.id_grn = state.id_grn + 2*state.r_dot_grn*state.dt;
%calucalte remaining mass of fuel
state.m_fuel = state.m_fuel - state.m_dot_fuel*state.dt;

end