%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main
% 2025/03/28
% Adrian Comisso
%
% Desccription: 
% This function computes any nessecary initial conditons and
% calls upon the Simulation_Loop function.
% 
% Inputs:
% state - a struct that stores the current state of each varible
% fuel - a struct that stores data on the type of fuel selected
%
% Outputs:
% output - a struct that stores all the past states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = Main(state,fuel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inital Condition Calcualtions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Simultion %%%%%
%time (s)
state.t = 0;

%%%%% Other Constants %%%%%
%ideal gas constant (J/mol*K)
state.R = 8.3145;
%acceleration due to gravity (m/s^2)
state.a_g = 9.81;

%%%%% Tank %%%%%
%test for which fluid is used
if strcmp(state.ox_name,'NitrousOxide')
    %define limits of saturation dome
    state.T_crit = 309.57;
    state.T_triple = 183;
    %specific gas constant (J/kg*K) from engineering toolbox
    state.R_sp = 188.91;
elseif strcmp(state.ox_name,'CarbonDioxide') && state.model_tank == 2
    %define limits of saturation dome
    state.T_crit = 304.25;
    state.T_triple = 216.738;
    %specific gas constant (J/kg*K) from engineering toolbox
    state.R_sp = 188.92;
else %err HRAP model with non nitrous oxide
    %notify user
    fprintf('Error you a have entered an invalid oxidizer or you have selected to use (HRAP tank method) which only works with NitrousOxide')
    %end simulation prematurley
    return
end    
%which tank method is used
if state.model_tank <= 1 
    %get oxidizer properties
    oxidixer = Nitrous_Properties(state.T_tank);
    %oxidizer liquid density (kg/m^3)
    state.rho_ox_l = oxidixer.rho_l;
    %oxidizer vapor density (kg/m^3)
    state.rho_ox_v = oxidixer.rho_v;
    %tank pressure (Pa)
    state.P_tank = oxidixer.P_v;
    %volume of liquid oxidizer (m^3)
    state.V_ox_l = (state.m_ox_total - (state.rho_ox_v*state.V_tank))/(state.rho_ox_l-state.rho_ox_v);
    %volume of vapor oxidizer (m^3)
    state.V_ox_v = state.V_tank - state.V_ox_l;
    %mass of liquid oxidizer (kg)
    state.m_ox_l = state.rho_ox_l * state.V_ox_l;
    %mass of vapor oxidizer (kg)
    state.m_ox_v = state.rho_ox_v * state.V_ox_v;
    %quality of oxidizer (percent of vapor by mass)
    state.quality_ox = state.m_ox_v/state.m_ox_total;
    %total internal energy of oxidizer (J)
    state.U_ox_total = 0;
    %liquid oxidizer specific enthalpy (J/kg)
    state.h_ox_l = 0;
    %vapor oxidizer specific enthalpy (J/kg)
    state.h_ox_v = 0;
else %tank_method == 2
    %oxidizer liquid density (kg/m^3)
    state.rho_ox_l = CoolProp.PropsSI('D','T',state.T_tank,'Q',0,state.ox_name);
    %oxidizer vapor density (kg/m^3)
    state.rho_ox_v = CoolProp.PropsSI('D','T',state.T_tank,'Q',1,state.ox_name);
    %quality of oxidizer (percent of vapor by mass)
    state.quality_ox = ((state.V_tank/state.m_ox_total)-(1/state.rho_ox_l))/((1/state.rho_ox_v)-(1/state.rho_ox_l));
    %tank pressure (Pa)
    state.P_tank = CoolProp.PropsSI('P','T',state.T_tank,'Q',state.quality_ox,state.ox_name);
    %volume of liquid oxidizer (m^3)
    state.V_ox_l = (state.m_ox_total - (state.rho_ox_v*state.V_tank))/(state.rho_ox_l-state.rho_ox_v);
    %volume of vapor oxidizer (m^3)
    state.V_ox_v = state.V_tank - state.V_ox_l;
    %mass of liquid oxidizer (kg)
    state.m_ox_l = state.rho_ox_l * state.V_ox_l;
    %mass of vapor oxidizer (kg)
    state.m_ox_v = state.rho_ox_v * state.V_ox_v;
    %total internal energy of oxidizer (J)
    state.U_ox_total = state.m_ox_total*CoolProp.PropsSI('U','T',state.T_tank,'Q',state.quality_ox,state.ox_name);
    %liquid oxidizer specific enthalpy (J/kg)
    state.h_ox_l = CoolProp.PropsSI('H','T',state.T_tank,'Q',0,state.ox_name);
    %vapor oxidizer specific enthalpy (J/kg)
    state.h_ox_v = CoolProp.PropsSI('H','T',state.T_tank,'Q',1,state.ox_name);
end
%set initial change in vapor mass (kg)
state.delta_m_ox_v = 0.001;
%set initial change in tank pressure (Pa)
state.delta_P_tank = 0;
%set initial change in tank Temperature (K)
state.delta_T_tank = 0;

%%%%% Vent %%%%%
%area of vent (m^2)
state.A_vent = pi*(state.d_vent/2)^2;
%delta P vent (Pa)
state.delta_P_vent = state.P_tank - state.P_atm;
%inital mass flow vent (kg/s)
state.m_dot_vent = 0;

%%%%% Injector %%%%%
%area of injector orifice (m^2)
state.A_inj = pi*(state.d_inj/2)^2;
%total area of injector (m^2)
state.A_inj_total = state.A_inj*state.N_inj;
%delta P injector (Pa)
state.delta_P_inj = state.P_tank - state.P_cmbr;
%initial mass flow of injector (kg/s)
state.m_dot_inj = 0;
%initial P2 critical (Pa) for HEM based injector modesl
state.P_2_crit = 0;

%%%%% Chamber %%%%%
%mass of fuel (kg)
state.m_fuel = state.rho_fuel*(state.l_grn*((pi*(state.od_grn/2)^2)-(pi*(state.id_grn/2)^2)));
%intial volume chamber with grains (kg)
state.V_cmbr_grn = state.V_cmbr_empty - 0.25*pi*(state.od_grn^2 - state.id_grn^2)*state.l_grn;
%intial mass of gas in the chamber (kg)
state.m_gas = state.rho_atm*state.V_cmbr_grn;
%initial shifting OF
state.OF = 0;
%C* cobustion efficiency
state.c_star = 0;
%initial fuel mass flow rate (kg/s)
state.m_dot_fuel = 0;
%initial regression rate (m/s)
state.r_dot_grn = 0;

%%%% Nozzle %%%%%
%initial nozzle mass flow (kg/s)
state.m_dot_noz = 0;
%initial combusted gas mass flow (kg/s)
state.m_dot_gas = 0;
%initital specific heat ratio
state.gamma = 0;
%initial exit mach number
state.M_exit = 0;
%initial exit pressure (Pa)
state.P_exit = state.P_atm;
%initial motor force (N)
state.F_rocket = 0;

%%%%% Center of Mass %%%%%
%find CG of oxidizer liquid (m)
state.cg_ox_l = ((state.V_ox_l/(0.25*pi*state.d_tank^2))/2) + state.l_tank_bottom;
%find CG of oxidizer vapor (m)
state.cg_ox_v = ((state.V_ox_v/(0.25*pi*state.d_tank^2))/2) + state.cg_ox_l;
%find CG of fuel (m)
state.cg_fuel = state.l_grn_bottom + (state.l_grn/2);
%find CG of gas in chamber (m)
state.cg_gas = state.cg_fuel;
%find total mass (kg)
state.m_motor = state.m_dry_motor + state.m_ox_l + state.m_ox_v + state.m_gas + state.m_fuel;
%find total CG from nozzle (m)
state.cg_motor = ((state.m_dry_motor*state.cg_dry_motor) + (state.m_ox_l*state.cg_ox_l) + (state.m_ox_v*state.cg_ox_v) ...
           + (state.m_gas*state.cg_gas) + (state.m_fuel*state.cg_fuel))/state.m_motor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%call loop
output = Simulation_Loop(state,fuel);