%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main_Rev3
% 2024/08/01
% Adrian Comisso
%
% Desccription: 
% This code manages and initiates the inputs calls the main sim loop and
% plots data using fully rebuilt functions based of HRAP
%
% To Do:*************************************************************
% Unit conversions
% plots
% ui?
% add notes
% report
% mass calcs
%
% Referneces:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Clear Data %%%%%
clc
clearvars -except L

%%%%% variable storage objects %%%%%
%static parameters
param = struct();
%current state
state = struct(); 
%units **************************************
%units = struct(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Simulation Parameters %%%%%
%max sim time (s)
t_max = 20;
%time step (s)
dt = 0.001;
%cold flow
cold_flow = true;

%%%%% Atmospheric Properties %%%%%
%atmospheric pressure (Pa)
P_atm = 101325;
%atmospheric temperature (K)
T_atm = 293.15; 
%atmosphric density (kg/m^3)
rho_atm = 1.225;

%%%%% Other Constants %%%%%
%ideal gas constant (J/mol*K)
R = 8.3145;
%acceleration due to gravity (m/s^2)
a_g = 9.81;

%%%%% Tank Properties %%%%%
%tank model (note: method 0 only works with nitrous oxide)
%0 = HRAP Method (Nitrous Run Tank Emptying PDF)
%1 = Internal Energy Method (Internal Energy Feed System PDF)
model_tank = 0;
%oxidizer tank fluid
%'NitrousOxide'
%'CarbonDioxide'
ox_name = 'NitrousOxide';
%oxidizer tank volume (m^3)
V_tank = 0.00585;
%initial tank temperature (K)
T_tank = 285.52;%294.04;
%initial mass of oxidizer (kg)
m_ox_total = 4.309;

%%%%% Vent Properties %%%%%
%vent diameter (m)
d_vent = 0.00211;
%vent coefficient of discharge
Cd_vent = 0.6;
%vent state
%0 = no vent
%1 = vent vapor to atmosphere
%2 = vent vapor to chamber
state_vent = 0;
%vapoor venting model
%SPI_v = single phase invompressible vapor
%PGC_v = perfect gass compresible vapor
model_vent = 'PGC_v';

%%%%% Injector Properties %%%%%
%number of injector holes
N_inj = 60;
%diameter of injector holes (m)
d_inj = 0.0015;
%injector coefficient of discharge
Cd_inj = 0.2;
%liquid injector model
%SPI = single phase incompressible
%HEM = homgenous equilibrium model
%Dyer = Dyer model
model_inj_l = 'SPI';
%vapoor injector model
%SPI_v = single phase invompressible vapor
%PGC_v = perfect gass compresible vapor
model_inj_v = 'PGC_v';

%%%%% Chamber Properties %%%%%
%empty chamber volume (m^3)
V_cmbr_empty = 0.0098;
%initial combustion chamber pressure (Pa)
P_cmbr = P_atm;
%initial grain inner diameter (m)
id_grn = 0.0711;
%grain outer diameter (m)
od_grn = 0.1016;
%grain lenght (m)
l_grn = 0.9144;
%c_star efficiency
c_star_eff = 0.8;
%fuel grain density (kg/m^3)
rho_fuel = 900;
%regression coeficients
a = 0.0876;
n = 0.3953;
m = 0;

%%%%% Nozzle Properties %%%%%
%nozzle throat diameter (m)
d_noz_throat = 0.0348;
%nozzle expansion ratio
exp_ratio_noz = 4.0293;
%nozzle efficiency
eff_noz = 0.95;
%nozzle coefficient of discharge
Cd_noz = 0.9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit Conversions (could possibly switch with below?)****
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%convert units of everything***** needs function*****************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcualtions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Tank %%%%%
%test for which fluid is used
if strcmp(ox_name,'NitrousOxide')
    %define limits of saturation dome
    T_crit = 309.57;
    T_triple = 183;
    %specific gas constant (J/kg*K) from engineering toolbox
    R_sp = 188.91;
elseif strcmp(ox_name,'CarbonDioxide') && model_tank == 1
    %define limits of saturation dome
    T_crit = 304.25;
    T_triple = 216.738;
    %specific gas constant (J/kg*K) from engineering toolbox
    R_sp = 188.92;
else %if fluid isnt nitrous and your trying to use HRAP model
    %notify user
    fprintf('Error you a have entered an invalid oxidizer or you have selected to use (HRAP method) which only works with NitrousOxide')
    %end simulation prematurley
    return
end    
%which tank method is used
if model_tank == 0
    %get oxidizer properties
    oxidixer = Nitrous_Properties_Rev1(T_tank);
    %oxidizer liquid density (kg/m^3)
    rho_ox_l = oxidixer.rho_l;
    %oxidizer vapor density (kg/m^3)
    rho_ox_v = oxidixer.rho_v;
    %tank pressure (Pa)
    P_tank = oxidixer.P_v;
    %volume of liquid oxidizer (m^3)
    V_ox_l = (m_ox_total - (rho_ox_v*V_tank))/(rho_ox_l-rho_ox_v);
    %volume of vapor oxidizer (m^3)
    V_ox_v = V_tank - V_ox_l;
    %mass of liquid oxidizer (kg)
    m_ox_l = rho_ox_l * V_ox_l;
    %mass of vapor oxidizer (kg)
    m_ox_v = rho_ox_v * V_ox_v;
    %quality of oxidizer (percent of vapor by mass)
    quality_ox = m_ox_v/m_ox_total;
    %total internal energy of oxidizer (J)
    U_ox_total = 0;
    %liquid oxidizer specific enthalpy (J/kg)
    h_ox_l = 0;
    %vapor oxidizer specific enthalpy (J/kg)
    h_ox_v = 0;
else %tank_method == 1
    %oxidizer liquid density (kg/m^3)
    rho_ox_l = CoolProp.PropsSI('D','T',T_tank,'Q',0,ox_name);
    %oxidizer vapor density (kg/m^3)
    rho_ox_v = CoolProp.PropsSI('D','T',T_tank,'Q',1,ox_name);
    %quality of oxidizer (percent of vapor by mass)
    quality_ox = ((V_tank/m_ox_total)-(1/rho_ox_l))/((1/rho_ox_v)-(1/rho_ox_l));
    %tank pressure (Pa)
    P_tank = CoolProp.PropsSI('P','T',T_tank,'Q',quality_ox,ox_name);
    %volume of liquid oxidizer (m^3)
    V_ox_l = (m_ox_total - (rho_ox_v*V_tank))/(rho_ox_l-rho_ox_v);
    %volume of vapor oxidizer (m^3)
    V_ox_v = V_tank - V_ox_l;
    %mass of liquid oxidizer (kg)
    m_ox_l = rho_ox_l * V_ox_l;
    %mass of vapor oxidizer (kg)
    m_ox_v = rho_ox_v * V_ox_v;
    %total internal energy of oxidizer (J)
    U_ox_total = m_ox_total*CoolProp.PropsSI('U','T',T_tank,'Q',quality_ox,ox_name);
    %liquid oxidizer specific enthalpy (J/kg)
    h_ox_l = CoolProp.PropsSI('H','T',T_tank,'Q',0,ox_name);
    %vapor oxidizer specific enthalpy (J/kg)
    h_ox_v = CoolProp.PropsSI('H','T',T_tank,'Q',1,ox_name);
end

%%%%% Vent %%%%%
%area of vent (m^2)
A_vent = pi*(d_vent/2)^2;
%delta P vent
delta_P_vent = P_tank - P_atm;

%%%%% Injector %%%%%
%area of injector orifice (m^2)
A_inj = pi*(d_inj/2)^2;
%total area of injector (m^2)
A_inj_total = A_inj*N_inj;
%delta P injector
delta_P_inj = P_tank - P_cmbr;

%%%%% Chamber %%%%%
%mass of fuel (kg)
m_fuel = rho_fuel*(l_grn*((pi*(od_grn/2)^2)-(pi*(id_grn/2)^2)));
%intial volume chamber with grains (kg)
V_cmbr = V_cmbr_empty - 0.25*pi*(od_grn^2 - id_grn^2)*l_grn;
%intial mass of gas in the chamber (kg)
m_gas = rho_atm*V_cmbr;

%%%%% Center of Mass %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign Varibles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Fuel File %%%%% *************needs to be changable********
path = 'C:\Users\adric\OneDrive\Desktop\4-CUInSpace\Hybrid Sim\HRAP_Decompiled_Matlab_Versions\CUInspace_HRAP_Rev1';
file = 'Paraffin_Rev1.mat';
load(fullfile(path, file));

%%%%% Static Parameters %%%%%
%simulation
param.t_max = t_max;
param.dt = dt;
%tank
param.cold_flow = cold_flow;
param.model_tank = model_tank;
param.V_tank = V_tank;
param.ox_name = ox_name;
param.T_crit = T_crit;
param.T_triple = T_triple;
param.R_sp = R_sp;
%vent
param.A_vent = A_vent;
param.Cd_vent = Cd_vent;
param.state_vent = state_vent;
param.model_vent = model_vent;
%injector
param.A_inj_total = A_inj_total;
param.Cd_inj = Cd_inj;
param.model_inj_l = model_inj_l;
param.model_inj_v = model_inj_v;
%chamber
param.V_cmbr = V_cmbr_empty;
param.od_grn = od_grn;
param.l_grn = l_grn;
param.c_star_eff = c_star_eff;
param.rho_fuel = rho_fuel;
param.a = a;
param.n = n;
param.m = m;
%nozzle
param.d_noz_throat = d_noz_throat;
param.exp_ratio_noz = exp_ratio_noz;
param.eff_noz = eff_noz;
param.Cd_noz = Cd_noz;
param.P_atm = P_atm;
param.T_atm = T_atm;
param.R = R;
param.a_g = a_g;

%%%%% Initial Conditions %%%%%
%tank
state.T_tank = T_tank;
state.P_tank = P_tank;
state.U_ox_total = U_ox_total;
state.rho_ox_l = rho_ox_l;
state.rho_ox_v = rho_ox_v;
state.m_ox_total = m_ox_total;
state.V_ox_l = V_ox_l;
state.V_ox_v = V_ox_v;
state.m_ox_l = m_ox_l;
state.m_ox_v = m_ox_v;
state.quality_ox = quality_ox;
state.delta_m_ox_v = 0.001; %set to non zero?*****
state.delta_T_tank = 0;
state.h_ox_l = h_ox_l;
state.h_ox_v = h_ox_v;
%injector
state.delta_P_vent = delta_P_vent;
state.m_dot_inj = 0;
state.delta_P_inj = 0;
state.m_dot_vent = 0;
state.P_2_crit = 0;
%chamber
state.m_fuel = m_fuel;
state.P_cmbr = P_cmbr;
state.OF = 0;
state.c_star = 0;
state.m_dot_fuel = 0;
state.id_grn = id_grn; 
state.r_dot_grn = 0;
%nozzle
state.m_dot_noz = 0;
state.m_dot_gas = 0;
state.m_gas = m_gas;
state.V_cmbr_grn = V_cmbr;
state.gamma = 0;
state.M_exit = 0;
state.P_exit = P_atm;
state.F_rocket = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%call loop
output = Simulation_Loop_Rev1(param,state,paraffin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%display end condition
disp(output.end_cond)

hold on
figure(1);
plot(output.t,output.P_tank)
hold on
figure(2);
plot(output.t,output.T_tank)
hold on
figure(3);
plot(output.t,output.m_dot_inj)
hold on
figure(4);
plot(output.t,output.m_ox_total)
hold on
figure(4);
plot(output.t,output.thrust)

%calcualte and give report isp, class motor avg thrust etc