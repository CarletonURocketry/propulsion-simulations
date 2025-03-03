%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nitrous oxide
% 2024/07/31
% Adrian Comisso
%
% Desccription: 
% This function is designed to return the properties of nitrous oxide given
% the temperature the curves are found using paywalled data from the
% british equivalent of NIST
%
% Inputs:
% Input 1 - Temperature
%
% Output:
% Output 1 - Object containing: vapor pressure,
%                              liquid desnity,
%                              vapor density,
%                              specific enthalpy of vaporization,
%                              specific heat capacity
%                              vapor compressiblilit factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output] = nitrous_properties(T)

%%%%% Molecular Properties %%%%%
%molar mass
%source (NIST): https://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Mask=4#Thermo-Phase
output.MM = 0.0440124;
%gas constant (J/kg*K)
%source (Eng Tools): https://www.engineeringtoolbox.com/individual-universal-gas-constant-d_588.html
output.R = 188.91;
%specific heat ratio
%source (Eng Tools): https://www.engineeringtoolbox.com/specific-heat-ratio-d_608.html
output.gamma = 1.31;

%%%%% Critical Point %%%%%
%Source (NIST): https://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Mask=4#Thermo-Phase
%critical pressure (Pa)
P_crit = 7251000; %7238000; temprarlily made to match HRAP***********
%critical temperature (K)
T_crit = 309.57; %309.56;
%critical density (kg/m^3)
rho_crit = 452; %453.328;

%%%%% Vapor Pressure %%%%%
%vapor pressure curve fit varibles
a1 = -6.71893;
a2 = 1.35966;
a3 = -1.3779;
a4 = -4.051;
%vapor pressure (Pa)
output.P_v = P_crit*exp((1/(T/T_crit))*(a1*(1-T/T_crit) + ...
    a2*(1-(T/T_crit))^(3/2) + a3*(1-(T/T_crit))^(5/2) + a4*(1-(T/T_crit))^5));

%%%%% Density %%%%%
%liquid density curve fit varibles
b1 = 1.72328;
b2 = -0.83950;
b3 = 0.51060;
b4 = -0.10412;
%liquid density (kg/m^3)
output.rho_l = rho_crit*exp(b1*(1-(T/T_crit))^(1/3) + ...
    b2*(1-(T/T_crit))^(2/3) + b3*(1-(T/T_crit)) + b4*(1-(T/T_crit))^(4/3));
%liquid density curve fit varibles
c1 = -1.00900;
c2 = -6.28792;
c3 = 7.50332;
c4 = -7.90463;
c5 = 0.629427;
%vapor density (kg/m^3)
output.rho_v = rho_crit*exp(c1*((T_crit/T)-1)^(1/3) + ...
    c2*((T_crit/T)-1)^(2/3) + c3*((T_crit/T)-1) + c4*((T_crit/T)-1)^(4/3) + ...
    c5*((T_crit/T)-1)^(5/3));

%%%%% Enthalpy of Vaporization %%%%%

%liquid enthalpy curve fit varibles
d1 = -200;
d2 = 116.043;
d3 = -917.225;
d4 = 794.779;
d5 = -589.587;
%vapor enthalpy curve fit varibles
e1 = -200;
e2 = 440.055;
e3 = -459.701;
e4 = 434.081;
e5 = -485.338;
%enthalpy of vaporization (J)
output.delta_H = (e1-d1) + (e2-d2)*(1-(T/T_crit))^(1/3) + (e3-d3)* ...
    (1-(T/T_crit))^(2/3) + (e4-d4)*(1-(T/T_crit)) + (e5-d5)*(1-(T/T_crit))^(4/3);

%%%%% Specific Heat Capacity of Saturated Liquid %%%%%
%heat capacity curve fit varibles
f1 = 2.49973;
f2 = 0.023454;
f3 = -3.80136;
f4 = 13.0945;
f5 = -14.5180;
%specific heat capacity (J/kg*K) 
output.Cp = f1*(1 + f2*(1-(T/T_crit))^(-1) + f3* ...
    (1-(T/T_crit)) + f4*(1-(T/T_crit))^2 + f5*(1-(T/T_crit))^3);

%%%%% Compressibility Factor %%%%%
%note: see 'Nitrous decomperssion' aspire space page 5 for simpler method
%compresibility factor
output.Z = (P_crit*exp((1/(T/T_crit))*(a1*(1-T/T_crit) + ...
    a2*(1-(T/T_crit))^(3/2) + a3*(1-(T/T_crit))^(5/2) + a4*(1-(T/T_crit))^5)))/ ...
    ((rho_crit*exp(c1*((T_crit/T)-1)^(1/3) + c2*((T_crit/T)-1)^(2/3) + c3*((T_crit/T)-1)...
    + c4*((T_crit/T)-1)^(4/3) + c5*((T_crit/T)-1)^(5/3)))*output.R*T);

end