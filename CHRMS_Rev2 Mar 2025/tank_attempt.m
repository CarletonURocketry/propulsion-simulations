

function [o] = tank_calcs(x)%*maybe this doent need to be function*****

%enter vapor only phase 
%should this trigger with liqid still in sys or negative liquid?
if x.m_l > m_l_unvaporized
    %burnout
end

%get oxidizer porperties
oxidizer = nitrous_properties(x.T_tank);
%find heat removed by the vaporising oxidizer(J)***might be sensative to initial value
%*** potnetially calcualte volme of injector manifold and find m given dens
%*** potentially bernoili out the mas flow as if injector ambient pressure
%*** calibrate startup with cold flow data
%of that area for first iteration*****
delta_Q = x.delta_m_v * oxidizer.delta_H;
%find change in tank temperature (K)
delta_T = -delta_Q/(x.m_l*oxidizer.Cp);%****say remaining liquid???***
%find the new oxidizer liquid temperature (K)
x.T_tank = x.T_tank - delta_T; %***seems like temp should equalize ?

%get oxidizer porperties %should i be getting this again?
oxidizer = nitrous_properties(x.T_tank);

%model 1 ******************************************
%calcualte loss factor ***see pag 4 nitrous run tank point 1 (will vary)***
D = s.K/((s.N_inj*s.A_inj)^2);
%calcualte delta_P (Pa)
x.delta_P = oxidizer.P_v - x.P_chamber;%***tank? chamber? oxidizer? which
%calucalte m_dot_oxidizer
x.m_dot_l = sqrt((2*oxidizer.rho_l*x.delta_P)/D);

%model 2 *******************************************
%model 3 ************************************

%calucalte mass that exited tank this timestep
delta_m_l = x.m_dot_l * s.delta_t;
%calcualte the reaming mass of liquid if none of it vaporized
m_l_unvaporized = x.m_l - delta_m_l;
%calcualte the remaining mass of oxidizer in the tank
x.m_total = x.m_total - delta_m_l;
%calcualte the remaing mass of liquid if some does vaporize
x.m_l = (s.V_tank - (x.m_total/oxidizer.rho_v))/((1/oxidizer.rho_l)-(1/oxidizer.rho_v));
%calcualte mass vaporized
x.delta_m_v = m_l_unvaporized - x.m_l;

%%%%%%if vapor phase
gamma_nitrous = 1.3;%****put this in a nitrous stucture or in the nitrous propeteis


%question about assumption in pg 5 first para******
%comapre using a textbook blowdown from t atmo*****

T1 = x.T_Tank
rho1 = x.rho_v
Z1 = oxidizer.Z
P1 = x.P_Tank

%for loop maksure there is a lcause about blowing up in case its unstable
%find Z0 ***z0 line is from critical point to p=0 Z= 1
%extrapolate new Z2 ***not sure if this is what was intended by paper***
%loop
    %find T2
    %find P2
    %find Z2
    %if delta Z small enough exit loop
    %if number of loops is to big exit code and say loop blew up

%calcualte rho2



end