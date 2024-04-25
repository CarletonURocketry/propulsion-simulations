% CU InSpace Hybrid Rocket Engineering Model
% Academic Year 2022-2023
% Last update February 2, 2023
% Quentin Alexander
% quentinalexander@cmail.carleton.ca
% 
% This script simulates the flight of a N2O-paraffin hybrid rocket, in
% particular the propulsion system
% 
% Returns propellant burn rate, pressure and temperature over time in the
% tank and combustion chamber, thrust, Isp, and thrust to weight ratio, and
% altitude, velocity, and acceleration
% 
% Warning: the simulation will fail if tank temperature exceeds
% 36 degrees C (critical temperature of N2O)

clear all; %#ok<CLALL> 
close all;

mDot_ox_noise = 0.00; % random noise in oxidizer mass flow, approximately
                      % as a fraction of mass flow
mode = 0;             % 0 - flight; 1 - static fire; 2 - cold flow
ras = true;           % are we putting this through RAS aero

%% Vehicle Loading
m_ox = 16;                % [kg] mass of oxidizer loaded
r_port = 3.4375*0.0254/2; % [m]  initial fuel grain inner radius
% ^ optimal diameter ~ 3 7/16"

%% Constants
g = 9.80665;               % [m/s2]    acceleration due to gravity
                           %           (ISO 80000)
R_univ = 8314.46261815324; % [J/kmolK] universal gas constant
R_air = R_univ/28.97;      % [J/kgK]   specific gas constant of air

%% Launch Site Conditions
T_amb = 36;     % [deg C] ambient temperature
T_tank = 20;%T_amb; % [deg C] tank temperature
P_amb = 85E3;   % [Pa]    ambient pressure
elev = 1400;    % [m]     elevation of launch site above sea level

P_sea = P_amb/(1-elev*2.25577E-5)^5.25588; % [Pa] pressure at sea level
                                           %      formula from engineer's
                                           %      toolbox

%% Vehicle Mass
m_tank = 4.45;%4.0;                % [kg] mass of the tank
m_dry = 91*0.45359237 + m_tank;    % [kg] dry mass of the rocket

%% Geometry
% aerostructure
r_vehicle = 6.25*0.0254/2; % [m] radius of vehicle
N_fin = 4;                 %     number of fins
l_fin = 0.11;              % [m] length of fins
w_fin = 0.0065;            % [m] width of fins
C_drag = 0.40;             %     drag coefficient of rocket {total guess}

A_vehicle = pi*r_vehicle^2 + N_fin*l_fin*w_fin;
% ^ [m2] vehicle cross sectional area

% oxidizer tank
H_tank = 1.25;%1.25;   % [m] tank height
R_tank = 0.0762; % [m] tank radius

V_tank = pi()*H_tank*R_tank^2; % [m3] tank volume

% injector
N_jet = 50;%50;%12;       %     number of injector holes {44 for burnell}
r_jet = 0.75E-3;%0.00155; % [m] injector hole radius

A_inj = N_jet*pi*r_jet^2; % [m2] injector flow area
% ^ set manually for other injector types (not ideal but its what we have)

Cd_inj = 0.5; % injector discharge coefficient {complete estimate, refine
              %                                 using cold flow}
Ccv_inj = 1;  % assumes no vena contracta formation outside injector due to
              % chamfering the inlet
C_bur = 0.15; % Burnell coefficient, 0.15 as suggested by Waxman in his
              % hybrid injector paper

% chamber and fuel grain dimensions
A_port = pi*r_port^2;     % [m2] initial fuel grain port area
L_fuel = 36*0.0254;       % [m]  fuel grain length
t_tube = 0.08*0.0254;     % [m]  thickness of fuel cardboard tube
r_chamber = 4.5*0.0254/2; % [m]  outer radius of combustion chamber
t_chamber = 0.125*0.0254; % [m]  thickness of chamber wall and fuel liner

r_grain = r_chamber - t_chamber; % [m] outer radius of fuel grain
r_fuel_end = r_grain - t_tube;   % [m] outer radius of fuel

% nozzle geometry
r_t = 1.37*0.0127;%1.42258*0.0127;% % [m] throat radius
r_e = r_t*2.2;%3.2441*0.0127;%  % [m] exit radius {original gave exit pressure
                  %                  much below atmospheric}
r_f = 2/1000;     % [m] throat fillet radius

A_t = pi*r_t^2; % [m2] throat area
A_e = pi*r_e^2; % [m2] exit area

alpha = 20*pi/180;         % deg to rad, cone half angle
lambda = (1+cos(alpha))/2; % correction coefficient (cosine losses)
length_nozzle = ((r_e)-(r_t)-(r_f)*(cos(alpha)-1))/(tan(alpha)); % [m]
length_throat = r_f*sin(alpha);                                  % [m]
length = length_nozzle+length_throat;                            % [m]

% launch rail
l_rail = 33*0.3048; % [m] length of the launch rail

%% Propellant Properties
% fuel properties
rho_f = 920;         % [kg/m3]  fuel density

% oxidizer properties
oxProps = readmatrix('N2O_fluidProperties_SI.csv'); % Couch et al
M_ox = 44.013;    % [kg/kmol] molar mass of N2O
P_super = 0;      % [Pa]      supercharging pressure {look into tank
                  %                                 pressurization systems}
                  % not well implemented yet

P_cr = 72.38e5; % [Pa] critical pressure of N2O
T_cr = 309.56;  % [K]  critical temperature of N2O
% ^ Ohgaki, via NIST

omega_ox = 0.142; % acentric factor of N2O, from Matheson Gas Databook
kappa_ox = 0.37464+1.54226*omega_ox-0.26992*omega_ox^2;

R_ox = R_univ/M_ox; % [J/kgK] oxidizer specific gas constant

% Peng-Robinson coefficients
a_ox = 0.45724*R_ox^2*T_cr^2/P_cr;
b_ox = 0.07780*R_ox*T_cr/P_cr;

B1 = -209.559; % specific enthalpy of N2O as an ideal gas [J/kg]
B2 =  61.3277; % coefficients
B3 = -52.5969; % valid for temperature range [183 K - 1000 K]
B4 =  249.352; % from IHS ESDU
B5 = -38.4368;

D1 = 0.2934e5 ; % molar heat capacity of N2O gas at constant pressure
D2 = 0.3236e5 ; % [J/(kmol*K)] coefficients
D3 = 1.1238e3 ; % valid for Temp range [100 K - 1500 K]
D4 = 0.2177e5 ; % from Perry's Chemical Engineers' Handbook
D5 = 479.4 ;

% ballistic coefficients
fuel_alpha = 0.0000876; % [m3/kg] estimate based on Liu et al
fuel_n = 0.3953;        %         estimate based on Liu et al

%% Vehicle Initial Conditions
alt = elev;     % [m]     initial altitude is launch site elevation
vel = 0;        % [m/s]   initial speed is zero {we aren't spinlaunch}

r_fuel_burn = r_port;                     % [m] set current fuel radius to
                                          %     initial radius

% first iteration needs initial chamber pressure estimate
P_c = 2.5E6; % [Pa]
% on all subsequent iterations, initial estimate will be pressure at
% previous iteration

% this just needs to be anything other than 1
x_tank = 0;
vap = false;

% for finding the maximum thrust and the time it happens
% there's probably a more clever way to implement this
thrust_max = 0;
time_thmax = 0;

% convert everything to Kelvin
T_amb = T_amb + 273.15;
T_tank = T_tank + 273.15;

% warning flags
pres_warn = false;
temp_warn = false;
volu_warn = false;

%% Evolution Over Time
rDot_ev = [];          % [mm/s]   fuel regression rate
rem_fuel_ev = [];      % [m]      remaining fuel
rem_fuel_mass_ev = []; % [kg]     remaining mass of fuel
m_ox_ev = [];          % [kg]     remaining mass of oxidizer
x_tank_ev = [];        %          vapour quality in tank
mDot_ox_ev = [];       % [kg/s]   oxidizer mass flow rate
G_ox_ev = [];          % [kg/m2s] oxidizer mass flux through fuel port
mDot_f_ev = [];        % [kg/s]   fuel mass flow rate
mDot_t_ev = [];        % [kg/s]   total propellant mass flow rate
OF_ev = [];            %          oxidizer/fuel mass ratio

T_tank_ev = [];        % [K]      oxidizer tank temperature
T_c_ev = [];           % [K]      combustion chamber temperature
P_tank_ev = [];        % [Pa]     oxidizer tank pressure
P_c_ev = [];           % [Pa]     combustion chamber pressure
P_e_ev = [];           % [Pa]     exhaust pressure
P_air_ev = [];         % [Pa]     atmospheric pressure

thrust_ev = [];        % [N]      thrust
Isp_ev = [];           % [s]      specific impulse
twr_ev = [];           %          thrust to weight ratio
D_ev = [];             % [N]      drag

alt_ev = [];           % [m]      vehicle altitude
vel_ev = [];           % [m/s]    vehicle speed
Mach_ev = [];          % [Mach]   vehicle speed
accel_ev = [];         % [m/s2]   vehicle acceleration
g_force = [];          % [g]      vehicle acceleration

slip_ev = [];

% time
timeStep = 0.01; % [s]
t_stop = 200;    % run for 200s. if it goes this long something is wrong
step_stop = t_stop/timeStep;

%% burn loop
for y = 0:step_stop
    sav = y;

    % set times
    rDot_ev(y+1,1) = y*timeStep;
    rem_fuel_ev(y+1,1) = y*timeStep;
    rem_fuel_mass_ev(y+1,1) = y*timeStep;
    m_ox_ev(y+1,1) = y*timeStep;
    x_tank_ev(y+1,1) = y*timeStep;
    mDot_ox_ev(y+1,1) = y*timeStep;
    G_ox_ev(y+1,1) = y*timeStep;
    mDot_f_ev(y+1,1) = y*timeStep;
    mDot_t_ev(y+1,1) = y*timeStep;
    OF_ev(y+1,1) = y*timeStep;

    T_tank_ev(y+1,1) = y*timeStep;
    T_c_ev(y+1,1) = y*timeStep;
    P_tank_ev(y+1,1) = y*timeStep;
    P_c_ev(y+1,1) = y*timeStep;
    P_e_ev(y+1,1) = y*timeStep;
    P_air_ev(y+1,1) = y*timeStep;

    thrust_ev(y+1,1) = y*timeStep;
    Isp_ev(y+1,1) = y*timeStep;
    twr_ev(y+1,1) = y*timeStep;
    D_ev(y+1,1) = y*timeStep;

    alt_ev(y+1,1) = y*timeStep;
    vel_ev(y+1,1) = y*timeStep;
    Mach_ev(y+1,1) = y*timeStep;
    accel_ev(y+1,1) = y*timeStep;
    g_force(y+1,1) = y*timeStep;

    %% tank pressure
    T_tank_ev(y+1,2) = T_tank;
    if x_tank == 1
        P_tank_ev(y+1,2) = P_tank;
    else
        P_vap = interpolate(T_tank_ev(y+1,2),1,2,oxProps);
        P_tank_ev(y+1,2) = P_vap+P_super;
    end

    % rho_ox = interpolate(T_tank_ev(y+1,2),1,3,oxProps);

    %% mass flows, combustion chamber conditions, and exhaust properties
    if mode == 2
        % cold flow, only need to calculate mass flow
        if x_tank == 1
            rho_ox = m_ox_ev(y,2)/V_tank;
            % ^ assumes incompressible flow between tank and injector

            % assumes perfect gas, from Waxman
            if P_amb/P_tank_ev(y+1,2) <= (2/(k_ox+1))^(k_ox/(k_ox-1))
                % critical flow
                % assumes perfect gas, from Waxman
                P_cr = P_tank_ev(y+1,2)*(2/(k_ox+1))^(k_ox/(k_ox-1));
                mDot_ox_ev(y+1,2) = A_inj*Ccv_inj*Cd_inj*rho_ox*sqrt(2*(c_v+R_ox)*T_tank_ev(y+1,2)*((P_cr/P_tank_ev(y+1,2))^(2/k_ox)-(P_cr/P_tank_ev(y+1,2))^((k_ox+1)/k_ox))) * (1 + mDot_ox_noise*randn);
            else
                % subcritical flow
                % assumes perfect gas, from Waxman
                mDot_ox_ev(y+1,2) = A_inj*Ccv_inj*Cd_inj*rho_ox*sqrt(2*(c_v+R_ox)*T_tank_ev(y+1,2)*((P_amb/P_tank_ev(y+1,2))^(2/k_ox)-(P_amb/P_tank_ev(y+1,2))^((k_ox+1)/k_ox))) * (1 + mDot_ox_noise*randn);
            end
        else
            % two phase mass flow
            % from Burnell via Waxman

            % assumes choking only depends on vapour phase
            k_ox = interpolate_ratio(T_tank_ev(y+1,2),1,21,20,oxProps);

            if P_amb/P_tank_ev(y+1,2) <= (2/(k_ox+1))^(k_ox/(k_ox-1))
                % critical flow
                % assumes converging-only oriface, Ma=1 at injector
                % outlet
                P_inj_e = P_tank_ev(y+1,2)*(2/(k_ox+1))^(k_ox/(k_ox-1));
            else
                % subcritical flow
                % injector outlet pressure equals back pressure
                P_inj_e = P_amb;
            end

            v_2l = interpolate(P_inj_e,2,4,oxProps);
            v_2g = interpolate(P_inj_e,2,16,oxProps);
            v_2ratio = interpolate_ratio(P_inj_e,2,4,16,oxProps);
            h_ox = interpolate(T_tank_ev(y+1,2),1,6,oxProps);
            h_2l = interpolate(P_inj_e,2,6,oxProps);
            h_2g = interpolate(P_inj_e,2,18,oxProps);

            s_ox = interpolate(T_tank_ev(y+1,2),1,7,oxProps);
            s_2l = interpolate(P_inj_e,2,7,oxProps);
            s_2g = interpolate(P_inj_e,2,19,oxProps);

            x2 = (s_ox-s_2l)/(s_2g-s_2l); % ox quality downstream of injector

            % below not used ***
            babitskiy = @(x) sqrt(2.*(h_ox-(x.*h_2g+(1-x).*h_2l)))-sqrt(2.*(x.*v_2g+(1-x).*v_2l).*(P_tank_ev(y+1,2)-P_amb));
            %babitskiy = @(x) h_ox-(x.*h_2g+(1-x).*h_2l)-(x.*v_2g+(1-x).*v_2l).*(P_tank_ev(y+1,2)-P_c);

            for j = 0:0.01:1
                if babitskiy(j)*babitskiy(0)<0
                    %x2 = fzero(babitskiy,[0,j])*1.21;
                    break
                end
            end
            % above not used ***

            slip = v_2ratio^(-1/3);
            % 1/3 suggested by Moody and recommended by Waxman
            % 1/2 suggested by Fauske but refuted by Waxman

            % Moody (non-homogenous equilibrium)
            mDot_ox_ev(y+1,2) = (A_inj*Ccv_inj*Cd_inj*slip*sqrt(2*(h_ox-x2*h_2g-(1-x2)*h_2l)/(x2*(slip^2-1)+1))/(v_2g*(x2+slip*(1-x2)*v_2ratio))) * (1 + mDot_ox_noise*randn);

            % Burnell (homogeneous frozen non-equilibrium)
            %rho_ox = interpolate(T_tank_ev(y+1,2),1,3,oxProps);
            %mDot_ox_ev(y+1,2) = A_inj*Ccv_inj*Cd_inj*sqrt(2*rho_ox*(P_tank_ev(y+1,2)-P_vap*(1-C_bur))) * (1 + mDot_ox_noise*randn);
        end
    else
        % burn
        % all of this is interdependent and so must be solved iteratively
        P_c_residual = [];
        residueSize = 6; % number of previous iterations compared
        residuePrecision = 1e-6; % allowable fractional difference

        for z = 1:1000
            % loop until values converge
            % cap of 1000 to prevent infinite loop. It should converge well
            % before this (<30)

            %% oxidizer mass flux
            if x_tank == 1
                rho_ox = m_ox_ev(y,2)/V_tank;
                % ^ assumes incompressible flow between tank and injector

                if P_c/P_tank_ev(y+1,2) <= (2/(k_ox+1))^(k_ox/(k_ox-1))
                    % critical flow
                    % assumes perfect gas, from Waxman
                    P_cr = P_tank_ev(y+1,2)*(2/(k_ox+1))^(k_ox/(k_ox-1));
                    mDot_ox = A_inj*Ccv_inj*Cd_inj*rho_ox*sqrt(2*(c_v+R_ox)*T_tank_ev(y+1,2)*((P_cr/P_tank_ev(y+1,2))^(2/k_ox)-(P_cr/P_tank_ev(y+1,2))^((k_ox+1)/k_ox))) * (1 + mDot_ox_noise*randn);
                else
                    % subcritical flow
                    % assumes perfect gas, from Waxman
                    mDot_ox = A_inj*Ccv_inj*Cd_inj*rho_ox*sqrt(2*(c_v+R_ox)*T_tank_ev(y+1,2)*((P_c/P_tank_ev(y+1,2))^(2/k_ox)-(P_c/P_tank_ev(y+1,2))^((k_ox+1)/k_ox))) * (1 + mDot_ox_noise*randn);
                end
            else
                % two phase mass flow
                % from Moody via Waxman

                % assumes choking only depends on vapour phase
                k_ox = interpolate_ratio(T_tank_ev(y+1,2),1,21,20,oxProps);
                
                if P_c/P_tank_ev(y+1,2) <= (2/(k_ox+1))^(k_ox/(k_ox-1))
                    % critical flow
                    % assumes converging-only oriface, Ma=1 at injector
                    % outlet
                    P_inj_e = P_tank_ev(y+1,2)*(2/(k_ox+1))^(k_ox/(k_ox-1));
                else
                    % subcritical flow
                    % injector outlet pressure equals back pressure
                    P_inj_e = P_c;
                end

                v_2l = interpolate(P_inj_e,2,4,oxProps);
                v_2g = interpolate(P_inj_e,2,16,oxProps);
                v_2ratio = interpolate_ratio(P_inj_e,2,4,16,oxProps);
                h_ox = interpolate(T_tank_ev(y+1,2),1,6,oxProps);
                h_2l = interpolate(P_inj_e,2,6,oxProps);
                h_2g = interpolate(P_inj_e,2,18,oxProps);

                s_ox = interpolate(T_tank_ev(y+1,2),1,7,oxProps);
                s_2l = interpolate(P_inj_e,2,7,oxProps);
                s_2g = interpolate(P_inj_e,2,19,oxProps);

                x2 = (s_ox-s_2l)/(s_2g-s_2l);
                % ^ ox quality downstream of injector

                slip = v_2ratio^(-1/3);
                % 1/3 suggested by Moody and recommended by Waxman
                % 1/2 suggested by Fauske but refuted by Waxman

                % Moody
                mDot_ox = (A_inj*Ccv_inj*Cd_inj*slip*sqrt(2*(h_ox-x2*h_2g-(1-x2)*h_2l)/(x2*(slip^2-1)+1))/(v_2g*(x2+slip*(1-x2)*v_2ratio))) * (1 + mDot_ox_noise*randn);

                % Burnell (homogeneous frozen non-equilibrium)
                % rho_ox = interpolate(T_tank_ev(y+1,2),1,3,oxProps);
                % mDot_ox = A_inj*Ccv_inj*Cd_inj*sqrt(2*rho_ox*(P_tank_ev(y+1,2)-P_vap*(1-C_bur))) * (1 + mDot_ox_noise*randn);
            end

            G_ox = mDot_ox/(pi*r_fuel_burn^2);

            %% fuel regression and mass flow
            rDot = 1000*fuel_alpha*(G_ox)^fuel_n;

            dr_fuel_burn = rDot*timeStep/1000;

            mDot_f = rho_f*(pi*(r_fuel_burn+dr_fuel_burn)^2*L_fuel-pi*r_fuel_burn^2*L_fuel)/timeStep;

            mDot_t = mDot_ox+mDot_f;
            OF = mDot_t/mDot_f;

            %% chamber temperature
            % molar mass, ratio of specific heats, and chamber temperature
            % are based on curves fit to data from RPA

            % ideally we'd implement a Gibbs free energy solver but I don't
            % know how to do that, plus it would likely increase
            % computation time

            T_A = P_c*3.807e-4-5013;  % chamber temperature coefficients
            T_B = P_c*-2.579e-4+3550; % valid for pressure range 1.8 to
            T_C = P_c*6.188e-5-582.3; % 2.9MPa
            T_D = 27.17/(-0.5729+exp(6.696e-8*(P_c+1.644e-3)))-17.01;
            T_E = -0.7608/(-0.6346+exp(6.498e-8*(P_c+0.5652)))+0.7391;

            T_c = T_A+T_B*OF+T_C*OF^2+T_D*OF^3+T_E*OF^4;
            % ^ [K] combustion chamber temperature
            %       valid for mixture ratio range 5.25 to 9.5

            K_A = P_c*1.408e-8-0.2563; % heat ratio coefficients
            K_B = P_c*-3.301e-8+2.190; % valid for pressure range 1.8 to
            K_C = 1.439;               % 2.9MPa
            K_D = P_c*1.542e-7+5.454;
            K_E = P_c*-3.598e-9+1.280;

            k_c = K_A/(K_B+exp(-K_C*(OF-K_D)))+K_E;
            % ^ ratio of specific heats of exhaust gasses
            %   valid for mixture ratio range 5.25 to 9.5

            M_A = -6.869;              % molar mass coefficients
            M_B = -5.122/(8.424+exp(1.305e-6*(P_c-0.4316)))-0.06811;
            M_C = -0.4484/(1.088+exp(-6.519e-7*(P_c+2.219e-7)))+0.06493;
            M_D = P_c*5.095e-7+2.297;  % valid for pressure range 1.8 to
            M_E = P_c*-5.524e-8+28.82; % 2.9MPa

            M_c = M_A/(M_B+exp(-M_C*(OF-M_D)))+M_E;

            R_c = R_univ/(M_c); % exhaust specific gas constant

            mach_area = @(M) (1./M).*( (2/(k_c+1)) .* (1+(k_c-1)/2*M.^2) ) .^ ((k_c+1)/(2*(k_c-1))) - A_e/A_t;
            % ^ area ratio mach relationship (isentropic flow)

            for zz = 1:0.1:5
                if mach_area(1)*mach_area(zz)<0
                    Ma_e = fzero(mach_area,[1,zz]); % exhaust mach number
                    break
                end
            end

            P_c = mDot_t*sqrt(k_c*R_c*T_c)/(A_t*k_c*sqrt((2/(k_c+1))^((k_c+1)/(k_c-1))));

            P_ratio = (1+((k_c-1)*Ma_e^2)/2)^(-k_c/(k_c-1));
            % ^ ratio of exhaust pressure to chamber pressure

            T_e = T_c/((1/P_ratio)^((k_c-1)/(k_c)));
            % ^ [K] exhaust temperature

            v_e = Ma_e*sqrt(k_c*R_c*T_e); % [m/s] exhaust velocity

            % residuals
            % this will end the loop when the current calculated pressure
            % differs from the average of the previous [residueSize]
            % calculated pressures by less than [residuePrecision*100]
            % percent
            P_c_residual(z) = P_c;

            if z > residueSize
                rollingAvDiff = mean(P_c_residual(z-residueSize:z-1))/P_c;

                if (rollingAvDiff <= 1+residuePrecision) && (rollingAvDiff >= 1-residuePrecision)
                    break;
                end
            end

            if z==1000
                disp(["WARNING: chamber pressure failed to converge after 1000 iterations at time ",num2str(y*timeStep),"s"])
            end
        end

        r_fuel_burn = r_fuel_burn + dr_fuel_burn;

        mDot_ox_ev(y+1,2) = mDot_ox;
        G_ox_ev(y+1,2) = G_ox;
        rDot_ev(y+1,2) = rDot;
        rem_fuel_ev(y+1,2) = r_fuel_end-r_fuel_burn;
        rem_fuel_mass_ev(y+1,2) = pi*(r_fuel_end^2 - r_fuel_burn^2)*L_fuel*rho_f;
        mDot_f_ev(y+1,2) = mDot_f;
        mDot_t_ev(y+1,2) = mDot_t;
        OF_ev(y+1,2) = OF;

        T_c_ev(y+1,2) = T_c;
        P_c_ev(y+1,2) = P_c;
        P_e_ev(y+1,2) = P_c_ev(y+1,2)*P_ratio;
    end

    %% tank temperature
    v_SP_oxl = interpolate(T_tank_ev(y+1,2),1,4,oxProps);
    v_SP_oxg = interpolate(T_tank_ev(y+1,2),1,16,oxProps);
    c_v_oxl = interpolate(T_tank_ev(y+1,2),1,8,oxProps);
    c_v_oxg = interpolate(T_tank_ev(y+1,2),1,20,oxProps);
    c_p_T = (4.8 + 0.00322*T_tank_ev(y+1,2))*155.239;
    h_ox_l = interpolate(T_tank_ev(y+1,2),1,6,oxProps);
    h_ox_g = interpolate(T_tank_ev(y+1,2),1,18,oxProps);
    h_vap = h_ox_g - h_ox_l;

    if y>0
        m_ox_g_prev = m_ox_g;
    end

    m_ox = m_ox - mDot_ox_ev(y+1,2)*timeStep;
    m_ox_g = (V_tank-m_ox*v_SP_oxl)/(v_SP_oxg-v_SP_oxl);
    m_ox_l = m_ox-m_ox_g;

    if y==0
        m_ox_g_prev = m_ox_g;
    end

    if m_ox*v_SP_oxl>V_tank
        volu_warn = true;
    end

    x_tank = m_ox_g/m_ox;

    if x_tank >= 1 || vap
        % tank contents have fully vapourized
        x_tank = 1;

        % equations of state
        c_v_ox = @(T) (D1 + D2.*((D3./T)./sinh(D3./T)).^2 + D4.*((D5./T)./cosh(D5./T)).^2)/M_ox - R_ox;
        h_ox = @(T_r) (B1+B2.*T_r.^.5+B3.*T_r+B4.*T_r.^1.5+B5.*T_r.^2).*1000;

        V_1 = V_tank*m_ox/(m_ox + mDot_ox_ev(y+1,2)*timeStep);

        c_v = c_v_ox(T_tank);
        k_ox = (c_v+R_ox)/c_v; % assumes ideal gas

        % assumes adiabatic expansion of ideal gas
        T_tank = T_tank*(V_tank\V_1)^(k_ox-1);

        % Peng-Robinson gas equation
        P_2 = R_ox*T_tank/(V_tank/m_ox+b_ox)-a_ox*(1+kappa_ox*(1-sqrt(T_tank/T_cr)))^2/((V_tank/m_ox)^2+2*b_ox*V_tank/m_ox-b_ox^2);

        % fudged to make tank pressure continuous. This is to match the
        % behaviour in Waterloo's static fire
        if ~vap
            P_fudge = P_tank_ev(y+1,2) - P_2;
        end

        P_tank = P_2 + P_fudge;

        vap = true;
    else
        % tank contents are mixed
        % adapted from method used by Fernandez
        dm_ox_g = m_ox_g - m_ox_g_prev;
        dm_ox_l = - dm_ox_g - mDot_ox_ev(y+1,2);

        a = (m_tank*c_p_T+m_ox_l*c_v_oxl+m_ox_g*c_v_oxg); 
        % ^ very far off ideal gas model. not unexpected, ideal gas is very
        %   bad for NOS

        b = P_tank_ev(y+1,2)*v_SP_oxl;
        e = - h_vap + P_tank_ev(y+1,2)*v_SP_oxg;

        dT = (b*dm_ox_l+e*dm_ox_g)/a;
        T_tank = T_tank+dT*timeStep;

        Z_tank = P_tank_ev(y+1,2)*v_SP_oxg/(R_ox*T_tank);
    end

    m_ox_ev(y+1,2) = m_ox;
    x_tank_ev(y+1,2) = m_ox_g/m_ox;

    if mode == 2
        if m_ox <= 0.1
            disp('Burnout: out of oxidizer')
            break;
        end
        if P_tank_ev(y+1,2) <= P_amb
            disp('Burnout: tank pressure less than atmospheric pressure')
            disp(['Remaining oxidizer: ',num2str(m_ox), 'kg'])
            break;
        end

        % skip the rest of the loop (not relevant)
        continue;
    end

    %% thrust and drag
    T_air = T_amb - (alt-elev)*9.8/1000;
    % ^ [deg C] temperature variation with altitude, assuming dry air
    P_air = P_sea*(1-alt*2.25577E-5)^5.25588;
    % ^ [Pa]    air pressure at altitude (engineer's toolbox)
    rho_air = P_air/(R_air*T_air);
    % ^ [kg/m3] ideal gas density of air

    P_air_ev(y+1,2) = P_air;

    thrust_ev(y+1,2) = lambda*mDot_t_ev(y+1,2)*v_e + (P_e_ev(y+1,2)-P_air_ev(y+1,2))*A_e;
    Isp_ev(y+1,2) = thrust_ev(y+1,2)/(g*mDot_t_ev(y+1,2));

    if thrust_ev(y+1,2) > thrust_max
        thrust_max = thrust_ev(y+1,2);
        time_thmax = thrust_ev(y+1,1);
    end

    m_vehicle = m_dry + m_ox_ev(y+1,2) + rem_fuel_mass_ev(y+1,2); % [kg] vehicle mass

    twr_ev(y+1,2) = thrust_ev(y+1,2)/(g*m_vehicle);

    D_ev(y+1,2) = (C_drag*rho_air*vel^2*A_vehicle)/2; % drag

    %% flight kinematics
    if mode==0
        accel_ev(y+1,2) = (thrust_ev(y+1,2)-D_ev(y+1,2))/m_vehicle - g;
        g_force(y+1,2) = accel_ev(y+1,2)/g;

        vel = vel + accel_ev(y+1,2)*timeStep; % speed
        vel_ev(y+1,2) = vel;

        Mach_ev(y+1,2) = vel/sqrt(1.4*R_air*T_air);

        alt = alt + vel*timeStep; % altitude
        alt_ev(y+1,2) = alt;
        
        if y ~= 0
            if alt_ev(y+1,2) - elev > l_rail && alt_ev(y,2) - elev < l_rail
                disp(['Speed at end of launch rail: ',num2str(vel_ev(y,2)/0.3048),'ft/s'])
            end
        end
    end

    %% physical stops
    if r_fuel_burn>=r_fuel_end
        disp('Burnout: out of fuel')
        disp(['Remaining oxidizer: ',num2str(m_ox), 'kg'])
        break;
    end
    if P_tank_ev(y+1,2) <= P_c_ev(y+1,2)
        disp('Burnout: tank pressure less than chamber pressure')
        disp(['Remaining oxidizer: ',num2str(m_ox), 'kg'])
        disp(['Remaining fuel:     ',num2str(rem_fuel_mass_ev(y+1,2)), 'kg'])
        break;
    end
    if m_ox <= 0
        disp('Burnout: out of oxidizer')
        disp(['Remaining fuel: ',num2str(rem_fuel_mass_ev(y+1,2)), 'kg'])
        break;
    end
    if mDot_ox_ev(y+1,2) <= 0
        disp('Burnout: oxidizer flow stopped')
        disp(['Remaining oxidizer: ',num2str(m_ox), 'kg'])
        disp(['Remaining fuel:     ',num2str(rem_fuel_mass_ev(y+1,2)), 'kg'])
        break;
    end
    if T_c_ev(y+1,2) <= 2500
        disp('Cutoff: combustion poorly characterized, risk of instability')
        disp(['Remaining oxidizer: ',num2str(m_ox), 'kg'])
        disp(['Remaining fuel:     ',num2str(rem_fuel_mass_ev(y+1,2)), 'kg'])
        break;
    end

    if P_c_ev(y+1,2)>3E6
        pres_warn = true;
    end
    if T_c_ev(y+1,2)>3250
        temp_warn = true;
    end
end

if mode == 2
    disp(' ')
    disp(['Flow duration: ',num2str(sav*timeStep),'s'])
    disp(['Maximum flow rate: ',num2str(max(mDot_ox_ev(:,2))),'kg/s'])
else
    disp(' ')
    disp(['Burn duration: ',num2str(sav*timeStep),'s'])
    disp(['Impulse: ',num2str(trapz(thrust_ev(:,2))*timeStep),'Ns'])
    disp(['Maximum thrust: ',num2str(thrust_max),'N'])
    disp(['Average thrust: ',num2str(trapz(thrust_ev(:,2))/sav),'N'])
    disp(['Average Isp: ',num2str(trapz(Isp_ev(:,2))/sav),'s'])
end

%% coast loop
if mode==0
    for y = sav:step_stop
        accel_ev(y+1,1) = y*timeStep;
        vel_ev(y+1,1) = y*timeStep;
        alt_ev(y+1,1) = y*timeStep;
        Mach_ev(y+1,1) = y*timeStep;
        g_force(y+1,1) = y*timeStep;
        D_ev(y+1,1) = y*timeStep;

        T_air = T_amb - (alt-elev)*9.8/1000;
        % ^ [deg C] temperature variation with altitude, assuming dry air
        P_air = P_sea*(1-alt*2.25577E-5)^5.25588;
        % ^ [Pa]    air pressure at altitude (engineer's toolbox)
        rho_air = P_air/(R_air*T_air);
        % ^ [kg/m3] ideal gas density of air

        D_ev(y+1,2) = (C_drag*rho_air*vel^2*A_vehicle)/2; % drag

        accel_ev(y+1,2) = -D_ev(y+1,2)/m_vehicle - g; % acceleration
        g_force(y+1,2) = accel_ev(y+1,2)/g;

        vel = vel + accel_ev(y+1,2)*timeStep; % speed
        vel_ev(y+1,2) = vel;
        Mach_ev(y+1,2) = vel/sqrt(1.4*R_air*T_air);

        alt = alt + vel*timeStep; % altitude
        alt_ev(y+1,2) = alt;

        if vel <= 0
            % stop at apogee
            break;
        end
    end

% altitude in feet
alt_ev(:,3) = alt_ev(:,2)/0.3048;

end

% convert tank temperature to deg C for easier understanding
T_tank_ev(:,2) = T_tank_ev(:,2) - 273.15;

%% Plots
if mode == 2
    figure(1)
        plot(m_ox_ev(:,1),m_ox_ev(:,2),'b')
        title('Remaining Oxidizer vs Time')
        xlabel('Time [s]')
        ylabel('Mass Remaining [kg]');
    
    figure(2)
        plot(mDot_ox_ev(:,1),mDot_ox_ev(:,2),'b')
        title('Mass Flow Rate vs Time')
        xlabel('Time [s]')
        ylabel('Mass Flow Rate [kg/s]');
    
    figure(3)
        plot(P_tank_ev(:,1),P_tank_ev(:,2),'b')
        title('Tank Pressure vs Time')
        xlabel('Time [s]')
        ylabel('Pressure [Pa]');

    P_delta = P_tank_ev(:,2) - P_amb;
    figure(4)
        plot(P_tank_ev(:,1),P_delta,'m')
        title('Pressure Drop Across Injector vs Time')
        xlabel('Time [s]')
        ylabel('Pressure Drop [Pa]');

    figure(5)
        plot(T_tank_ev(:,1),T_tank_ev(:,2),'b')
        title('Tank Temperature vs Time')
        xlabel('Time [s]')
        ylabel('Temperature [C]');
else
    % propellant plots
    figure(1)
        plot(rDot_ev(:,1),rDot_ev(:,2))
        title('Fuel Regression Rate vs Time')
        xlabel('Time [s]')
        ylabel('Regression Rate [mm/s]');

    figure(2)
        plot(m_ox_ev(:,1),m_ox_ev(:,2),'b',rem_fuel_mass_ev(:,1),rem_fuel_mass_ev(:,2),'r')
        title('Remaining Propellant vs Time')
        xlabel('Time [s]')
        ylabel('Mass Remaining [kg]')
        legend('oxidizer','fuel');

    figure(3);
        plot(mDot_ox_ev(:,1),mDot_ox_ev(:,2),'b',mDot_f_ev(:,1),mDot_f_ev(:,2),'r',mDot_t_ev(:,1),mDot_t_ev(:,2),'m')
        title('Mass Flow Rate vs Time')
        xlabel('Time [s]')
        ylabel('Mass Flow Rate [kg/s]')
        legend('oxidizer mass flow rate','fuel mass flow rate','total mass flow rate');

    figure(4);
        plot(G_ox_ev(:,1),G_ox_ev(:,2),'b')
        title('Oxidizer Mass Flux Through Fuel Port vs Time')
        xlabel('Time [s]')
        ylabel('Mass Flux [kg/m2s]');

    figure(5)
        plot(OF_ev(:,1),OF_ev(:,2))
        title('Mixture Ratio vs Time')
        xlabel('Time [s]')
        ylabel('Mixture Ratio O/F');

    % pressure and temperature plots
    figure(6)
        plot(P_tank_ev(:,1),P_tank_ev(:,2),'b',P_c_ev(:,1),P_c_ev(:,2),'r')
        title('Pressure vs Time (tank and chamber)')
        xlabel('Time [s]')
        ylabel('Pressure [Pa]')
        legend('tank pressure','combustion chamber pressure');

    P_delta = P_tank_ev(:,2) - P_c_ev(:,2);
    figure(7)
        plot(P_tank_ev(:,1),P_delta,'m')
        title('Pressure Drop Across Injector vs Time')
        xlabel('Time [s]')
        ylabel('Pressure Drop [Pa]');

    figure(8)
        plot(P_e_ev(:,1),P_e_ev(:,2),'y',P_air_ev(:,1),P_air_ev(:,2),'c')
        title('Pressure vs Time (exhaust and atmosphere)')
        xlabel('Time [s]')
        ylabel('Pressure [Pa]')
        legend('exhaust pressure','air pressure');

    figure(9)
        plot(T_tank_ev(:,1),T_tank_ev(:,2),'b')
        title('Tank Temperature vs Time')
        xlabel('Time [s]')
        ylabel('Temperature [C]');

    figure(10)
        plot(T_c_ev(:,1),T_c_ev(:,2),'r')
        title('Combustion Chamber Temperature vs Time')
        xlabel('Time [s]')
        ylabel('Temperature [K]');

    % performance plots
    figure(11)
        plot(thrust_ev(:,1),thrust_ev(:,2))
        title('Thrust vs Time')
        xlabel('Time [s]')
        ylabel('Thrust [N]');

    figure(12)
        plot(Isp_ev(:,1),Isp_ev(:,2))
        title('Specific Impulse vs Time')
        xlabel('Time [s]')
        ylabel('Specific Impulse [s]');

    figure(13)
        plot(twr_ev(:,1),twr_ev(:,2))
        title('Thrust to Weight Ratio vs Time')
        xlabel('Time [s]')
        ylabel('Thrust to Weight Ratio');

    if mode==0
        figure(14)
            plot(D_ev(:,1),D_ev(:,2))
            title('Drag vs Time')
            xlabel('Time [s]')
            ylabel('Drag [N]');
    end

    % kinematic plots
    if mode==0
        figure(15)
            plot(alt_ev(:,1),alt_ev(:,2))
            title('Altitude vs Time')
            xlabel('Time [s]')
            ylabel('Altitude [m]');

        figure(16)
            plot(alt_ev(:,1),alt_ev(:,3))
            title('Altitude vs Time')
            xlabel('Time [s]')
            ylabel('Altitude [ft]');

        figure(17)
            plot(vel_ev(:,1),vel_ev(:,2))
            title('Velocity vs Time')
            xlabel('Time [s]')
            ylabel('Velocity [m/s]');

        figure(18)
            plot(Mach_ev(:,1),Mach_ev(:,2))
            title('Velocity vs Time')
            xlabel('Time [s]')
            ylabel('Velocity [Mach]');

        figure(19)
            plot(accel_ev(:,1),accel_ev(:,2))
            title('Acceleration vs Time')
            xlabel('Time [s]')
            ylabel('Acceleration [m/s2]');
    
        figure(20)
            plot(g_force(:,1),g_force(:,2))
            title('Acceleration vs Time')
            xlabel('Time [s]')
            ylabel('Acceleration [g]');
    end

    if ras
        ds_length = height(thrust_ev(:,1)); % length of the time matrix
        ds_freq = ds_length/31;             % sampling frequency required
        warning('off','MATLAB:colon:nonIntegerIndex')
        % ^ suppress colon opperator integer warning
        ThrustCurveRA = thrust_ev(2:ds_freq:end,:);
        % ^ thrust curve with 31 elements needed for rasp.eng file format

        % final point
        ThrustCurveRA(32,2) = 0;
        ThrustCurveRA(32,1) = ds_length*timeStep;

        for i = 1:31
            if time_thmax < ThrustCurveRA(i+1,1)
                % specify maximum thrust
                ThrustCurveRA(i,1) = time_thmax;
                ThrustCurveRA(i,2) = thrust_max;
                break
            end
        end

        disp(' ')
        disp('RAS input:')

        for i = 1:32
            disp([num2str(round(ThrustCurveRA(i,1),2)),' ',num2str(round(ThrustCurveRA(i,2),2))])
        end
    end
end

% warnings
if pres_warn
    disp('DANGER: chamber pressure exceeded 3MPa')
end
if temp_warn
    disp('DANGER: chamber temperature exceeded 3250K')
end
if volu_warn
    disp('WARNING: loaded oxidizer volume greater than tank volume')
end

%% Functions
% interpolate a value on a given table, given known value A, horizontal
% index of known value i, and horizontal index of desired value j
function [B] = interpolate (A, i, j, table)
    B = interp1(table(:,i),table(:,j),A,'spline');
end

function [B] = interpolate_ratio (A, i, j, k, table)
    B = interp1(table(:,i),table(:,j)./table(:,k),A,'spline');
end