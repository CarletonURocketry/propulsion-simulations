% CU InSpace Hybrid Rocket Engineering Model
% Academic Year 2022-2023
% Last update Nov 11, 2022
% Quentin Alexander
% quentinalexander@cmail.carleton.ca
% 
% This script simulates the flight of a N2O-paraffin hybrid rocket given
% launch site conditions and vehicle mass and geometry
% 
% Returns propellant burn rate, pressure and temperature over time in the
% tank and combustion chamber, thrust, Isp, and thrust to weight ratio, and
% altitude, velocity (in m/s and Mach), and acceleration
% 
% Warning: the simulation will fail if ambient temperature is higher than
% 36 degrees C (critical temperature of N2O)

clear all; %#ok<CLALL> 
close all;

mDot_ox_noise = 0.00; % random noise in oxidizer mass flow, approximately
                      % as a fraction of mass flow
flight = true;        % true - flight ; false - static fire
ras = false;          % are we putting this through RAS

%% Constants
g = 9.81;                  % [m/s2]    acceleration due to gravity
R_univ = 8314.46261815324; % [J/kmolK] universal gas constant
R_air = R_univ/28.97;      % [J/kgK]   specific gas constant of air

%% Launch Site Conditions
T_amb = 36;    % [deg C] ambient temperature
P_amb = 85E3;  % [Pa]    ambient pressure
elev = 1400;   % [m]     elevation of launch site above sea level

P_sea = P_amb/(1-elev*2.25577E-5)^5.25588; % [Pa] pressure at sea level
                                           %      formula from engineer's
                                           %      toolbox

%% Vehicle Mass
m_ox = 13;                 % [kg] mass of oxidizer loaded in the tank
m_tank = 3.2;              % [kg] mass of the tank
m_dry = 102.74*0.45359237; % [kg] dry mass of the rocket (includes tank)

%% Geometry
% aerostructure
r_vehicle = 6.25*0.0254/2; % [m] radius of vehicle
N_fin = 4;                 %     number of fins
l_fin = 0.17;              % [m] length of fins
w_fin = 0.0065;            % [m] width of fins
C_drag = 0.40;             %     drag coefficient of rocket {total guess}

% oxidizer tank
H_tank = 1;   % [m] tank height
R_tank = 0.0762; % [m] tank radius

V_tank = pi()*H_tank*R_tank^2; % [m3] tank volume

% injector
N_jet = 50;%12;      %     number of injector holes
r_jet = 0.75E-3;%0.00155; % [m] injector hole radius

A_jet = pi*r_jet^2; % [m2] injector hole area

Cd_inj = 0.5; % injector discharge coefficient {complete estimate, refine
              %                                 using cold flow}
Ccv_inj = 1;  % assumes no vena contracta formation outside injector due to
              % chamfering the inlet
C_bur = 0.15; % Burnell coefficient, 0.15 as suggested by Waxman in his
              % hybrid injector paper

% chamber and fuel grain dimensions
r_port = 2.375*25.4/2000; % [m]  initial fuel cylinder port radius
A_port = pi*r_port^2;     % [m2] initial fuel cylinder port area
L_fuel = 1.35;            % [m]  fuel grain lenght
d_casing = 0.098;         % [m]  outer diameter of engine
t_casing = 0.0025;        % [m]  thickness of engine casing
t_liner = 0.005;          % [m]  thickness of fuel liner

% nozzle geometry
r_t = 34.48/2000; % [m] throat radius
r_e = 71/2000;    % [m] exit radius {original gave exit pressure
                  %                  much below atmospheric}
r_f = 2/1000;     % [m] throat fillet radius

A_t = pi*r_t^2; % [m2] throat area
A_e = pi*r_e^2; % [m2] exit area

alpha = 20*pi/180;         % deg to rad, cone half angle
lambda = (1+cos(alpha))/2; % correction coefficient (cosine losses)
length_nozzle = ((r_e)-(r_t)-(r_f)*(cos(alpha)-1))/(tan(alpha)); % [m]
length_throat = r_f*sin(alpha);                                  % [m]
length = length_nozzle+length_throat;                            % [m]

%% Propellant Properties
% fuel properties
rho_f = 920;         % [kg/m3]  fuel density
M_fmon = 0.01402668; % [kg/mol] molar mass of fuel monomer {I hope using
                     %                                      this is valid}

% oxidizer properties
props = readmatrix('N2O_fluidProperties_SI.csv'); % Couch et al
M_ox = 0.044013; % [kg/mol] molar mass of N2O
P_super = 0;     % [Pa]     supercharging pressure {look into tank
                 %                                  pressurization systems}

% ballistic coefficients
fuel_alpha = 0.0000876; % [m3/kg] estimate based on Liu et al
fuel_n = 0.3953;      %         estimate based on Liu et al

%% Vehicle Initial Conditions
alt = elev;     % [m]     initial altitude is launch site elevation
vel = 0;        % [m/s]   initial speed is zero
T_tank = T_amb; % [deg C] tank temperature, assumed to start at ambient
                %         temperature

r_fuel_burn = r_port;                     % [m] set current fuel radius to
                                          %     initial radius
r_fuel_end = d_casing/2-t_casing-t_liner; % [m] radius of the combustion
                                          %     chamber once all the fuel
                                          %     is burnt

%% Evolution Over Time
rDot_ev = [];       % [mm/s] fuel regression rate
rem_fuel = [];      % [m]    remaining fuel
rem_fuel_mass = []; % [kg]   remaining mass of fuel
m_ox_ev = [];       % [kg]   remaining mass of oxidizer
mDot_ox_ev = [];    % [kg/s] oxidizer mass flow rate
mDot_f_ev = [];     % [kg/s] fuel mass flow rate
mDot_t_ev = [];     % [kg/s] total propellant mass flow rate

T_tank_ev = [];     % [K]    oxidizer tank temperature
T_c_ev = [];        % [K]    combustion chamber temperature
P_tank_ev = [];     % [Pa]   oxidizer tank pressure
P_c_ev = [];        % [Pa]   combustion chamber pressure
P_e_ev = [];        % [Pa]   exhaust pressure
P_air_ev = [];      % [Pa]   atmospheric pressure

thrust_ev = [];     % [N]    thrust
Isp_ev = [];        % [s]    specific impulse
twr_ev = [];        %        thrust to weight ratio
D_ev = [];          % [N]    drag

alt_ev = [];        % [m]    vehicle altitude
vel_ev = [];        % [m/s]  vehicle speed
Mach_ev = [];       % [Mach] vehicle speed
accel_ev = [];      % [m/s2] vehicle acceleration
g_force = [];       % [g]    vehicle acceleration

% time
timeStep = 0.01; % [s]
t_stop = 100;    % run for 100s. if it goes this long something is wrong
step_stop = t_stop/timeStep;

%% burn loop
for y = 0:step_stop
    % set times
    T_tank_ev(y+1,1) = y*timeStep;
    P_tank_ev(y+1,1) = y*timeStep;
    mDot_ox_ev(y+1,1) = y*timeStep;
    m_ox_ev(y+1,1) = y*timeStep;
    rDot_ev(y+1,1) = y*timeStep;
    rem_fuel(y+1,1) = y*timeStep;
    rem_fuel_mass(y+1,1) = y*timeStep;
    mDot_f_ev(y+1,1) = y*timeStep;
    mDot_t_ev(y+1,1) = y*timeStep;
    P_c_ev(y+1,1) = y*timeStep;
    P_e_ev(y+1,1) = y*timeStep;
    T_c_ev(y+1,1) = y*timeStep;
    thrust_ev(y+1,1) = y*timeStep;
    Isp_ev(y+1,1) = y*timeStep;
    twr_ev(y+1,1) = y*timeStep;
    P_air_ev(y+1,1) = y*timeStep;
    accel_ev(y+1,1) = y*timeStep;
    g_force(y+1,1) = y*timeStep;
    vel_ev(y+1,1) = y*timeStep;
    alt_ev(y+1,1) = y*timeStep;
    Mach_ev(y+1,1) = y*timeStep;
    D_ev(y+1,1) = y*timeStep;

    %% tank pressure
    T_tank_ev(y+1,2) = T_tank;
    P_vap = oxProp(T_tank_ev(y+1,2),2,props);
    P_tank_ev(y+1,2) = P_vap+P_super;
    rho_ox = oxProp(T_tank_ev(y+1,2),3,props);

    %% oxidizer mass flux
    v_jet = (1/rho_ox)*Ccv_inj*Cd_inj*sqrt(2*rho_ox*(P_tank_ev(y+1,2) - P_vap*(1-C_bur)));
    mDot_ox_ev(y+1,2) = rho_ox*N_jet*A_jet*v_jet + rho_ox*N_jet*A_jet*v_jet*mDot_ox_noise*randn;

    G_ox = mDot_ox_ev(y+1,2)/(pi*r_fuel_burn^2);
    
    %% fuel regression and mass flow
    rDot_ev(y+1,2) = 1000*fuel_alpha*(G_ox)^fuel_n;
    rem_fuel(y+1,2) = r_fuel_end-r_fuel_burn;
    rem_fuel_mass(y+1,2) = pi*(r_fuel_end^2 - r_fuel_burn^2)*L_fuel*rho_f;

    r_fuel_burn = r_fuel_burn+rDot_ev(y+1,2)*timeStep/1000;

    mDot_f_ev(y+1,2) = rho_f*(rDot_ev(y+1,2)/1000)*2*pi*r_fuel_burn*L_fuel;

    mDot_t_ev(y+1,2) = mDot_ox_ev(y+1,2)+mDot_f_ev(y+1,2);

    %% tank temperature
    v_SP_oxl = oxProp(T_tank_ev(y+1,2),4,props);
    v_SP_oxg = oxProp(T_tank_ev(y+1,2),16,props);
    c_v_oxl = oxProp(T_tank_ev(y+1,2),8,props);
    c_v_oxg = oxProp(T_tank_ev(y+1,2),20,props);
    c_p_T = (4.8 + 0.00322*(T_tank_ev(y+1,2)+273.15))*155.239;
    h_ox_l = oxProp(T_tank_ev(y+1,2),6,props);
    h_ox_g = oxProp(T_tank_ev(y+1,2),18,props);
    h_vap = h_ox_g - h_ox_l;

    m_ox = m_ox - mDot_ox_ev(y+1,2)*timeStep;
    m_ox_g = (V_tank-m_ox*v_SP_oxl)/(v_SP_oxg-v_SP_oxl);
    m_ox_l = m_ox-m_ox_g;

    m_ox_ev(y+1,2) = m_ox;

    dm_ox_g = mDot_ox_ev(y+1,2)*v_SP_oxl/(v_SP_oxg-v_SP_oxl);
    dm_ox_l = - dm_ox_g - mDot_ox_ev(y+1,2);

    a = (m_tank*c_p_T+m_ox_l*c_v_oxl+m_ox_g*c_v_oxg); 
    % ^ very far off ideal gas model. not unexpected, ideal gas is very bad
    %   for NOS

    b = P_tank_ev(y+1,2)*v_SP_oxl;
    e = - h_vap + P_tank_ev(y+1,2)*v_SP_oxg;

    dT = (b*dm_ox_l+e*dm_ox_g)/a;
    T_tank = T_tank+dT*timeStep;

    %% chamber and exhaust pressure
    M_c = mDot_t_ev(y+1,2)/(mDot_ox_ev(y+1,2)/M_ox + 2*mDot_f_ev(y+1,2)/M_fmon);
    % ^ kg/mol, approximation of exhaust molar mass
    R_c = R_univ/(M_c*1000); % exhaust specific gas constant
    k_c = 1.23; % sp. heat ratio in the chamber, from Brendan
    % k_c is assumed to stay approximately constant

    fcn = @(M) (1./M).*( (2/(k_c+1)) .* (1+(k_c-1)/2*M.^2) ) .^ ((k_c+1)/(2*(k_c-1))) - A_e/A_t;
    % ^ area ratio mach relationship (isentropic flow). redefined each loop
    %   in case I do get an expression for k_c

    Ma_e = fzero(fcn,3); % exhaust mach number

    P_ratio = (1+((k_c-1)*Ma_e^2)/2)^(-k_c/(k_c-1));
    % ^ ratio of exhaust pressure to chamber pressure

    C_F = sqrt(k_c*(2/(k_c+1))^((k_c+1)/(k_c-1)))*sqrt((2*k_c/(k_c-1))*(1-P_ratio^((k_c-1)/k_c))) + A_e*P_ratio/A_t;
    % ^ coefficient of thrust (Zolla et al)

    T_c_ev(y+1,2) = 2750; % calculation needed (gibbs free energy?)

    T_e = T_c_ev(y+1,2)/((1/P_ratio)^((k_c-1)/(k_c)));
    % ^ [K] exhaust temperature

    v_e = Ma_e*sqrt(k_c*R_c*T_e); % [m/s] exhaust velocity

    cstar = v_e/C_F; % [m/s] characteristic velocity

    P_c_ev(y+1,2) = mDot_t_ev(y+1,2)*cstar/A_t;
    P_e_ev(y+1,2) = P_c_ev(y+1,2)*P_ratio;

    %% thrust and drag
    T_air = T_amb - (alt-elev)*9.8/1000;
    % ^ [deg C] temperature variation with altitude, assuming dry air
    P_air = P_sea*(1-alt*2.25577E-5)^5.25588;
    % ^ [Pa]    air pressure at altitude (engineer's toolbox)
    rho_air = P_air/(R_air*(T_air+273.15));
    % ^ [kg/m3] ideal gas density of air

    P_air_ev(y+1,2) = P_air;

    thrust_ev(y+1,2) = lambda*mDot_t_ev(y+1,2)*v_e + (P_e_ev(y+1,2)-P_air_ev(y+1,2))*A_e;
    Isp_ev(y+1,2) = thrust_ev(y+1,2)/(g*mDot_t_ev(y+1,2));

    m_vehicle = m_dry + m_ox_ev(y+1,2) + rem_fuel_mass(y+1,2); % [kg] vehicle mass

    twr_ev(y+1,2) = thrust_ev(y+1,2)/(g*m_vehicle);

    D_ev(y+1,2) = (C_drag*rho_air*vel^2*(pi*r_vehicle^2 + N_fin*l_fin*w_fin))/2; % drag

    %% flight kinematics
    if flight
        accel_ev(y+1,2) = (thrust_ev(y+1,2)-D_ev(y+1,2))/m_vehicle - g;
        g_force(y+1,2) = accel_ev(y+1,2)/g;

        vel = vel + accel_ev(y+1,2)*timeStep; % speed
        vel_ev(y+1,2) = vel;

        Mach_ev(y+1,2) = vel/sqrt(1.4*R_air*(T_air+273.15));

        alt = alt + vel*timeStep; % altitude
        alt_ev(y+1,2) = alt;
    end

    %% physical stops
    if r_fuel_burn>=r_fuel_end
        disp('burnout: out of fuel')
        disp(['remaining oxidizer: ',num2str(m_ox), 'kg'])
        sav = y;
        break;
    end
    if P_tank_ev(y+1,2) <= P_c_ev(y+1,2)
        disp('burnout: tank pressure less than chamber pressure')
        disp(['remaining oxidizer: ',num2str(m_ox), 'kg'])
        disp(['remaining fuel:     ',num2str(rem_fuel(y+1,2)*1000), 'mm'])
        sav = y;
        break;
    end
    if m_ox <= 0
        disp('burnout: out of oxidizer')
        disp(['remaining fuel: ',num2str(rem_fuel(y+1,2)*1000), 'mm'])
        sav = y;
        break;
    end

    if P_c_ev(y+1,2)>2.75E6
        disp('warning: chamber pressure has exceeded design constraint')
    end
    if T_c_ev(y+1,2)>2750
        disp('warning: chamber temperature has exceeded design constraint')
    end
end

%% coast loop
if flight
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
        rho_air = P_air/(R_air*(T_air+273.15));
        % ^ [kg/m3] ideal gas density of air

        D_ev(y+1,2) = (C_drag*rho_air*vel^2*(pi*r_vehicle^2 + N_fin*l_fin*w_fin))/2; % drag

        accel_ev(y+1,2) = -D_ev(y+1,2)/m_vehicle - g; % acceleration
        g_force(y+1,2) = accel_ev(y+1,2)/g;

        vel = vel + accel_ev(y+1,2)*timeStep; % speed
        vel_ev(y+1,2) = vel;
        Mach_ev(y+1,2) = vel/sqrt(1.4*R_air*(T_air+273.15));

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

% convert tank temperature to K
%T_tank_ev(:,2) = T_tank_ev(:,2) + 273.15;

disp(['Impulse: ',num2str(trapz(thrust_ev(:,2))*timeStep),'Ns'])

%% Plots
% propellant plots
figure(1)
    plot(rDot_ev(:,1),rDot_ev(:,2))
    title('Fuel Regression Rate vs Time')
    xlabel('Time [s]')
    ylabel('Regression Rate [mm/s]');

figure(2)
    plot(m_ox_ev(:,1),m_ox_ev(:,2),'b',rem_fuel_mass(:,1),rem_fuel_mass(:,2),'r')
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

% pressure and temperature plots
figure(4)
    plot(P_tank_ev(:,1),P_tank_ev(:,2),'b',P_c_ev(:,1),P_c_ev(:,2),'r')
    title('Pressure vs Time (tank and chamber)')
    xlabel('Time [s]')
    ylabel('Pressure [Pa]')
    legend('tank pressure','combustion chamber pressure');

figure(5)
    plot(P_e_ev(:,1),P_e_ev(:,2),'y',P_air_ev(:,1),P_air_ev(:,2),'c')
    title('Pressure vs Time (exhaust and atmosphere)')
    xlabel('Time [s]')
    ylabel('Pressure [Pa]')
    legend('exhaust pressure','air pressure');

figure(6)
    plot(T_tank_ev(:,1),T_tank_ev(:,2),'b')
    title('Tank Temperature vs Time')
    xlabel('Time [s]')
    ylabel('Temperature [C]');

figure(7)
    plot(T_c_ev(:,1),T_c_ev(:,2),'r')
    title('Combustion Chamber Temperature vs Time')
    xlabel('Time [s]')
    ylabel('Temperature [K]');

% performance plots
figure(8)
    plot(thrust_ev(:,1),thrust_ev(:,2))
    title('Thrust vs Time')
    xlabel('Time [s]')
    ylabel('Thrust [N]');

figure(9)
    plot(Isp_ev(:,1),Isp_ev(:,2))
    title('Specific Impulse vs Time')
    xlabel('Time [s]')
    ylabel('Specific Impulse [s]');

figure(10)
    plot(twr_ev(:,1),twr_ev(:,2))
    title('Thrust to Weight Ratio vs Time')
    xlabel('Time [s]')
    ylabel('Thrust to Weight Ratio');

if flight
    figure(11)
        plot(D_ev(:,1),D_ev(:,2))
        title('Drag vs Time')
        xlabel('Time [s]')
        ylabel('Drag [N]');
end

% kinematic plots
if flight
    figure(12)
        plot(alt_ev(:,1),alt_ev(:,2))
        title('Altitude vs Time')
        xlabel('Time [s]')
        ylabel('Altitude [m]');

    figure(13)
        plot(alt_ev(:,1),alt_ev(:,3))
        title('Altitude vs Time')
        xlabel('Time [s]')
        ylabel('Altitude [ft]');

    figure(14)
        plot(vel_ev(:,1),vel_ev(:,2))
        title('Velocity vs Time')
        xlabel('Time [s]')
        ylabel('Velocity [m/s]');

    figure(15)
        plot(Mach_ev(:,1),Mach_ev(:,2))
        title('Velocity vs Time')
        xlabel('Time [s]')
        ylabel('Velocity [Mach]');

    figure(16)
        plot(accel_ev(:,1),accel_ev(:,2))
        title('Acceleration vs Time')
        xlabel('Time [s]')
        ylabel('Acceleration [m/s2]');
    
    figure(17)
        plot(g_force(:,1),g_force(:,2))
        title('Acceleration vs Time')
        xlabel('Time [s]')
        ylabel('Acceleration [g]');
end

if ras
    ThrustCurve_big = [round(thrust_ev(:,1),2),round(thrust_ev(:,2),2)]; %Thrust curve for export
    ds_length = height(thrust_ev(:,1)); %length of the time matrix
    ds_freq = ds_length/31; %sampling frequency required
    ThrustCurveRA = round(ThrustCurve_big(1:ds_freq:end,:),3); %Thrust curve matrix with 30 elements needed for RASaero II
    ThrustCurveRA(32,2) = 0;
    ThrustCurveRA(32,1) = ds_length*timeStep;
    ThrustCurveRA(1,1) = ThrustCurveRA(2,1)/2;
    ThrustCurveRA(1,2) = ThrustCurveRA(2,2)/2;
end

%% Functions
% look up the property at a given temperature in the steam table
% where T is a temperature in degrees C and j is the horizontal index of
% the desired property
function [prop] = oxProp (T, j, props)
    i = 2*T+181;

    prop = props(ceil(i),j)*(1-(ceil(i)-i))+props(floor(i),j)*(ceil(i)-i);
end