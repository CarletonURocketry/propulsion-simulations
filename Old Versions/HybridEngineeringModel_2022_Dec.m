% CU InSpace Hybrid Rocket Engineering Model
% Academic Year 2022-2023
% Last update Dec 1, 2022
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
flight = true;        % true - flight ; false - static fire
ras = true;           % are we putting this through RAS aero

%% Constants
g = 9.80665;               % [m/s2]    acceleration due to gravity
                           %           (ISO 80000)
R_univ = 8314.46261815324; % [J/kmolK] universal gas constant
R_air = R_univ/28.97;      % [J/kgK]   specific gas constant of air

%% Launch Site Conditions
T_amb = 36;     % [deg C] ambient temperature
T_tank = T_amb; % [deg C] tank temperature
P_amb = 85E3;   % [Pa]    ambient pressure
elev = 1400;    % [m]     elevation of launch site above sea level

P_sea = P_amb/(1-elev*2.25577E-5)^5.25588; % [Pa] pressure at sea level
                                           %      formula from engineer's
                                           %      toolbox

%% Vehicle Mass
m_ox = 13.5;                         % [kg] mass of oxidizer loaded
m_tank = 4.45;%4.0;                      % [kg] mass of the tank
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
H_tank = 1.4;%1.25;   % [m] tank height
R_tank = 0.0762; % [m] tank radius

V_tank = pi()*H_tank*R_tank^2; % [m3] tank volume

% injector
N_jet = 44;%50;%12;      %     number of injector holes
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
r_port = 3.5*0.0254/2;%2.375*25.4/2000; % [m]  initial fuel grain inner radius {optimal around 0.0336}
A_port = pi*r_port^2;     % [m2] initial fuel grain port area
L_fuel = 36*0.0254;       % [m]  fuel grain lenght
d_casing = 0.098;         % [m]  outer diameter of engine
t_casing = 0.0025;        % [m]  thickness of engine casing
t_liner = 0.08*0.0254;    % [m]  thickness of fuel liner

r_fuel_end = 4.25*0.0254/2 - t_liner;%d_casing/2-t_casing-t_liner; [3.75]
% ^ [m] outer radius of the fuel grain

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

% launch rail
l_rail = 33*0.3048; % [m] length of the launch rail

%% Propellant Properties
% fuel properties
rho_f = 920;         % [kg/m3]  fuel density

% oxidizer properties
oxProps = readmatrix('N2O_fluidProperties_SI.csv'); % Couch et al
M_ox = 44.013; % [kg/kmol] molar mass of N2O {currently unused}
P_super = 0;   % [Pa]      supercharging pressure {look into tank
               %                                   pressurization systems}

% ballistic coefficients
fuel_alpha = 0.0000876; % [m3/kg] estimate based on Liu et al
fuel_n = 0.3953;      %         estimate based on Liu et al

%% Vehicle Initial Conditions
alt = elev;     % [m]     initial altitude is launch site elevation
vel = 0;        % [m/s]   initial speed is zero {we aren't spinlaunch}

r_fuel_burn = r_port;                     % [m] set current fuel radius to
                                          %     initial radius

% first iteration needs initial chamber pressure estimate
P_c = 2.5E6; % [Pa]
% on all subsequent iterations, initial estimate will be pressure at
% previous iteration

thrust_max = 0;

% convert everything to Kelvin
T_amb = T_amb + 273.15;
T_tank = T_tank + 273.15;

% warning flags
pres_warn = false;
temp_warn = false;
volu_warn = false;

%% Evolution Over Time
rDot_ev = [];       % [mm/s] fuel regression rate
rem_fuel_ev = [];      % [m]    remaining fuel
rem_fuel_mass_ev = []; % [kg]   remaining mass of fuel
m_ox_ev = [];       % [kg]   remaining mass of oxidizer
mDot_ox_ev = [];    % [kg/s] oxidizer mass flow rate
G_ox_ev = [];       % [kg/m2s] oxidizer mass flux through fuel port
mDot_f_ev = [];     % [kg/s] fuel mass flow rate
mDot_t_ev = [];     % [kg/s] total propellant mass flow rate
OF_ev = [];         %        oxidizer/fuel mass ratio

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
    rem_fuel_ev(y+1,1) = y*timeStep;
    rem_fuel_mass_ev(y+1,1) = y*timeStep;
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
    OF_ev(y+1,1) = y*timeStep;
    G_ox_ev(y+1,1) = y*timeStep;

    %% tank pressure
    T_tank_ev(y+1,2) = T_tank;
    P_vap = interpolate(T_tank_ev(y+1,2),1,2,oxProps);
    P_tank_ev(y+1,2) = P_vap+P_super;
    rho_ox = interpolate(T_tank_ev(y+1,2),1,3,oxProps);

    %% mass flows, combustion chamber conditions, and exhaust properties
    % all of this is interdependent and so must be solved iteratively
    for z = 1:15
        % loop 15 times
        % this should get reasonably close to converging without tanking
        % compile time too hard

        %% oxidizer mass flux
        v_jet = (1/rho_ox)*Ccv_inj*Cd_inj*sqrt(2*rho_ox*(P_tank_ev(y+1,2) - P_vap*(1-C_bur)));
        mDot_ox = rho_ox*A_inj*v_jet + rho_ox*A_inj*v_jet*mDot_ox_noise*randn;

        G_ox = mDot_ox/(pi*r_fuel_burn^2);
    
        %% fuel regression and mass flow
        rDot = 1000*fuel_alpha*(G_ox)^fuel_n;

        dr_fuel_burn = rDot*timeStep/1000;

        mDot_f = rho_f*(rDot/1000)*2*pi*(r_fuel_burn + dr_fuel_burn)*L_fuel;

        mDot_t = mDot_ox+mDot_f;
        OF = mDot_t/mDot_f;

        %% chamber temperature
        % molar mass, ratio of specific heats, and chamber temperature are
        % based on curves fit to data from RPA
        
        % molar mass
        a_k = -7.27889933*(P_c/1E6) + 30.5345935;
        c_k = -0.115438083*(P_c/1E6) + 1.28230255;
        b_k = 1.16276656;

        k_c = a_k*exp(-c_k*OF) + b_k;
        % ^ ratio of specific heats of exhaust gasses

        % heat ratio
        a_M = -0.870565643*(P_c/1E6) - 27.5173943;
        c_M = 0.354710276;
        b_M = 0.133090595*(P_c/1E6) + 27.9812312;

        M_c = a_M*exp(-c_M*OF) + b_M;
        % ^ [kg/kmol] molar mass of exhaust gasses

        % temperature
        a_T = -123437.783*exp(-0.537956683*(P_c/1E6)) - 64150.1198;
        c_T = -0.00284431818*(P_c/1E6)^3 + 0.0298272727*(P_c/1E6)^2 - 0.146185925*(P_c/1E6) + 1.26928851;
        b_T = 3.06138889*(P_c/1E6)^3 - 31.0249405*(P_c/1E6)^2 + 141.377475*(P_c/1E6) + 3056.33191;

        T_c = a_T*exp(-c_T*OF) + b_T;
        % [K] temperature in combustion chamber

        R_c = R_univ/(M_c); % exhaust specific gas constant

        fcn = @(M) (1./M).*( (2/(k_c+1)) .* (1+(k_c-1)/2*M.^2) ) .^ ((k_c+1)/(2*(k_c-1))) - A_e/A_t;
        % ^ area ratio mach relationship (isentropic flow)

        Ma_e = fzero(fcn,3); % exhaust mach number

        P_ratio = (1+((k_c-1)*Ma_e^2)/2)^(-k_c/(k_c-1));
        % ^ ratio of exhaust pressure to chamber pressure

        C_F = sqrt(k_c*(2/(k_c+1))^((k_c+1)/(k_c-1)))*sqrt((2*k_c/(k_c-1))*(1-P_ratio^((k_c-1)/k_c))) + A_e*P_ratio/A_t;
        % ^ coefficient of thrust (Zolla et al)

        T_e = T_c/((1/P_ratio)^((k_c-1)/(k_c)));
        % ^ [K] exhaust temperature

        v_e = Ma_e*sqrt(k_c*R_c*T_e); % [m/s] exhaust velocity
        cstar = v_e/C_F;              % [m/s] characteristic velocity

        P_c = mDot_t*cstar/A_t;
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

    m_ox_ev(y+1,2) = m_ox;

    dm_ox_g = m_ox_g - m_ox_g_prev;%mDot_ox_ev(y+1,2)*v_SP_oxl/(v_SP_oxg-v_SP_oxl);
    dm_ox_l = - dm_ox_g - mDot_ox_ev(y+1,2);

    a = (m_tank*c_p_T+m_ox_l*c_v_oxl+m_ox_g*c_v_oxg); 
    % ^ very far off ideal gas model. not unexpected, ideal gas is very bad
    %   for NOS

    b = P_tank_ev(y+1,2)*v_SP_oxl;
    e = - h_vap + P_tank_ev(y+1,2)*v_SP_oxg;

    dT = (b*dm_ox_l+e*dm_ox_g)/a;
    T_tank = T_tank+dT*timeStep;

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
    end

    m_vehicle = m_dry + m_ox_ev(y+1,2) + rem_fuel_mass_ev(y+1,2); % [kg] vehicle mass

    twr_ev(y+1,2) = thrust_ev(y+1,2)/(g*m_vehicle);

    D_ev(y+1,2) = (C_drag*rho_air*vel^2*A_vehicle)/2; % drag

    %% flight kinematics
    if flight
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
        sav = y;
        break;
    end
    if P_tank_ev(y+1,2) <= P_c_ev(y+1,2)
        disp('Burnout: tank pressure less than chamber pressure')
        disp(['Remaining oxidizer: ',num2str(m_ox), 'kg'])
        disp(['Remaining fuel:     ',num2str(rem_fuel_mass_ev(y+1,2)), 'kg'])
        sav = y;
        break;
    end
    if m_ox <= 0
        disp('Burnout: out of oxidizer')
        disp(['Remaining fuel: ',num2str(rem_fuel_mass_ev(y+1,2)), 'kg'])
        sav = y;
        break;
    end

    if P_c_ev(y+1,2)>3E6
        pres_warn = true;
    end
    if T_c_ev(y+1,2)>3250
        temp_warn = true;
    end
end

disp(' ')
disp(['Burn duration: ',num2str(sav*timeStep),'s'])
disp(['Impulse: ',num2str(trapz(thrust_ev(:,2))*timeStep),'Ns'])
disp(['Maximum thrust: ',num2str(thrust_max),'N'])
disp(['Average thrust: ',num2str(trapz(thrust_ev(:,2))/sav),'N'])
disp(['Average Isp: ',num2str(trapz(Isp_ev(:,2))/sav),'s'])

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

P_delta = P_tank_ev(:,2) - P_c_ev(:,2);

%% Plots
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
    ylabel('Mass Flux [kg/m2s]')

figure(5)
    plot(OF_ev(:,1),OF_ev(:,2))
    title('Mixture Ratio vs Time')
    xlabel('Time [s]')
    ylabel('Mixture Ratio O/F')

% pressure and temperature plots
figure(6)
    plot(P_tank_ev(:,1),P_tank_ev(:,2),'b',P_c_ev(:,1),P_c_ev(:,2),'r')
    title('Pressure vs Time (tank and chamber)')
    xlabel('Time [s]')
    ylabel('Pressure [Pa]')
    legend('tank pressure','combustion chamber pressure');

figure(7)
    plot(P_e_ev(:,1),P_delta,'m')
    title('Pressure Drop Across Injector vs Time')
    xlabel('Time [s]')
    ylabel('Pressure Drop [Pa]')

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

if flight
    figure(14)
        plot(D_ev(:,1),D_ev(:,2))
        title('Drag vs Time')
        xlabel('Time [s]')
        ylabel('Drag [N]');
end

% kinematic plots
if flight
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
    ThrustCurve_big = [round(thrust_ev(:,1),2),round(thrust_ev(:,2),2)]; %Thrust curve for export
    ds_length = height(thrust_ev(:,1)); %length of the time matrix
    ds_freq = ds_length/31; %sampling frequency required
    warning('off','MATLAB:colon:nonIntegerIndex') % suppress colon opperator integer warning
    ThrustCurveRA = round(ThrustCurve_big(1:ds_freq:end,:),3); %Thrust curve matrix with 30 elements needed for RASaero II
    ThrustCurveRA(32,2) = 0;
    ThrustCurveRA(32,1) = ds_length*timeStep;
    ThrustCurveRA(1,1) = ThrustCurveRA(2,1)/2;
    ThrustCurveRA(1,2) = ThrustCurveRA(2,2)/2;

    disp(' ')
    disp('RAS input:')

    for i = 1:32
        disp([num2str(ThrustCurveRA(i,1)),' ',num2str(ThrustCurveRA(i,2))])
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
    B = interp1(table(:,i),table(:,j),A,'pchip');
end