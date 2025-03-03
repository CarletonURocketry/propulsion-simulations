% Downstream Pressure vs Mass Flow Rate Comparison for SPI and HEM methods
% Adrian Comiss

% Define constants and initial conditions
param.ox_name = 'N2O';   % Oxidizer (Nitrous Oxide in this case)
state.T_tank = 304;      % Tank temperature (293)
state.P_tank = CoolProp.PropsSI('P','T',state.T_tank,'Q',0,param.ox_name);
state.rho_ox_l = CoolProp.PropsSI('D','T',state.T_tank,'Q',0,param.ox_name);    % Density of liquid N2O (kg/m^3)
state.rho_ox_v = CoolProp.PropsSI('D','T',state.T_tank,'Q',1,param.ox_name);     % Density of vapor N2O (kg/m^3)
state.t = 0; %set time to zero to not trip up HEM mdoel

Cd = 0.7;                % Coefficient of discharge (assumed value)
A = 0.0001;              % Orifice area (m^2)
delta_P_range = linspace(1e5, 4.9e6, 100); % Range of pressure drops (Pa)

% Allocate arrays for storing results
m_dot_SPI = zeros(size(delta_P_range));
m_dot_HEM = zeros(size(delta_P_range));
m_dot_Dyer = zeros(size(delta_P_range));

% Loop over different pressure drops
for i = 1:length(delta_P_range)
    delta_P = delta_P_range(i);
    
    % SPI method
    m_dot_SPI(i) = Injector_Mass_Flow_Rev_1(param, state, 'SPI', A, Cd, delta_P);
    
    % HEM method
    m_dot_HEM(i) = Injector_Mass_Flow_Rev_1(param, state, 'HEM', A, Cd, delta_P);

    % Dyer method
    m_dot_Dyer(i) = 0.5*m_dot_HEM(i) + 0.5*m_dot_SPI(i);
    %m_dot_Dyer(i) = Injector_Mass_Flow_Rev_1(param, state, 'Dyer', A, Cd, delta_P);
end

% Plotting results
figure;
plot(delta_P_range / 1e6, m_dot_SPI, 'b', 'LineWidth', 2); hold on;
plot(delta_P_range / 1e6, m_dot_HEM, 'r--', 'LineWidth', 2);
plot(delta_P_range / 1e6, m_dot_Dyer, 'g.', 'LineWidth', 2);
xlabel('Downstream Pressure Drop (MPa)');
ylabel('Mass Flow Rate (kg/s)');
legend('SPI Method', 'HEM Method','Dyer Method');
title('Comparison of Mass Flow Rate vs Pressure Drop for SPI and HEM');
grid on;
