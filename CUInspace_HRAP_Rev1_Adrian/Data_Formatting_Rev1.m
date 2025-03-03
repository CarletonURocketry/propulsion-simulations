%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_Formatting 
% 2024/11/06
% Adrian Comisso
%
% Desccription: 
% prgram that displaz test data and to generally format data with
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%which data flow data would you like to import
%CF1_1 = Cold flow 1.1(first ever cold flow)
%CF1_1 = Cold flow 2.1(quarry cold flow  that didnt happen only vapor phase)
%CF1_1 = Cold flow 2.2
%CF1_1 = Cold flow 3.1(october 23 2024 colf flow )
%CF1_1 = Cold flow 3.2
data = 'CF3_1';

%file paths for cold flows
CF1_1_path = 'C:\Users\adric\OneDrive\Desktop\4-CUInSpace\Hybrid Sim\HRAP_Decompiled_Matlab_Versions\CUInspace_HRAP_Rev1\Cold Flow Data\CF1\CF1_1';
CF2_1_path = 'C:\Users\adric\OneDrive\Desktop\4-CUInSpace\Hybrid Sim\HRAP_Decompiled_Matlab_Versions\CUInspace_HRAP_Rev1\Cold Flow Data\CF2\CF2_1';
CF2_2_path = 'C:\Users\adric\OneDrive\Desktop\4-CUInSpace\Hybrid Sim\HRAP_Decompiled_Matlab_Versions\CUInspace_HRAP_Rev1\Cold Flow Data\CF2\CF2_2';
CF3_1_path = 'C:\Users\adric\OneDrive\Desktop\4-CUInSpace\Hybrid Sim\HRAP_Decompiled_Matlab_Versions\CUInspace_HRAP_Rev1\Cold Flow Data\CF3\CF3_1';
CF3_2_path = 'C:\Users\adric\OneDrive\Desktop\4-CUInSpace\Hybrid Sim\HRAP_Decompiled_Matlab_Versions\CUInspace_HRAP_Rev1\Cold Flow Data\CF3\CF3_2';

%load data
if data == 'CF1_1'
    file = 'CF1_1_time.mat';
    load(fullfile(CF1_1_path, file));
    file = 'CF1_1_m_ox.mat';
    load(fullfile(CF1_1_path, file));
    file = 'CF1_1_P_tank.mat';
    load(fullfile(CF1_1_path, file));
    file = 'CF1_1_T_tank.mat';
    load(fullfile(CF1_1_path, file));
    file = 'CF1_1_T_tank_calculated.mat';
    load(fullfile(CF1_1_path, file));
    file = 'CF1_1_mdot.mat';
    load(fullfile(CF1_1_path, file));
elseif data == 'CF2_1'
    % didnt bother to write this cuz only vapor phase
    disp('neglected to fill out import staments for this colf flow')
elseif data == 'CF2_2'
    % didnt bother to write this cuz only vapor phase
    disp('neglected to fill out import staments for this colf flow')
elseif data == 'CF3_1'
    file = 'CF3_1_time.mat';
    load(fullfile(CF3_1_path, file));
    file = 'CF3_1_m_ox.mat';
    load(fullfile(CF3_1_path, file));
    file = 'CF3_1_P_tank.mat';
    load(fullfile(CF3_1_path, file));
    file = 'CF3_1_P_inj.mat';
    load(fullfile(CF3_1_path, file));
    file = 'CF3_1_T_tank.mat';
    load(fullfile(CF3_1_path, file));
    file = 'CF3_1_mdot.mat';
    load(fullfile(CF3_1_path, file));
    file = 'CF3_1_P_source.mat';
    load(fullfile(CF3_1_path, file));
    file = 'CF3_1_mdot_ideal.mat';
    load(fullfile(CF3_1_path, file));
    file = 'CF3_1_Cd.mat';
    load(fullfile(CF3_1_path, file));
elseif data == 'CF3_2'
    file = 'CF3_2_time.mat';
    load(fullfile(CF3_2_path, file));
    file = 'CF3_2_m_ox.mat';
    load(fullfile(CF3_2_path, file));
    file = 'CF3_2_P_tank.mat';
    load(fullfile(CF3_2_path, file));
    file = 'CF3_2_P_inj.mat';
    load(fullfile(CF3_2_path, file));
    file = 'CF3_2_T_tank.mat';
    load(fullfile(CF3_2_path, file));
    file = 'CF3_2_mdot.mat';
    load(fullfile(CF3_2_path, file));
    file = 'CF3_2_P_source.mat';
    load(fullfile(CF3_2_path, file));
    file = 'CF3_2_mdot_ideal.mat';
    load(fullfile(CF3_1_path, file));
    file = 'CF3_2_Cd.mat';
    load(fullfile(CF3_1_path, file));
end
    
hold on
figure(1);
plot(CF3_1_time,CF3_1_P_tank)
hold on
figure(2);
plot(CF3_1_time,CF3_1_T_tank)
hold on
figure(3);
plot(CF3_1_time,CF3_1_mdot)
hold on
figure(4);
plot(CF3_1_time,CF3_1_m_ox)
hold on
