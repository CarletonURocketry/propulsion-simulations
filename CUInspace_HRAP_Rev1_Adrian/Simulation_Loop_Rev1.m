%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation_Loop_Rev1
% 2024/08/01
% Adrian Comisso
%
% Desccription: 
% This function is the main time loop for the simulation it calls the
% calculation functions tests if the burn is complete and returns all the  
%
% Inputs:
% param - stores the static parameters of the scenario (eg. tank volume)
% state - stores the current state of each varible
%
% Outputs:
% output - stores the past state of all varibles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = Simulation_Loop(param,state,fuel)

%initilaize time
state.t = 0;

%initialize output arrays
output.t = [];
output.T_tank = [];
output.P_tank = [];
output.U_ox_total = [];
output.delta_P_tank = [];
output.rho_ox_l = [];
output.rho_ox_v = [];
output.m_ox_total = [];
output.V_ox_l = [];
output.V_ox_v = [];
output.m_ox_l = [];
output.m_ox_v = [];
output.quality_ox = [];
output.delta_m_ox_v = [];
output.delta_T_tank = [];
output.m_dot_inj = [];
output.delta_P_inj = [];
output.P_cmbr = [];
output.OF = [];
output.m_dot_fuel = [];
output.id_grn = []; 
output.r_dot_grn = [];
output.F_rocket = [];

%main loop
while true  

    %*************debug ***************
    % i = state.t/param.dt;
    % if state.t >= 11.19%10215%15730 %8257 %5466
    %     fprintf("AC time---------")
    %     disp(i)
    %     fprintf("T_tank")
    %     disp(state.T_tank)
    %     fprintf("P_tank")
    %     disp(state.P_tank)
    %     fprintf("m_ox_total")
    %     disp(state.m_ox_total)
    %     fprintf("m_dot_inj")
    %     disp(state.m_dot_inj)
    %     fprintf("P_cmbr----------")
    %     disp(state.P_cmbr)
    %     disp("")
    %     %output.end_cond = 'artificial end';
    %     %break
    % end
    %***************debug ***********************
    
    %Tank calcualtions
    if param.model_tank == 0
        [state] = Tank_Rev2(param,state,output);
    else % param.tank_method == 1
        [state] = Energy_Tank_Rev1(param,state,output);
    end
    
    if param.cold_flow == false
        %regression calculations
        [state] = Regression_Rev1(param,state);
        
        %combustion calculations
        [state] = Combustion_Rev1(param,state,fuel);
        
        %chamber calculations
        [state] = Chamber_Rev1(param,state);
    
        %nozzle calculations
        [state] = Nozzle_Rev1(param,state);
    end
    
    %mass calucaltions ******************

    %iterate time
    state.t = state.t + param.dt;

    %test for burnout conditions
    if state.id_grn >= param.od_grn
        output.end_cond = 'Fuel Depleted';
        break
    elseif state.m_ox_total <= 0.1
        output.end_cond = 'Oxidizer Depleted';
        break
    elseif state.t >= param.t_max
        output.end_cond = 'Max Simulation Time Reached';
        break
    elseif state.P_cmbr < param.P_atm
        output.end_cond = 'Burn Complete';
        break
    elseif state.T_tank <= param.T_triple || state.T_tank >= param.T_crit
        output.end_cond = 'temperarture exited expected region';
        break
    end
    
    %append data to output ***** should move up above t???**************
    output.t = [output.t; state.t];
    output.T_tank = [output.T_tank; state.T_tank];
    output.P_tank = [output.P_tank; state.P_tank];
    output.U_ox_total = [output.U_ox_total; state.U_ox_total];
    output.delta_P_tank = [output.delta_P_tank; state.delta_P_tank];
    output.rho_ox_l = [output.rho_ox_l; state.rho_ox_l];
    output.rho_ox_v = [output.rho_ox_v; state.rho_ox_v];
    output.m_ox_total = [output.m_ox_total; state.m_ox_total];
    output.V_ox_l = [output.V_ox_l; state.V_ox_l];
    output.V_ox_v = [output.V_ox_v; state.V_ox_v];
    output.m_ox_l = [output.m_ox_l; state.m_ox_l];
    output.m_ox_v = [output.m_ox_v; state.m_ox_v];
    output.quality_ox = [output.quality_ox; state.quality_ox];
    output.delta_m_ox_v = [output.delta_m_ox_v; state.delta_m_ox_v];
    output.delta_T_tank = [output.delta_T_tank; state.delta_T_tank];
    output.m_dot_inj = [output.m_dot_inj; state.m_dot_inj];
    output.delta_P_inj = [output.delta_P_inj; state.delta_P_inj];
    output.P_cmbr = [output.P_cmbr; state.P_cmbr];
    output.OF = [output.OF; state.OF];
    output.m_dot_fuel = [output.m_dot_fuel; state.m_dot_fuel];
    output.id_grn = [output.id_grn; state.id_grn]; 
    output.r_dot_grn = [output.r_dot_grn; state.r_dot_grn];
    output.F_rocket = [output.F_rocket; state.F_rocket];

end
end