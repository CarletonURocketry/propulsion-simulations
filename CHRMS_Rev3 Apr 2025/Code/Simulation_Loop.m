%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation_Loop
% 2025/03/28
% Adrian Comisso
%
% Desccription: 
% This function is the main loop of the simulation. It calls calculation
% functions tests, if the burn is complete and returns all the results  
%
% Inputs:
% state - a struct that stores the current state of each varible
%
% Outputs:
% output - a struct that stores all the past states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = Simulation_Loop(state,fuel)

%initialize output arrays and capture initial state
output.t = [state.t];
output.T_tank = [state.T_tank];
output.P_tank = [state.P_tank];
output.U_ox_total = [state.U_ox_total];
output.delta_P_tank = [state.delta_P_tank];
output.rho_ox_l = [state.rho_ox_l];
output.rho_ox_v = [state.rho_ox_v];
output.m_ox_total = [state.m_ox_total];
output.V_ox_l = [state.V_ox_l];
output.V_ox_v = [state.V_ox_v];
output.m_ox_l = [state.m_ox_l];
output.m_ox_v = [state.m_ox_v];
output.quality_ox = [state.quality_ox];
output.delta_m_ox_v = [state.delta_m_ox_v];
output.delta_T_tank = [state.delta_T_tank];
output.m_dot_inj = [state.m_dot_inj];
output.delta_P_inj = [state.delta_P_inj];
output.P_cmbr = [state.P_cmbr];
output.OF = [state.OF];
output.m_fuel = [state.m_fuel];
output.m_dot_fuel = [state.m_dot_fuel];
output.id_grn = [state.id_grn]; 
output.r_dot_grn = [state.r_dot_grn];
output.m_dot_noz = [state.m_dot_noz];
output.m_gas = [state.m_gas];
output.F_rocket = [state.F_rocket];
output.m_motor = [state.m_motor];
output.cg_motor = [state.cg_motor]; 

%main loop
while true  

    %Tank calcualtions
    if state.model_tank == 0
        [state] = HRAP_Tank(state,output);
    elseif state.model_tank == 1
        [state] = HRAP_Tank_Remastered(state,output);
    else %state.tank_method == 2
        [state] = Energy_Tank(state,output);
    end
    
    if state.cold_flow == false
        %regression calculations
        [state] = Regression(state);
        %combustion calculations
        [state] = Combustion(state,fuel);
        %chamber calculations
        [state] = Chamber(state);
        %nozzle calculations
        [state] = Nozzle(state);
    end

    %center of mass calcualtions
    [state] = Mass(state);

    %iterate time
    state.t = state.t + state.dt;
    
    %append data to output
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
    output.m_fuel = [output.m_fuel; state.m_fuel];
    output.m_dot_fuel = [output.m_dot_fuel; state.m_dot_fuel];
    output.id_grn = [output.id_grn; state.id_grn]; 
    output.r_dot_grn = [output.r_dot_grn; state.r_dot_grn];
    output.m_dot_noz = [output.m_dot_noz; state.m_dot_noz];
    output.m_gas = [output.m_gas; state.m_gas];
    output.F_rocket = [output.F_rocket; state.F_rocket];
    output.m_motor = [output.m_motor; state.m_motor];
    output.cg_motor = [output.cg_motor; state.cg_motor]; 

    %test for burnout conditions
    if state.m_ox_total <= 0.1
        output.end_cond = 'Oxidizer Depleted';
        break
    elseif state.t >= state.t_max
        output.end_cond = 'Max Simulation Time Reached';
        break
    elseif state.T_tank <= state.T_triple || state.T_tank >= state.T_crit
        output.end_cond = 'temperarture exited expected region';
        break
    elseif state.P_tank <= state.P_cmbr
        output.end_cond = 'motor backfired into tank (or cold flow ended)';
        break
    end

    if state.cold_flow == false
        if state.id_grn >= state.od_grn
            output.end_cond = 'Fuel Depleted';
            break
        elseif state.P_cmbr <= state.P_atm
            output.end_cond = 'Burn Complete';
            break
        end
    end

end
end