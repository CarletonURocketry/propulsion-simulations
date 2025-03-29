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

%initialize output arrays and capture initial state
% (can make this one large array with titles?**)
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
output.m_dot_fuel = [state.m_dot_fuel];
output.id_grn = [state.id_grn]; 
output.r_dot_grn = [state.r_dot_grn];
output.F_rocket = [state.F_rocket];
output.m_motor = [state.m_motor];
output.cg_motor = [state.cg_motor]; 

%main loop
while true  

    %Tank calcualtions
    if param.model_tank == 0
        [state] = CRAP_HRAP_Tank_Rev1(param,state,output);
    elseif param.model_tank == 1
        [state] = CRAP_HRAP_Tank_Remastered_Rev1(param,state,output);
    else %param.tank_method == 2
        [state] = CRAP_Energy_Tank_Rev1(param,state,output);
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

    %center of mass calcualtions
    [state] = CRAP_COM_Rev1(param,state);

    %iterate time
    state.t = state.t + param.dt;
    
    % *************debug ***************
    i = state.t/param.dt;
    if state.t >= (8.276 - 0.0000001)%stupid floating point bulshit
        fprintf("AC time---------")
        disp(state.t)
        fprintf("T_tank")
        disp(state.T_tank)
        fprintf("P_tank")
        disp(state.P_tank)
        fprintf("m_ox_total")
        disp(state.m_ox_total)
        fprintf("m_dot_inj")
        disp(state.m_dot_inj)
        fprintf("P_cmbr----------")
        disp(state.P_cmbr)
        disp("")
        %output.end_cond = 'artificial end';
        %break
    end
    % ***************debug ***********************

    %append data to output ***** should move up above t???**************
    %would go faster if full list was initilized first********
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
    output.m_motor = [output.m_motor; state.m_motor];
    output.cg_motor = [output.cg_motor; state.cg_motor]; 


    %test for burnout conditions
    if state.m_ox_total <= 0.1
        output.end_cond = 'Oxidizer Depleted';
        break
    elseif state.t >= param.t_max
        output.end_cond = 'Max Simulation Time Reached';
        break
    elseif state.T_tank <= param.T_triple || state.T_tank >= param.T_crit
        output.end_cond = 'temperarture exited expected region';
        break
    elseif state.P_tank <= state.P_cmbr
        output.end_cond = 'motor backfired into tank (or cold flow ended)';
        break
    end

    if param.cold_flow == false
        if state.id_grn >= param.od_grn
            output.end_cond = 'Fuel Depleted';
            break
        elseif state.P_cmbr <= param.P_atm
            output.end_cond = 'Burn Complete';
            break
        end
    end

end
end