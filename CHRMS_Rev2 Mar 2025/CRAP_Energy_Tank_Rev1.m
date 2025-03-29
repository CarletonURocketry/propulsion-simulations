%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy Tank Rev1
% 2024/09/02
% Adrian Comisso
%
% Desccription: 
% This function is calucaltes the new state of the run tank based on the
% previous iteration. The function tests weather the tank is still
% saturated or has reached a point where only vapor is left then performs
% the corresponding calculations to get tank pres and mass flow.  
%
% Inputs:
% param - stores the static parameters of the scenario (eg. tank volume)
% state - stores the current state of each varible
%
% Outputs:
% new_state - returns the updated state of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state] = Energy_Tank(param,state,output)
%***needs compersilbility too*********should use similar method to hrap
%check if any liquid remains
%********************************fluid*********************
if state.m_ox_l > 0
    %calcuatle pressure drop across injector (Pa)
    state.delta_P_inj = state.P_tank - state.P_cmbr;
    %calucalte mass flow of injector (kg/s)
    [m_dot,state] = Injector_Mass_Flow_Rev_1(param,state,param.model_inj_l,param.A_inj_total,param.Cd_inj,state.delta_P_inj);
    state.m_dot_inj = m_dot;
    %calcualte mass flow of vent (kg/s)
    if param.state_vent == 0
        state.delta_P_vent = 0;
        state.m_dot_vent = 0;
    elseif param.state_vent == 1
        state.delta_P_vent = state.P_tank - param.P_atm;
        [m_dot,state] = Injector_Mass_Flow_Rev_1(param,state,param.model_vent,param.A_vent,param.Cd_vent,state.delta_P_vent);
        state.m_dot_vent = m_dot;
    elseif param.state_vent == 2
        state.delta_P_vent = state.P_tank - param.P_cmbr;
        [m_dot,state] = Injector_Mass_Flow_Rev_1(param,state,param.model_vent,param.A_vent,param.Cd_vent,state.delta_P_vent);
        state.m_dot_vent = m_dot;
    end

    %%%%% New Mass/Energy %%%%%
    %calcualte new total mass (kg)
    state.m_ox_total = state.m_ox_total - state.m_dot_inj*param.dt - state.m_dot_vent*param.dt;
    %calculate new total internal energy (J)
    state.U_ox_total = state.U_ox_total - (state.m_dot_inj*param.dt*state.h_ox_l) - (state.m_dot_vent*param.dt*state.h_ox_v);
    
    %%%%% Initial Loop States %%%%%
    %loop counter
    i = 0;
    %initial delta T loop (gets imidiatley overwitten)
    delta_T_loop = 1;
    
    %%%%% Solve Temperature %%%%%
    %test if change meets required tolerance
    while abs(delta_T_loop) >= 0.000001
        %notify user if tank temperature is no longer saturated liquid
        %************************not specific to liquid type *******************************
        if state.T_tank >=  param.T_crit || state.T_tank <= param.T_triple
            fprintf("Temparute outside of saturation dome")
            break
        end
        %liquid oxidizer specific internal energy (J/kg)  
        u_ox_l = CoolProp.PropsSI('U','T',state.T_tank,'Q',0,param.ox_name);
        %liquid oxidizer specific internal energy (J/kg)
        u_ox_v = CoolProp.PropsSI('U','T',state.T_tank,'Q',1,param.ox_name);
        %liquid oxidizer density (kg/m^3)
        state.rho_ox_l = CoolProp.PropsSI('D','T',state.T_tank,'Q',0,param.ox_name);
        %vapor oxidizer density (kg/m^3)
        state.rho_ox_v = CoolProp.PropsSI('D','T',state.T_tank,'Q',1,param.ox_name);
        %quality
        state.quality_ox = ((state.U_ox_total/state.m_ox_total) - u_ox_l)/(u_ox_v-u_ox_l);
        if state.quality_ox > 1
            state.quality_ox = 1;
        end
        %difference from real tank volume (m^3)
        V_diff = param.V_tank - state.m_ox_total*(((1-state.quality_ox)/state.rho_ox_l) + (state.quality_ox/state.rho_ox_v));
        %delta T for this loop (K)
        delta_T_loop = state.T_tank*(1-V_diff) - state.T_tank;
        %T tank (K)
        state.T_tank = state.T_tank + delta_T_loop;
        
        %iterate counter
        i = i+1;
        %check if loop has blown up
        if i > 10000
            fprintf('loop blew up')
            break
        end
    end
    
    %%%%% New State %%%%%
    %liquid oxidizer specific enthalpy (J/kg)
    state.h_ox_l = CoolProp.PropsSI('H','T',state.T_tank,'Q',0,param.ox_name);
    %vapor oxidizer specific enthalpy (J/kg)
    state.h_ox_v = CoolProp.PropsSI('H','T',state.T_tank,'Q',1,param.ox_name);
    %change in tank pressure (Pa)
    state.delta_P_tank = CoolProp.PropsSI('P','T',state.T_tank,'Q',state.quality_ox,param.ox_name) - state.P_tank;
    %tank pressure (Pa)
    state.P_tank = CoolProp.PropsSI('P','T',state.T_tank,'Q',state.quality_ox,param.ox_name);
    %mass of liquid oxidizer (kg)
    state.m_ox_l = (1-state.quality_ox)*state.m_ox_total;
    %mass of vapor oxidizer (kg)
    state.m_ox_v = state.quality_ox*state.m_ox_total;
else %state.m_ox_l <= 0
    %correct negative mass of liquid****some mass will disapear*****
    state.m_ox_l = 0;
    
    %calcuatle pressure drop across injector (Pa)
    state.delta_P_inj = state.P_tank - state.P_cmbr;
    %calucalte m_dot_oxidizer (kg/s)
    [m_dot,state] = Injector_Mass_Flow_Rev_1(param,state,param.model_inj_v,param.A_inj_total,param.Cd_inj,state.delta_P_inj);
    state.m_dot_inj = m_dot;
    %calcualte mass flow of vent (kg/s)
    if param.state_vent == 0
        state.delta_P_vent = 0;
        state.m_dot_vent = 0;
    elseif param.state_vent == 1
        state.delta_P_vent = state.P_tank - param.P_atm;
        [m_dot,state] = Injector_Mass_Flow_Rev_1(param,state,param.model_vent,param.A_vent,param.Cd_vent,state.delta_P_vent);
        state.m_dot_vent = m_dot;
    elseif param.state_vent == 2
        state.delta_P_vent = state.P_tank - param.P_cmbr;
        [m_dot,state] = Injector_Mass_Flow_Rev_1(param,state,param.model_vent,param.A_vent,param.Cd_vent,state.delta_P_vent);
        state.m_dot_vent = m_dot;
    end

    %%%%% Initiate Varibles %%%%%
    %mass at beginng of time step (kg)
    m_ox_1 = state.m_ox_total;
    %calcualte new total mass (kg)
    state.m_ox_total = state.m_ox_total - state.m_dot_inj*param.dt - state.m_dot_vent*param.dt;
    %calculate new total internal energy (J)
    state.U_ox_total = state.U_ox_total - (state.m_dot_inj*param.dt*state.h_ox_v) - (state.m_dot_vent*param.dt*state.h_ox_v);
    %mass at end of time step (kg)
    m_ox_2 = state.m_ox_total;
    %mass of vapor oxidizer (kg)
    state.m_ox_v = state.m_ox_total;
    
    %find Cp (J/kg*K)
    Cp = CoolProp.PropsSI('CP0MASS','T',293.15,'P',101325,param.ox_name);
    %find Cv (J/kg*K)
    Cv = CoolProp.PropsSI('CVMASS','T',293.15,'P',101325,param.ox_name);
    %find specific heat ratio gamma
    gamma = Cp/Cv;

    blowdown_method = "real gas"; %********make this a parameter ********

    if blowdown_method == "ideal gas"
        %clauclate new temperature (K)
        T_2 = state.T_tank * (m_ox_2/m_ox_1)^(gamma-1);
        %clauclate new pressure (Pa)
        P_2 = state.P_tank * (m_ox_2/m_ox_1)^gamma;
        %change in tank pressure (Pa)
        state.delta_P_tank = P_2 - state.P_tank;
        %%%%% New State %%%%%
        %new tank temp (K)
        state.T_tank = T_2;
        %tank pressure (Pa)
        state.P_tank = P_2;
        %vapor oxidizer density (kg/m^3)
        state.rho_ox_v = state.m_ox_total / param.V_tank;%CoolProp.PropsSI('D','T',state.T_tank,'P',state.P_tank,param.ox_name);
        %vapor oxidizer specific enthalpy (J/kg)
        state.h_ox_v = CoolProp.PropsSI('H','T',state.T_tank,'P',state.P_tank,param.ox_name);
    elseif blowdown_method == "real gas"
        %compressibility factor at beginng of time step
        Z_1 = CoolProp.PropsSI('Z','T',state.T_tank,'Q',1,param.ox_name);
        %initial compressibility factor guess
        Z_2 = Z_1;
    
        %%%%% Initial Loop States %%%%%
        %loop counter
        i = 0;
        %initial delta T loop (gets imidiatley overwitten)
        delta_Z_loop = 1;
    
        %%%%% Solve Compressiblilty %%%%%
        while delta_Z_loop >= 0.000001 %can get smaller test sensitivity ******
            %calcualte T2 
            T_2 = state.T_tank*(((Z_2*m_ox_2)/(Z_1*m_ox_1))^(gamma-1));
            %claucalte P2
            P_2 = state.P_tank*((T_2/state.T_tank)^(gamma/(gamma-1)));
            %calculate change in Z
            delta_Z_loop = abs(Z_2 - CoolProp.PropsSI('Z','T',T_2,'Q',1,param.ox_name));
            %assing values for next iteration
            Z_2 = (Z_2 + CoolProp.PropsSI('Z','T',T_2,'Q',1,param.ox_name))/2;
    
            %iterate loop counter
            i = i+1;

            %notify user if loop did not convererge
            if i>1000000
                disp("vapor loop blew up")
                break
            end

            %notify user number of loops
            fprintf("number of loops:" + i)
            
        end

        %%%%% New State %%%%%
        %new tank pressure (Pa)
        state.P_tank = P_2;%P_2 = CoolProp.PropsSI('P','T',state.T_tank,'Q',1,param.ox_name);
        %new tank temp (K)
        state.T_tank = T_2;
        %change in tank pressure (Pa)
        state.delta_P_tank = P_2 - state.P_tank;
        %tank pressure (Pa)
        state.P_tank = P_2;
        %vapor oxidizer density (kg/m^3)
        state.rho_ox_v = CoolProp.PropsSI('D','T',state.T_tank,'Q',1,param.ox_name);
        %vapor oxidizer specific enthalpy (J/kg)
        state.h_ox_v = CoolProp.PropsSI('H','T',state.T_tank,'Q',1,param.ox_name);
    end

    %vapor phase notification
    disp("vapor phase")
end
end
