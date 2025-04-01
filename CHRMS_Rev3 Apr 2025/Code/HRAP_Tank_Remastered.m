%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HRAP Tank Remastered
% 2025/03/28
% Adrian Comisso
%
% Desccription: 
% This function calucaltes the new state of the run tank based on the
% previous iteration. This version uses a revised copy of the HRAP model.
% For more information see the attached Theory of Operation docuemnt.
%
% Inputs:
% state - a struct that stores the current state of each varible
% output - a struct that stores all the past states
%
% Outputs:
% state - returns the updated state of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state] = HRAP_Tank_Remastered(state,output)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store State at beggining of iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get oxidizer porperties
oxidizer = Nitrous_Properties(state.T_tank);
%store change in pressure from last step (Pa)
state.delta_P_tank = oxidizer.P_v - state.P_tank;
%assign new tank pressure (Pa)
state.P_tank = oxidizer.P_v;
%assign new density of liquid (kg/m^3)
state.rho_ox_l = oxidizer.rho_l;
%assign new density of vapor (kg/m^3)
state.rho_ox_v = oxidizer.rho_v;

%make sure initial time step enters correct tank phase
if state.t == 0
    state.m_l_unvaporized = state.m_ox_l + 0.0001;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phases of Empyting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check the state of the tank
if state.m_ox_l < state.m_l_unvaporized && state.m_ox_l > 0

    %%%%% Find Mass Flow Rates %%%%%
    %calcuatle pressure drop across injector (Pa)
    state.delta_P_inj = state.P_tank - state.P_cmbr;
    %calucalte m_dot_inj (kg/s)
    [m_dot,state] = Mass_Flow(state,state.model_inj_l,state.A_inj_total,state.Cd_inj,state.delta_P_inj);
    state.m_dot_inj = m_dot;
    %calcualte m_dot_vent (kg/s)
    if state.state_vent == 0
        state.delta_P_vent = 0;
        state.m_dot_vent = 0;
    elseif state.state_vent == 1
        state.delta_P_vent = state.P_tank - state.P_atm;
        [m_dot,state] = Mass_Flow(state,state.model_vent,state.A_vent,state.Cd_vent,state.delta_P_vent);
        state.m_dot_vent = m_dot;
    elseif state.state_vent == 2
        state.delta_P_vent = state.P_tank - state.P_cmbr;
        [m_dot,state] = Mass_Flow(state,state.model_vent,state.A_vent,state.Cd_vent,state.delta_P_vent);
        state.m_dot_vent = m_dot;
    end
    
    %%%%% Update & Rebalance Tank Masses %%%%% 
    %calcualte the reaming mass of liquid if none of it vaporized (kg)
    state.m_l_unvaporized = state.m_ox_l - (state.m_dot_inj*state.dt);
    %calcualte the remaining mass of oxidizer in the tank (kg)
    state.m_ox_total = state.m_ox_total - (state.m_dot_inj*state.dt) - (state.m_dot_vent*state.dt);
    %calcualte the remaing mass of liquid if some does vaporize (kg)
    state.m_ox_l = (state.V_tank - (state.m_ox_total/oxidizer.rho_v))/((1/oxidizer.rho_l)-(1/oxidizer.rho_v));
    %calcualte mass vaporized (kg)
    state.delta_m_ox_v = state.m_l_unvaporized - state.m_ox_l;
    
    %%%%% Rebalance Tank Temperature & Pressure %%%%%
    %find heat removed by the vaporising oxidizer(J)
    delta_Q = state.delta_m_ox_v * oxidizer.delta_H;
    %find change in tank temperature (K)
    state.delta_T_tank = -delta_Q/(state.m_ox_l*oxidizer.Cp);
    %find the new oxidizer liquid temperature (K)
    state.T_tank = state.T_tank + state.delta_T_tank;

elseif state.m_ox_l >= state.m_l_unvaporized && state.m_ox_l > 0

    %%%%% Find Mass Flow Rates %%%%% 
    %calcuatle pressure drop across injector (Pa)
    state.delta_P_inj = state.P_tank - state.P_cmbr;
    %calucalte m_dot_inj (kg/s)
    [m_dot,state] = Mass_Flow(state,state.model_inj_l,state.A_inj_total,state.Cd_inj,state.delta_P_inj);
    state.m_dot_inj = m_dot;
    %calcualte m_dot_vent (kg/s)
    if state.state_vent == 0
        state.delta_P_vent = 0;
        state.m_dot_vent = 0;
    elseif state.state_vent == 1
        state.delta_P_vent = state.P_tank - state.P_atm;
        [m_dot,state] = Mass_Flow(state,state.model_vent,state.A_vent,state.Cd_vent,state.delta_P_vent);
        state.m_dot_vent = m_dot;
    elseif state.state_vent == 2
        state.delta_P_vent = state.P_tank - state.P_cmbr;
        [m_dot,state] = Mass_Flow(state,state.model_vent,state.A_vent,state.Cd_vent,state.delta_P_vent);
        state.m_dot_vent = m_dot;
    end
    
    %%%%% Update & Rebalance Tank Masses %%%%%
    %calcualte the remaining mass of oxidizer in the tank (kg)
    state.m_ox_total = state.m_ox_total - state.m_dot_inj*state.dt - state.m_dot_vent*state.dt;
    %set unvaporized mass (kg)
    state.m_l_unvaporized = 0;
    %calcualte the remaing mass of liquid (kg) 
    state.m_ox_l = (state.V_tank - (state.m_ox_total/oxidizer.rho_v))/ ...
                    ((1/oxidizer.rho_l)-(1/oxidizer.rho_v));

    %%%%% Rebalance Tank Temperature & Pressure %%%%%
    %clauclate average pressure drop (Pa)
    delta_P_tank_avg = mean(output.delta_P_tank(1:sum(output.delta_P_tank<0)));
    %temporary P_tank (Pa)
    P_tank_temp = state.P_tank + delta_P_tank_avg;
    %temporary vapor pressure equation of state
    vp = @(T) 7251000*exp((1/(T/309.57))*...
        (-6.71893*(1-T/309.57) + 1.35966*(1-(T/309.57))^(3/2) + -1.3779*...
        (1-(T/309.57))^(5/2) + -4.051*(1-(T/309.57))^5)) - P_tank_temp;
    %solve for the tank temperature given the temporary P_tank (K)
    state.T_tank = fzero(vp,state.T_tank);

else %vapor phase
    
    %%%%% Find Mass Flow Rates %%%%% 
    %calcuatle pressure drop across injector (Pa)
    state.delta_P_inj = state.P_tank - state.P_cmbr;
    %correct remaining liquid (kg)
    state.m_ox_l = 0;
    %calucalte m_dot_inj (kg/s)
    [m_dot,state] = Mass_Flow(state,state.model_inj_v,state.A_inj_total,state.Cd_inj,state.delta_P_inj);
    state.m_dot_inj = m_dot;
    %calcualte mass flow of vent (kg/s)
    if state.state_vent == 0
        state.delta_P_vent = 0;
        state.m_dot_vent = 0;
    elseif state.state_vent == 1
        state.delta_P_vent = state.P_tank - state.P_atm;
        [m_dot,state] = Mass_Flow(state,state.model_vent,state.A_vent,state.Cd_vent,state.delta_P_vent);
        state.m_dot_vent = m_dot;
    elseif state.state_vent == 2
        state.delta_P_vent = state.P_tank - state.P_cmbr;
        [m_dot,state] = Mass_Flow(state,state,state.model_vent,state.A_vent,state.Cd_vent,state.delta_P_vent);
        state.m_dot_vent = m_dot;
    end

    %%%%% Numerical Loop Varibles %%%%%
    %mass at beginng of time step (kg)
    m_ox_1 = state.m_ox_total;
    %calcualte m_ox_total (kg)
    state.m_ox_total = state.m_ox_total - state.m_dot_inj*state.dt - state.m_dot_vent*state.dt;
    %mass at end of time step (kg)
    m_ox_2 = state.m_ox_total;
    %init tank temp (K)
    T_1 = state.T_tank;
    %init tank pressure (Pa)
    P_1 = state.P_tank;
    %init Z guess
    Z_guess = oxidizer.Z;
    Z_old = Z_guess;
    %loop variables
    i = 0;
    delta_Z = 1;

    %%%%% Numerical Vapor Phase Loop %%%%%
    while delta_Z >= 0.000001 %can get smaller test sensitivity ******
        %calcualte T2 (K)
        T_2 = T_1*(((Z_guess*m_ox_2)/(Z_old*m_ox_1))^(oxidizer.gamma - 1));
        %calcualte P2 (Pa)
        P_2 = P_1*((T_2/T_1)^(oxidizer.gamma/(oxidizer.gamma - 1))); 
        %get nitrous properties at T2
        oxidizer = Nitrous_Properties(T_2);
        %calculate change in Z
        delta_Z = abs(Z_guess - oxidizer.Z);
        %assing values for next iteration
        Z_guess = (Z_guess + oxidizer.Z)/2;

        %iterate counter
        i = i+1;
        %test if loop has blown up
        if i>1000000
            fprintf("\nERROR: numerical vapor loop did not converge on time step: t = %g \nSee HRAP_Tank_Remastered.m line 193", state.t)
            break
        end
    end
    
    %%%%% Update Tank Temperature & Pressure %%%%%
    %store change in pressure (Pa)
    state.delta_P_tank = P_2 - state.P_tank;
    %store loop results
    state.T_tank = T_2;
    state.P_tank = P_2;
    

    %perhaps assign null or zero to values that arent used here?****
end

end