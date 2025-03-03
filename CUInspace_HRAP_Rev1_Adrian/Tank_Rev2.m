%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tank_Rev2
% 2024/08/01
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

function [state] = Tank(param,state,output)

%get oxidizer porperties at beginning of time step
oxidizer = Nitrous_Properties_Rev1(state.T_tank);
%store change in pressure from last step
state.delta_P_tank = oxidizer.P_v - state.P_tank;
%assign new tank pressure (Pa)
state.P_tank = oxidizer.P_v;
%assign new density of liquid (kg/m^3)
state.rho_ox_l = oxidizer.rho_l;
%assign new density of vapor (kg/m^3)
state.rho_ox_v = oxidizer.rho_v;

%make sure initial time step enters correct category
if state.t == 0
    state.m_l_unvaporized = state.m_ox_l + 0.0001;
end

%%%%% Test for Only Vapor %%%%%
%***should this trigger with liquid still in system or not?*****
if state.m_ox_l < state.m_l_unvaporized && state.m_ox_l > 0
    %calcuatle pressure drop across injector (Pa)
    state.delta_P_inj = state.P_tank - state.P_cmbr;
    %calucalte m_dot_oxidizer
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
        
    %calcualte the reaming mass of liquid if none of it vaporized
    state.m_l_unvaporized = state.m_ox_l - (state.m_dot_inj*param.dt);
    %calcualte the remaining mass of oxidizer in the tank
    state.m_ox_total = state.m_ox_total - (state.m_dot_inj*param.dt) - (state.m_dot_vent*param.dt);
    %calcualte the remaing mass of liquid if some does vaporize
    state.m_ox_l = (param.V_tank - (state.m_ox_total/oxidizer.rho_v))/((1/oxidizer.rho_l)-(1/oxidizer.rho_v));
    %calcualte mass vaporized
    state.delta_m_ox_v = state.m_l_unvaporized - state.m_ox_l;
    
    %find heat removed by the vaporising oxidizer(J)***might be sensative to initial value
    %*** potnetially calcualte volme of injector manifold and find m given dens
    %*** potentially bernoili out the mas flow as if injector ambient pressure
    %*** calibrate startup with cold flow data
    %of that area for first iteration*****
    delta_Q = state.delta_m_ox_v * oxidizer.delta_H;
    %find change in tank temperature (K)
    state.delta_T_tank = -delta_Q/(state.m_ox_l*oxidizer.Cp);%****say remaining liquid???***
    %find the new oxidizer liquid temperature (K)
    state.T_tank = state.T_tank + state.delta_T_tank; %***seems like temp should equalize ?

elseif state.m_ox_l >= state.m_l_unvaporized && state.m_ox_l > 0%starting to really not likethis method****this might not actually be needed?*****
    %calcualte delta_P (Pa)
    state.delta_P_inj = state.P_tank - state.P_cmbr;
    %calucalte m_dot_oxidizer
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

    %calcualte m_ox_total
    state.m_ox_total = state.m_ox_total - state.m_dot_inj*param.dt - state.m_dot_vent*param.dt;

    delta_P_tank_avg = mean(output.delta_P_tank(1:sum(output.delta_P_tank<0)));
    
    %weird temp ptank only for next step
    P_tank_temp = state.P_tank + delta_P_tank_avg;
    
    %dunno what this does****************************************
    vp = @(T) 7251000*exp((1/(T/309.57))*...
        (-6.71893*(1-T/309.57) + 1.35966*(1-(T/309.57))^(3/2) + -1.3779*...
        (1-(T/309.57))^(5/2) + -4.051*(1-(T/309.57))^5)) - P_tank_temp;
    
    state.T_tank = fzero(vp,state.T_tank);
    
    oxidizer = Nitrous_Properties_Rev1(state.T_tank);
    

    state.m_ox_l = (param.V_tank - (state.m_ox_total/oxidizer.rho_v))/ ...
                    ((1/oxidizer.rho_l)-(1/oxidizer.rho_v));
    

    state.m_l_unvaporized = 0;
else
    
    %calcualte delta_P (Pa)
    state.delta_P_inj = state.P_tank - state.P_cmbr;


    %*****this section is a bandaid for an HRAP bug ***********************
    if state.m_ox_l == 0
        [m_dot,state] = Injector_Mass_Flow_Rev_1(param,state,param.model_inj_v,param.A_inj_total,param.Cd_inj,state.delta_P_inj);
        state.m_dot_inj = m_dot;
    else
        [m_dot,state] = Injector_Mass_Flow_Rev_1(param,state,param.model_inj_l,param.A_inj_total,param.Cd_inj,state.delta_P_inj);
        state.m_dot_inj = m_dot;
    end

    if state.m_ox_l < 0
        state.m_ox_l = 0;
    end
    %******************************************************************
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


    %question about assumption in pg 5 first para******
    %comapre using a textbook blowdown from t atmo*****

    %%%%% initiate varibles for numerical loop solution %%%%%
    %mass at beginng of time step (kg)
    m_ox_1 = state.m_ox_total;
    %calcualte m_ox_total
    state.m_ox_total = state.m_ox_total - state.m_dot_inj*param.dt - state.m_dot_vent*param.dt;
    %mass at end of time step (kg)
    m_ox_2 = state.m_ox_total;
    %init tank temp (K)
    T_1 = state.T_tank;
    %init tank pressure (Pa)
    P_1 = state.P_tank;
    %init Z guess
    Z_guess = oxidizer.Z;
    Z_old = Z_guess;
    
    %find Cp (J/kg*K)
    Cp = CoolProp.PropsSI('CP0MASS','T',293.15,'P',101325,param.ox_name);
    %find Cv (J/kg*K)
    Cv = CoolProp.PropsSI('CVMASS','T',293.15,'P',101325,param.ox_name);
    %find specific heat ratio gamma
    gamma = Cp/Cv;

    %for loop maksure there is a lcause about blowing up in case its unstable
    %temp iteration counter***************
    i = 0;
    delta_Z = 1;
    %numerical vapor phase loop
    while delta_Z >= 0.000001 %can get smaller test sensitivity ******
        %calcualte T2 *** is this supposed to have z old in it
        T_2 = T_1*(((Z_guess*m_ox_2)/(Z_old*m_ox_1))^(0.3));%<--- should be oxidizer.gamma - 1 but not cuz trying to match hrap******
        %calcualte P2
        P_2 = P_1*((T_2/T_1)^(1.3/0.3)); %<--- should be oxidizer.gamma - 1
        %get nitrous properties at T2
        oxidizer = Nitrous_Properties_Rev1(T_2);
        %calculate change in Z
        delta_Z = abs(Z_guess - oxidizer.Z);

        %assing values for next iteration
        Z_guess = (Z_guess + oxidizer.Z)/2;

        %iterate temp i *************
        i = i+1;
        %test if loop has blow up
        if i>1000000
            disp("vapor loop blew up")
            break
        end
    end
    
    %notify user
    disp("vapor loop took")
    disp(i)
    disp("loops to complete")
    
    %get density of vapor
    state.rho_ox_v = oxidizer.rho_v;

    %assign loop results
    state.T_tank = T_2;
    state.P_tank = P_2;

    %perhaps assign null or zero to values that arent used here?****
end

end