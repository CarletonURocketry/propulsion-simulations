%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Injector Mass Flow 
% 2024/09/29
% Adrian Comisso
%
% Desccription: 
% Computes injector mass depending on selcted method accoridng to 
% "AN INVESTIGATION OF INJECTORS FOR USE WITH HIGH VAPOR PRESSURE 
% PROPELLANTS WITH APPLICATIONS TO HYBRID ROCKETS"
%
% Inputs:
% param - stores the static parameters of the scenario (eg. tank volume)
% state - stores the current state of each varible
% model - selcted m_dot calcuation model
% Cd - the coeffcicet of discharge of the specific vent or injector
% delta_P - the pressure drop across the specivic vent or injector
%
% Outputs:
% m_dot - mass flow rate through orifice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m_dot,state] = Injector_Mass_Flow(param,state,model,A,Cd,delta_P)

%****maybe just calculate tons of fluid properties here instead of in iff
%statment ot mabe even in thank model before mdot call?*****

%liquid------------------
%compressiblility?
%firction?
%vena contracta?
%Babitsky
%Moody
%Fauske
%Burnell
%Zaloudek
%Proposed
%vapor--------------------
%SPI_v
%IGC ideal gas compressible
%RGC real gas compresible

P_2 = state.P_tank - delta_P;%why is this needed just us p chamber?

if model == "SPI"
    %calcualte m_dot_inj (kg/s)
    m_dot = Cd*A*sqrt(2*state.rho_ox_l*delta_P);
elseif model == "HEM"
    %calculate the tank specific enthalpy of liquid (J/kg)
    h_ox_l = CoolProp.PropsSI('H','T',state.T_tank,'Q',0,param.ox_name);
    %clauclate the tank entropy of liquid (J/kg*K)
    s_ox_l = CoolProp.PropsSI('S','T',state.T_tank,'Q',0,param.ox_name);
    %initial loop coniditions
    m_dot_crit = 0;
    m_dot_old = 0;
    delta_m_dot = 1;
    %if t = 0 find P_2_crit long method otherwise P_2_crit will be less than last time step P_2_crit
    if state.t == 0
        %set initial P_2_crit
        state.P_2_crit = state.P_tank;
        %find P_2_crit (Pa)
        while delta_m_dot >= 0
            %iterate the downstream pressure (Pa)
            state.P_2_crit = state.P_2_crit - 25;
            %calculate the downstream specific enthalpy (J/kg)
            h_2 = CoolProp.PropsSI('H', 'P', state.P_2_crit, 'S', s_ox_l, param.ox_name);
            %clacualte the downstream denstiy (kg/m^3)
            rho_2 = CoolProp.PropsSI('D', 'P', state.P_2_crit, 'S', s_ox_l, param.ox_name);
            %calculate the mass flow rate (kg/s)
            m_dot_crit = Cd * A * rho_2 * sqrt(2 * (h_ox_l - h_2));
            %calucalte change in mass flow rate for this loop (kg/s)
            delta_m_dot = m_dot_crit - m_dot_old;
            %assign m_dot_old (kg/s)
            m_dot_old = m_dot_crit;
        end
    else
        %find P_2_crit (Pa)
        while delta_m_dot >= 0
            %iterate the downstream pressure (Pa)
            state.P_2_crit = state.P_2_crit - 20;
            %calculate the downstream specific enthalpy (J/kg)
            h_2 = CoolProp.PropsSI('H', 'P', state.P_2_crit, 'S', s_ox_l, param.ox_name);
            %clacualte the downstream denstiy (kg/m^3)
            rho_2 = CoolProp.PropsSI('D', 'P', state.P_2_crit, 'S', s_ox_l, param.ox_name);
            %calculate the mass flow rate (kg/s)
            m_dot_crit = Cd * A * rho_2 * sqrt(2 * (h_ox_l - h_2));
            %calucalte change in mass flow rate for this loop (kg/s)
            delta_m_dot = m_dot_crit - m_dot_old;
            %assign m_dot_old (kg/s)
            m_dot_old = m_dot_crit;
        end
    end
    %test if P_2 is below P_2_crit
    if  P_2 <= state.P_2_crit
        %asign choked mass flow:
        m_dot = m_dot_crit;
    else
        %calculate the downstream specific enthalpy (J/kg)
        h_2 = CoolProp.PropsSI('H', 'P', P_2, 'S', s_ox_l, param.ox_name);
        %clacualte the downstream denstiy (kg/m^3)
        rho_2 = CoolProp.PropsSI('D', 'P', P_2, 'S', s_ox_l, param.ox_name);
        %calculate the mass flow rate (kg/s)
        m_dot = Cd * A * rho_2 * sqrt(2 * (h_ox_l - h_2));
    end
elseif model == "Dyer"
    %calcualte m_dot_SPI (kg/s)
    [m_dot_SPI,state] = Injector_Mass_Flow_Rev_1(param,state,'SPI',A,Cd,delta_P);
    %calcualte m_dot_HEM (kg/s)
    [m_dot_HEM,state] = Injector_Mass_Flow_Rev_1(param,state,'HEM',A,Cd,delta_P);
    %calcualte non equilibrium parameter k
    k = sqrt((state.P_tank - P_2)/(state.P_tank - P_2));
    %calcualte m_dot (kg/s)
    m_dot = (k/(1+k))*m_dot_SPI + (1/(k+1))*m_dot_HEM;
elseif model == "SPI_v"
    %calcualte m_dot_inj (kg/s)
    m_dot = Cd*A*sqrt(2*state.rho_ox_v*delta_P);
elseif model == "PGC_v"
    %find Cp (J/kg*K)
    Cp = CoolProp.PropsSI('CP0MASS','T',state.T_tank,'D',state.rho_ox_v,param.ox_name);
    %find Cv (J/kg*K)
    Cv = CoolProp.PropsSI('CVMASS','T',state.T_tank,'D',state.rho_ox_v,param.ox_name);
    %find specific heat ratio gamma
    gamma = Cp/Cv;
    %find P2/P1 pressure ratio
    PR = P_2/state.P_tank;
    %find critical P2/P1 pressure ratio
    PR_choked = (2/(gamma+1))^(gamma/(gamma-1));
    %test if choked
    if PR <= PR_choked
        %calculate m_dot (kg/s)
        m_dot = Cd*A*sqrt(gamma*state.rho_ox_v*state.P_tank*((2/(gamma+1))^((gamma+1)/(gamma-1))));
    else %PR > PR_choked
        %clacualte m_dot (kg/s)
        m_dot = Cd*A*state.rho_ox_v*sqrt(2*Cp*state.T_tank*((PR^(2/gamma))-(PR^((gamma+1)/gamma))));
    end
end
    
end