%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%second attempt at main
%2024/07/30
%Adrian Comisso
%
%Desccription: 
%This code is degined to hyjack HRAPs functions and test
%functionality with the aim of reverse engineering it and creating a
%seperate simultor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% storage structures (to pass to functions) %%%%%
%inputs
i = struct();
%current state (durnig sim)
x = struct(); 
%output
o = struct();
%units ********************************
u = struct(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% tank %%%%%
i.tnk_V = 0.00585;
x.T_tnk = 293.15;
x.ox_props                  = NOX(x.T_tnk);
x.m_o                   = 4.31;
x.P_tnk                     = x.ox_props.Pv;
x.mdot_o                    = 0;
x.mLiq_new                  = (i.tnk_V - (x.m_o/x.ox_props.rho_v))/((1/x.ox_props.rho_l)-(1/x.ox_props.rho_v));
x.mLiq_old                  = x.mLiq_new + 1;
%%%%% injector %%%%%
i.inj_CdA = 0.2*(0.25*pi*(0.0015)^2);
i.inj_N = 60;
%%%%% camber %%%%%
i.cmbr_V = 0.0098;
i.grn_OD = 0.1016;
i.grn_ID = 0.0711;
i.grn_L = 0.9144;
i.cstar_eff = 0.8;

i.regression_model = @(s,x) shift_OF(s,x);
i.prop_Reg(1) = 0.0876;
i.prop_Reg(2) = 0.3953;
i.prop_Reg(3) = 0;
i.prop_Rho = 900;
i.const_OF = 1.5;
x.P_cmbr                    = 101325;
x.OF                        = 0;
x.mdot_f                    = 0;
x.mdot_n                    = 0;
x.rdot                      = 0;
x.grn_ID                    = i.grn_ID;

path = pwd;
file = 'Paraffin_Rev1.mat';
L = importdata(fullfile(path, file));
i.prop_nm = L.prop_nm; %#ok<*ADPROPLC> 
i.prop_k = L.prop_k;
i.prop_M = L.prop_M;
i.prop_OF = L.prop_OF;
i.prop_Pc = L.prop_Pc;
i.prop_T = L.prop_T;
i.opt_OF = L.opt_OF;
%%%%% nozzle %%%%%
i.noz_Cd = 0.9;
i.noz_thrt = 0.0348;
i.noz_ER = 4.0293;
i.noz_eff = 0.95;
%%%%% center of mass %%%%%
x.m_f                       = 0.25*pi*(i.grn_OD^2 - i.grn_ID^2)*i.prop_Rho*i.grn_L;
x.m_g                       = 1.225*(i.cmbr_V - 0.25*pi*(i.grn_OD^2 - i.grn_ID^2)*i.grn_L);
i.mp_calc = 0;
%%%%% other %%%%%
i.mtr_nm = "a";
i.tmax = 20;
i.tburn = 20;
i.dt = 0.001;
i.Pa = 101325;
t                           = 0;

i.vnt_S = 0;



o.t                         = zeros(1,i.tmax/i.dt + 1);
o.m_o                       = zeros(1,i.tmax/i.dt + 1);
o.P_tnk                     = zeros(1,i.tmax/i.dt + 1);
o.P_cmbr                    = zeros(1,i.tmax/i.dt + 1);
o.mdot_o                    = zeros(1,i.tmax/i.dt + 1);
o.mdot_f                    = zeros(1,i.tmax/i.dt + 1);
o.OF                        = zeros(1,i.tmax/i.dt + 1);
o.grn_ID                    = zeros(1,i.tmax/i.dt + 1);
o.mdot_n                    = zeros(1,i.tmax/i.dt + 1);
o.rdot                      = zeros(1,i.tmax/i.dt + 1);
o.m_f                       = zeros(1,i.tmax/i.dt + 1);
o.F_thr                     = zeros(1,i.tmax/i.dt + 1);
o.dP                        = zeros(1,i.tmax/i.dt + 1);

o.m_o(1)                    = x.m_o;
o.P_tnk(1)                  = x.P_tnk;
o.P_cmbr(1)                 = x.P_cmbr;
o.mdot_o(1)                 = x.mdot_o;
o.mdot_f(1)                 = x.mdot_f;
o.OF(1)                     = x.OF;
o.grn_ID(1)                 = x.grn_ID;
o.mdot_n(1)                 = x.mdot_n;
o.rdot(1)                   = x.rdot;
o.m_f(1)                    = x.m_f;



if i.mp_calc == 1
    o.m_t                   = zeros(1,i.tmax/i.dt + 1);
    o.cg                    = zeros(1,i.tmax/i.dt + 1);

    mp                          = mass(i,x);

    o.m_t(1)                = mp(1);
    o.cg(1)                 = mp(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[i,x,o,t] = sim_loop(i,x,o,t);

plot(o.t,o.P_cmbr)
disp(o.sim_end_cond)

%clear o i x u t



