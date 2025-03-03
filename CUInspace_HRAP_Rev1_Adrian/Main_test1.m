%atteptt at main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = struct(); % input structure
o = struct(); % output structure
x = struct(); % state structure
u = struct(); % selected units structure

s.grn_OD = 0.1016;
s.grn_ID = 0.0711;
s.grn_L = 0.9144;
s.cstar_eff = 0.8;

s.mtr_nm = "a";
s.tmax = 20;
s.tburn = 15;
s.dt = 0.001;
s.Pa = 101325;

s.tnk_V = 0.00585;

s.cmbr_V = 0.0098;

s.regression_model = @(s,x) shift_OF(s,x);
s.prop_Reg(1) = 0.0876;%0.12;
s.prop_Reg(2) = 0.3953;%0.49;
s.prop_Reg(3) = 0;

s.noz_Cd = 0.9;
s.noz_thrt = 0.0348;
s.noz_ER = 4.0293;

s.noz_eff = 0.95;
s.inj_CdA = 0.25*pi*(0.0015)^2*0.2;
s.inj_N = 60;

s.vnt_S = 0;

s.mp_calc = 0;

x.T_tnk = 285.52;

x.ox_props                  = NOX(x.T_tnk);
s.prop_Rho = 900;
x.m_o                   = 4.309;

x.P_tnk                     = x.ox_props.Pv;
x.P_cmbr                    = 101325;
x.mdot_o                    = 0;
x.mLiq_new                  = (s.tnk_V - (x.m_o/x.ox_props.rho_v))/((1/x.ox_props.rho_l)-(1/x.ox_props.rho_v));
x.mLiq_old                  = x.mLiq_new + 1;
x.m_f                       = 0.25*pi*(s.grn_OD^2 - s.grn_ID^2)*s.prop_Rho*s.grn_L;
x.m_g                       = 1.225*(s.cmbr_V - 0.25*pi*(s.grn_OD^2 - s.grn_ID^2)*s.grn_L);

x.OF                        = 0;

x.mdot_f                    = 0;
x.mdot_n                    = 0;
x.rdot                      = 0;
x.grn_ID                    = s.grn_ID;

s.prop_file = "C:\Users\adric\Desktop\4-CUInSpace\Hybrid Sim\HRAP_Decompiled_Matlab_Versions\CUInspace_HRAP_Rev1\propellant_configs\Paraffin.mat";
s.prop_nm = L.prop_nm; %#ok<*ADPROPLC> 
s.prop_k = L.prop_k;
s.prop_M = L.prop_M;
s.prop_OF = L.prop_OF;
s.prop_Pc = L.prop_Pc;
%s.prop_Reg = L.prop_Reg;
%s.prop_Rho = L.prop_Rho;
s.prop_T = L.prop_T;
s.opt_OF = L.opt_OF;

t                           = 0;

o.t                         = zeros(1,s.tmax/s.dt + 1);
o.m_o                       = zeros(1,s.tmax/s.dt + 1);
o.P_tnk                     = zeros(1,s.tmax/s.dt + 1);
o.P_cmbr                    = zeros(1,s.tmax/s.dt + 1);
o.mdot_o                    = zeros(1,s.tmax/s.dt + 1);
o.mdot_f                    = zeros(1,s.tmax/s.dt + 1);
o.OF                        = zeros(1,s.tmax/s.dt + 1);
o.grn_ID                    = zeros(1,s.tmax/s.dt + 1);
o.mdot_n                    = zeros(1,s.tmax/s.dt + 1);
o.rdot                      = zeros(1,s.tmax/s.dt + 1);
o.m_f                       = zeros(1,s.tmax/s.dt + 1);
o.F_thr                     = zeros(1,s.tmax/s.dt + 1);
o.dP                        = zeros(1,s.tmax/s.dt + 1);

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



if s.mp_calc == 1
    o.m_t                   = zeros(1,s.tmax/s.dt + 1);
    o.cg                    = zeros(1,s.tmax/s.dt + 1);

    mp                          = mass(s,x);

    o.m_t(1)                = mp(1);
    o.cg(1)                 = mp(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[s,x,o,t] = sim_loop(s,x,o,t);

plot(o.t,o.F_thr)

clear o s x u t



