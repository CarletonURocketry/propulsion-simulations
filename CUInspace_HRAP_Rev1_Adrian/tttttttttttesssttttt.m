% P_atm = 101325;%pa
% CF3_1_mdot_ideal = [];
% A = 0.00010603;%m^2
% 
% CF3_1_Cd = [];
% 
% % Iterate through each value in m_dot and sum them up
% for i = 1:length(CF3_1_mdot)
%     if i > 40
%         CF3_1_mdot_ideal = [CF3_1_mdot_ideal;0];
%         CF3_1_Cd = [CF3_1_Cd;0];
%     else
%         delta_P = CF3_1_P_tank(i) - P_atm;
%         rho_l = CoolProp.PropsSI('D','P',CF3_1_P_tank(i),'Q',0,'CarbonDioxide');
%         mdot_ideal = A*sqrt(2*rho_l*delta_P);
%         CF3_1_mdot_ideal = [CF3_1_mdot_ideal;mdot_ideal];
%         CF3_1_Cd = [CF3_1_Cd;(CF3_1_mdot(i)/CF3_1_mdot_ideal(i))];
%     end
% 
% end

avg = 0;
num = 0;

for i = 18:33
    avg = avg + CF3_1_mdot(i);
    num = num +1;
end

avg = avg/num;

% hold on
% figure(3)
% plot(CF3_1_time,CF3_1_mdot_ideal)
% figure(5)
% plot(CF3_1_time,CF3_1_Cd)
