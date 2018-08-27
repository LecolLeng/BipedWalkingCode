clc
clear

%% 设定人体的固有参数，身高、体重、平均步速
Human_heigh=1.750; %单位m
Human_weight=70;%单位Kg
Human_Velocity=1.4;%步行平均速度，单位m/s
Human_sex=0;%男 0 女1
%% 依据人体步行时固有参数，推导出
Vr=Human_Velocity/Human_heigh;%相对速度
Gait_frequence=64.8*(Vr^0.57)/60;%人体步频 参考文献：Grieve D W, Gear R J. The relationships between length of stride, step frequency, time of swing and speed of walking for children and adults[J]. Ergonomics, 1966, 9(5): 379-399.
Gait_cycle_time=1/Gait_frequence;%步态周期：
Gait_step_length=Human_Velocity*Gait_cycle_time; %步幅 Stride length，参考文献：
Gait_step_half_length=Gait_step_length/2; %步长，参考文献：
%***********腿长(m)*********参考文献：Roebuck, J.A., Kroemer, K.H.E., Thomson, W.G.: Engineering Anthropometry Methods. Wiley-interscience, New York (1975)
Leg_real_length=0.53*Human_heigh;
%***********脚长(cm)*********参考文献：Agnihotri A K, Purwar B, Googoolye K, et al. Estimation of stature by foot length[J]. Journal of forensic and legal medicine, 2007, 14(5): 279-283.
if Human_sex>0
    Foot_length=(Human_heigh*100-77.059)/3.536;
else
    Foot_length=(Human_heigh*100-68.586)/4.0636;
end
%*****双腿支撑期占整个步态周期的比例*********参考文献：Grieve D W, Gear R J. The relationships between length of stride, step frequency, time of swing and speed of walking for children and adults[J]. Ergonomics, 1966, 9(5): 379-399.
% a=0.07:0.01:0.34;
% num=size(a,2);
% for i=1:num
%   b(i)=(0.5-a(i))*(Human_heigh^0.5)/1.87;
%   tswing(i)=a(i)*Gait_cycle_time+b(i);
%   tstance(i)=1-tswing(i);
%   Double_leg_support_time_rate(i)=tstance(i)-tswing(i);
%
% end
% plot(a,Double_leg_support_time_rate,'-r','LineWidth',2);
% xlabel('a');
% ylabel('Double leg support time rate');
% set(gca,'FontSize',30,'Fontname', 'Times New Roman','LineWidth',1.5);
a=0.205;  
b=(0.5-a)*(Human_heigh^0.5)/1.87;
tswing=a*Gait_cycle_time+b;
tstance=1-tswing;
Double_leg_support_time_rate=tstance-tswing;

%*******人体质心振幅，单位m*****************参考文献：Hoover J, Meguid S A. Performance assessment of the suspended-load backpack[J]. International Journal of Mechanics and Materials in Design, 2011, 7(2): 111-121.
v=Human_Velocity;
Com_amplitude=(Leg_real_length/2)*(1-(1-(0.936*v/(Leg_real_length*2*1.504*(v/Leg_real_length)^0.57))^2)^0.5);


