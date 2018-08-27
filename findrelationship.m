clc
clear

%% �趨����Ĺ��в�������ߡ����ء�ƽ������
Human_heigh=1.750; %��λm
Human_weight=70;%��λKg
Human_Velocity=1.4;%����ƽ���ٶȣ���λm/s
Human_sex=0;%�� 0 Ů1
%% �������岽��ʱ���в������Ƶ���
Vr=Human_Velocity/Human_heigh;%����ٶ�
Gait_frequence=64.8*(Vr^0.57)/60;%���岽Ƶ �ο����ף�Grieve D W, Gear R J. The relationships between length of stride, step frequency, time of swing and speed of walking for children and adults[J]. Ergonomics, 1966, 9(5): 379-399.
Gait_cycle_time=1/Gait_frequence;%��̬���ڣ�
Gait_step_length=Human_Velocity*Gait_cycle_time; %���� Stride length���ο����ף�
Gait_step_half_length=Gait_step_length/2; %�������ο����ף�
%***********�ȳ�(m)*********�ο����ף�Roebuck, J.A., Kroemer, K.H.E., Thomson, W.G.: Engineering Anthropometry Methods. Wiley-interscience, New York (1975)
Leg_real_length=0.53*Human_heigh;
%***********�ų�(cm)*********�ο����ף�Agnihotri A K, Purwar B, Googoolye K, et al. Estimation of stature by foot length[J]. Journal of forensic and legal medicine, 2007, 14(5): 279-283.
if Human_sex>0
    Foot_length=(Human_heigh*100-77.059)/3.536;
else
    Foot_length=(Human_heigh*100-68.586)/4.0636;
end
%*****˫��֧����ռ������̬���ڵı���*********�ο����ף�Grieve D W, Gear R J. The relationships between length of stride, step frequency, time of swing and speed of walking for children and adults[J]. Ergonomics, 1966, 9(5): 379-399.
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

%*******���������������λm*****************�ο����ף�Hoover J, Meguid S A. Performance assessment of the suspended-load backpack[J]. International Journal of Mechanics and Materials in Design, 2011, 7(2): 111-121.
v=Human_Velocity;
Com_amplitude=(Leg_real_length/2)*(1-(1-(0.936*v/(Leg_real_length*2*1.504*(v/Leg_real_length)^0.57))^2)^0.5);


