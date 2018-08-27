%%%����˫��֧���ڣ�Y���������˶���0��pi-aerfa*pi����
function Lyq_CacuCosModQ1()
clc
clear all
%����������ز���,�ֱ�Ϊ1)�ȳ�m��2)�������m��3)��̬����s��4)��������kg��5)�������ٶȡ�6)�벽��(�ο�ֵ)��7)��ߵ����ĺ����ٶȡ�8)˫��֧�������
%Subs = [1.6,0.02,1.25,70,9.8,0.7,0.9,0.24  �ٷֱ�];
load('par');
L = Subs(1);
A = Subs(2);
T = Subs(3);
w = 2*2*pi/T;
m = Subs(4);
g = Subs(5);
d = Subs(6);
aerfa=Subs(8);
tdeta=Subs(9);
Q2percent=Subs(10);
ttQ10=0;
ttQ1=0:tdeta:((1-Q2percent*2)*T/4);

%%
%Q1����Legoneĩ�˵�x�����ϳ�ʼ��ֵ����������
LoneVxQ10 = Subs(7);
LoneLxQ1=zeros(1,length(ttQ1));
LoneLxQ1(1) = 0;
LoneVxQ1=zeros(1,length(ttQ1));
LoneVxQ1(1)=LoneVxQ10;
LoneFxQ1=zeros(1,length(ttQ1));
LoneAxQ1=zeros(1,length(ttQ1));
%Q1����Legoneĩ�˵�Y�����ϳ�ʼ��ֵ����������
LoneLyQ1=zeros(1,length(ttQ1));
LoneLyQ1(1)=L;
LoneVyQ1=zeros(1,length(ttQ1));
LoneVyQ1(1)=0 ;
LoneFyQ1=zeros(1,length(ttQ1));
LoneFyQ1(1)=-A*w^2*m*cos(w*ttQ10)+m*g+MSbp*g+MLbp*g-FM*cos(w*ttQ10-forcefai);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���Ա���ʱ��Ҫ�仯
LoneAyQ1=zeros(1,length(ttQ1));
LoneAyQ1(1)=-A*w^2*cos(w*ttQ10);
%Q1����Legoneĩ�˵�xY�����ϳ�ʼ��ֵ����������
LoneLxyQ1=zeros(1,length(ttQ1));
LoneVxyQ1=zeros(1,length(ttQ1));
LoneFxyQ1=zeros(1,length(ttQ1));
LoneAxyQ1=zeros(1,length(ttQ1));
%Q1�������ģ�X����ϵ��,������Ӧ��������
ComLxQ1=zeros(1,length(ttQ1));
ComLxQ1(1)=0;
ComVxQ1=zeros(1,length(ttQ1));
ComVxQ1(1)=Subs(7);
ComFxQ1=zeros(1,length(ttQ1));
ComFxQ1(1)=0;
ComAxQ1=zeros(1,length(ttQ1));
ComAxQ1(1)=0;
%Q1�������ģ�Y����ϵ��,������Ӧ��������
ComLyQ1=zeros(1,length(ttQ1));
ComVyQ1=zeros(1,length(ttQ1));
ComFyQ1=zeros(1,length(ttQ1));
ComAyQ1=zeros(1,length(ttQ1));
%Q1�������ģ�XY����ϵ��,������Ӧ��������
ComLxyQ1=zeros(1,length(ttQ1));
ComVxyQ1=zeros(1,length(ttQ1));
ComFxyQ1=zeros(1,length(ttQ1));
ComAxyQ1=zeros(1,length(ttQ1));
%%
for i=1:length(ttQ1)  %��0����һʱ�̿�ʼ
    
    if i>1
        %%
        %Q1����x�����ϵ���Ӧ����
        tspan = [0 ttQ1(i)];
        ComLVxQ10 = [LoneLxQ1(1) LoneVxQ1(1)];%x0��ʾ��һ��ʱ���ʼ����λ�ƣ�v0��ʾ��ʼ�����ٶ�
        [t,y] = ode113(@(t,y) odefunddQ1(t,y,A,L,w,m,g,MSbp,MLbp,FM,forcefai), tspan, ComLVxQ10);
       
        ComLxQ1(i)=y(end,1);
        ComVxQ1(i)=y(end,2);
        ComFxQ1(i)=(ComVxQ1(i)-ComVxQ1(i-1))/(ttQ1(i)-ttQ1(i-1))*(m+MSbp+MLbp);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���Ա���ʱ����Ҫ���Ǳ�������
        ComAxQ1(i)=ComFxQ1(i)/(m+MSbp+MLbp);
        %Q1����legoneĩ�˵㣬y�����ϵ���Ӧ����
        LoneLyQ1(i)=L-A+A*cos(w*ttQ1(i));
        LoneVyQ1(i)=-A*w*sin(w*ttQ1(i));
        LoneFyQ1(i)=-A*w^2*cos(w*ttQ1(i))*m+m*g++MSbp*g+MLbp*g-FM*cos(w*ttQ1(i)-forcefai);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���Ա���ʱ��Ҫ�仯;
        LoneAyQ1(i)=-A*w^2*cos(w*ttQ1(i));
    end
    %Q1����Legoneĩ�˵㣬X����ϵ��
    LoneLxQ1(i)=ComLxQ1(i);
    LoneVxQ1(i)=ComVxQ1(i);
    LoneFxQ1(i)=ComFxQ1(i);
    LoneAxQ1(i)=ComAxQ1(i);
    %Q1����Legoneĩ�˵�xY������
    LoneLxyQ1(i)=sqrt(LoneLxQ1(i)^2+LoneLyQ1(i)^2);
    LoneVxyQ1(i)=sqrt(LoneVxQ1(i)^2+LoneVyQ1(i)^2);
    LoneFxyQ1(i)=sqrt(LoneFxQ1(i)^2+LoneFyQ1(i)^2);
    LoneAxyQ1(i)=sqrt(LoneAxQ1(i)^2+LoneAyQ1(i)^2);
    %Q1�������ģ�Y����ϵ��
    ComLyQ1(i)=LoneLyQ1(i);
    ComVyQ1(i)=LoneVyQ1(i);
    ComFyQ1(i)=LoneFyQ1(i);
    ComAyQ1(i)=LoneAyQ1(i);
    %Q1����������xY������
    ComLxyQ1(i)=LoneLxyQ1(i);
    ComVxyQ1(i)=LoneVxyQ1(i);
    ComFxyQ1(i)=LoneFxyQ1(i);
    ComAxyQ1(i)=LoneAxyQ1(i);
    
end
% 
% figure(1)
% plot(ttQ1/T,LoneLxQ1,'-k',ttQ1/T,LoneLyQ1,'-.r',ttQ1/T,LoneLxyQ1,'--b','LineWidth',2)
% title('Lone-Com Displacement')
% legend('H','V','Com');
% xlabel('Gait cycle')
% ylabel('Displacement��m��')
% xlim([0,ttQ1(end)/T])
% figure(2)
% plot(ttQ1/T,LoneVxQ1,'-k',ttQ1/T,LoneVyQ1,'-.r',ttQ1/T,LoneVxyQ1,'--b','LineWidth',2)
% title('Lone-Com Velocity')
% legend('H','V','Com');
% xlabel('Gait cycle')
% ylabel('Velocity��m/s��')
% xlim([0,ttQ1(end)/T])
% figure(3)
% plot(ttQ1/T,LoneFxQ1,'-k',ttQ1/T,LoneFyQ1,'-.r',ttQ1/T,LoneFxyQ1,'--b','LineWidth',2)
% title('Lone-Com Force')
% legend('H','V','Com');
% xlabel('Gait cycle')
% ylabel('Force��N��')
% xlim([0,ttQ1(end)/T])
% figure(4)
% plot(ttQ1/T,LoneAxQ1,'-k',ttQ1/T,LoneAyQ1,'-.r',ttQ1/T,LoneAxyQ1,'--b','LineWidth',2)
% title('Lone-Com Acceleration')
% legend('H','V','Com');
% xlabel('Gait cycle')
% ylabel('Acceleration��m/s^2��')
% xlim([0,ttQ1(end)/T])

save('Q1','ttQ1','Subs','LoneLxQ1','LoneVxQ1','LoneFxQ1','LoneAxQ1','LoneLyQ1','LoneVyQ1','LoneFyQ1','LoneAyQ1','LoneLxyQ1','LoneVxyQ1','LoneFxyQ1','LoneAxyQ1'...
   ,'ComLxQ1','ComVxQ1','ComFxQ1','ComAxQ1','ComLyQ1','ComVyQ1','ComFyQ1','ComAyQ1','ComLxyQ1','ComVxyQ1','ComFxyQ1','ComAxyQ1' );
