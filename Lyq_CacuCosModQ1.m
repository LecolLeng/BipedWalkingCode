%%%考虑双腿支撑期，Y方向余弦运动，0到pi-aerfa*pi区间
function Lyq_CacuCosModQ1()
clc
clear all
%定义人体相关参数,分别为1)腿长m、2)质心振幅m、3)步态周期s、4)人体质量kg、5)重力加速度、6)半步长(参考值)、7)最高点质心横向速度、8)双腿支撑起比例
%Subs = [1.6,0.02,1.25,70,9.8,0.7,0.9,0.24  百分比];
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
%Q1区域Legone末端点x方向上初始赋值及变量设置
LoneVxQ10 = Subs(7);
LoneLxQ1=zeros(1,length(ttQ1));
LoneLxQ1(1) = 0;
LoneVxQ1=zeros(1,length(ttQ1));
LoneVxQ1(1)=LoneVxQ10;
LoneFxQ1=zeros(1,length(ttQ1));
LoneAxQ1=zeros(1,length(ttQ1));
%Q1区域Legone末端点Y方向上初始赋值及变量设置
LoneLyQ1=zeros(1,length(ttQ1));
LoneLyQ1(1)=L;
LoneVyQ1=zeros(1,length(ttQ1));
LoneVyQ1(1)=0 ;
LoneFyQ1=zeros(1,length(ttQ1));
LoneFyQ1(1)=-A*w^2*m*cos(w*ttQ10)+m*g+MSbp*g+MLbp*g-FM*cos(w*ttQ10-forcefai);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%弹性背包时需要变化
LoneAyQ1=zeros(1,length(ttQ1));
LoneAyQ1(1)=-A*w^2*cos(w*ttQ10);
%Q1区域Legone末端点xY方向上初始赋值及变量设置
LoneLxyQ1=zeros(1,length(ttQ1));
LoneVxyQ1=zeros(1,length(ttQ1));
LoneFxyQ1=zeros(1,length(ttQ1));
LoneAxyQ1=zeros(1,length(ttQ1));
%Q1区域质心，X坐标系下,质心相应变量设置
ComLxQ1=zeros(1,length(ttQ1));
ComLxQ1(1)=0;
ComVxQ1=zeros(1,length(ttQ1));
ComVxQ1(1)=Subs(7);
ComFxQ1=zeros(1,length(ttQ1));
ComFxQ1(1)=0;
ComAxQ1=zeros(1,length(ttQ1));
ComAxQ1(1)=0;
%Q1区域质心，Y坐标系下,质心相应变量设置
ComLyQ1=zeros(1,length(ttQ1));
ComVyQ1=zeros(1,length(ttQ1));
ComFyQ1=zeros(1,length(ttQ1));
ComAyQ1=zeros(1,length(ttQ1));
%Q1区域质心，XY坐标系下,质心相应变量设置
ComLxyQ1=zeros(1,length(ttQ1));
ComVxyQ1=zeros(1,length(ttQ1));
ComFxyQ1=zeros(1,length(ttQ1));
ComAxyQ1=zeros(1,length(ttQ1));
%%
for i=1:length(ttQ1)  %从0的下一时刻开始
    
    if i>1
        %%
        %Q1区域x方向上的相应计算
        tspan = [0 ttQ1(i)];
        ComLVxQ10 = [LoneLxQ1(1) LoneVxQ1(1)];%x0表示这一段时间初始横向位移，v0表示初始横向速度
        [t,y] = ode113(@(t,y) odefunddQ1(t,y,A,L,w,m,g,MSbp,MLbp,FM,forcefai), tspan, ComLVxQ10);
       
        ComLxQ1(i)=y(end,1);
        ComVxQ1(i)=y(end,2);
        ComFxQ1(i)=(ComVxQ1(i)-ComVxQ1(i-1))/(ttQ1(i)-ttQ1(i-1))*(m+MSbp+MLbp);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%弹性背包时，需要考虑背包质量
        ComAxQ1(i)=ComFxQ1(i)/(m+MSbp+MLbp);
        %Q1区域legone末端点，y方向上的相应计算
        LoneLyQ1(i)=L-A+A*cos(w*ttQ1(i));
        LoneVyQ1(i)=-A*w*sin(w*ttQ1(i));
        LoneFyQ1(i)=-A*w^2*cos(w*ttQ1(i))*m+m*g++MSbp*g+MLbp*g-FM*cos(w*ttQ1(i)-forcefai);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%弹性背包时需要变化;
        LoneAyQ1(i)=-A*w^2*cos(w*ttQ1(i));
    end
    %Q1区域Legone末端点，X坐标系下
    LoneLxQ1(i)=ComLxQ1(i);
    LoneVxQ1(i)=ComVxQ1(i);
    LoneFxQ1(i)=ComFxQ1(i);
    LoneAxQ1(i)=ComAxQ1(i);
    %Q1区域Legone末端点xY方向上
    LoneLxyQ1(i)=sqrt(LoneLxQ1(i)^2+LoneLyQ1(i)^2);
    LoneVxyQ1(i)=sqrt(LoneVxQ1(i)^2+LoneVyQ1(i)^2);
    LoneFxyQ1(i)=sqrt(LoneFxQ1(i)^2+LoneFyQ1(i)^2);
    LoneAxyQ1(i)=sqrt(LoneAxQ1(i)^2+LoneAyQ1(i)^2);
    %Q1区域质心，Y坐标系下
    ComLyQ1(i)=LoneLyQ1(i);
    ComVyQ1(i)=LoneVyQ1(i);
    ComFyQ1(i)=LoneFyQ1(i);
    ComAyQ1(i)=LoneAyQ1(i);
    %Q1区域质心在xY方向上
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
% ylabel('Displacement（m）')
% xlim([0,ttQ1(end)/T])
% figure(2)
% plot(ttQ1/T,LoneVxQ1,'-k',ttQ1/T,LoneVyQ1,'-.r',ttQ1/T,LoneVxyQ1,'--b','LineWidth',2)
% title('Lone-Com Velocity')
% legend('H','V','Com');
% xlabel('Gait cycle')
% ylabel('Velocity（m/s）')
% xlim([0,ttQ1(end)/T])
% figure(3)
% plot(ttQ1/T,LoneFxQ1,'-k',ttQ1/T,LoneFyQ1,'-.r',ttQ1/T,LoneFxyQ1,'--b','LineWidth',2)
% title('Lone-Com Force')
% legend('H','V','Com');
% xlabel('Gait cycle')
% ylabel('Force（N）')
% xlim([0,ttQ1(end)/T])
% figure(4)
% plot(ttQ1/T,LoneAxQ1,'-k',ttQ1/T,LoneAyQ1,'-.r',ttQ1/T,LoneAxyQ1,'--b','LineWidth',2)
% title('Lone-Com Acceleration')
% legend('H','V','Com');
% xlabel('Gait cycle')
% ylabel('Acceleration（m/s^2）')
% xlim([0,ttQ1(end)/T])

save('Q1','ttQ1','Subs','LoneLxQ1','LoneVxQ1','LoneFxQ1','LoneAxQ1','LoneLyQ1','LoneVyQ1','LoneFyQ1','LoneAyQ1','LoneLxyQ1','LoneVxyQ1','LoneFxyQ1','LoneAxyQ1'...
   ,'ComLxQ1','ComVxQ1','ComFxQ1','ComAxQ1','ComLyQ1','ComVyQ1','ComFyQ1','ComAyQ1','ComLxyQ1','ComVxyQ1','ComFxyQ1','ComAxyQ1' );
