%%%考虑双腿支撑期，Y方向余弦运动，0到pi-aerfa*pi区间
function Lyq_CacuCosModQ4()
clc
clear all
%载入人体相关参数,分别为1)腿长m、2)质心振幅m、3)步态周期s、4)人体质量kg、5)重力加速度、6)半步长、7)最高点质心横向速度、8)双腿支撑起比例
%载入Q3区计算结果
load('par.mat');
load('Q1.mat');
load('Q2.mat');
load('Q3.mat');
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
Q3percent=aerfa-Q2percent;
ttQ400=(1+Q3percent*2)*T/4;
ttQ40=ttQ400-rem(ttQ400,tdeta); 
ttQ4=ttQ40:tdeta:T/2;

%%



%Q4区域x方向上初始赋值及变量设置
LtwoLxQ4=zeros(1,length(ttQ4));
LtwoLxQ4(1)=LtwoLxQ3(end);
LtwoVxQ4=zeros(1,length(ttQ4));
LtwoVxQ4(1)=LtwoVxQ3(end);
LtwoFxQ4=zeros(1,length(ttQ4));
LtwoFxQ4(1)=LtwoFxQ3(end);
LtwoAxQ4=zeros(1,length(ttQ4));
LtwoAxQ4(1)=LtwoAxQ3(end);
%Q4区域Y方向上初始赋值及变量设置
LtwoLyQ4=zeros(1,length(ttQ4));
LtwoLyQ4(1)=LtwoLyQ3(end);
LtwoVyQ4=zeros(1,length(ttQ4));
LtwoVyQ4(1)=LtwoVyQ3(end);
LtwoFyQ4=zeros(1,length(ttQ4));
LtwoFyQ4(1)=LtwoFyQ3(end);
LtwoAyQ4=zeros(1,length(ttQ4));
LtwoAyQ4(1)=LtwoAyQ3(end);
%Q4区域Legtwo,xY方向上初始赋值及变量设置
LtwoLxyQ4=zeros(1,length(ttQ4));
LtwoLxyQ4(1)=LtwoLxyQ3(end);
LtwoVxyQ4=zeros(1,length(ttQ4));
LtwoVxyQ4(1)=LtwoVxyQ3(end);
LtwoFxyQ4=zeros(1,length(ttQ4));
LtwoFxyQ4(1)=LtwoFxyQ3(end);
LtwoAxyQ4=zeros(1,length(ttQ4));
LtwoAxyQ4(1)=LtwoAxyQ3(end);
%Q4区域质心,x方向,质心相应变量设置
ComLxQ4=zeros(1,length(ttQ4));
ComLxQ4(1)=ComLxQ3(end);
ComVxQ4=zeros(1,length(ttQ4));
ComVxQ4(1)=ComVxQ3(end);
ComFxQ4=zeros(1,length(ttQ4));
ComFxQ4(1)=ComFxQ3(end);
ComAxQ4=zeros(1,length(ttQ4));
ComAxQ4(1)=LoneAxQ3(end);
%Q4区域质心,y方向,质心相应变量设置
ComLyQ4=zeros(1,length(ttQ4));
ComLyQ4(1)=ComLyQ3(end);
ComVyQ4=zeros(1,length(ttQ4));
ComVyQ4(1)=ComVyQ3(end);
ComFyQ4=zeros(1,length(ttQ4));
ComFyQ4(1)=ComFyQ3(end);
ComAyQ4=zeros(1,length(ttQ4));
ComAyQ4(1)=LoneAyQ3(end);
%Q4区域质心,xy方向,质心相应变量设置
ComLxyQ4=zeros(1,length(ttQ4));
ComVxyQ4=zeros(1,length(ttQ4));
ComFxyQ4=zeros(1,length(ttQ4));
ComAxyQ4=zeros(1,length(ttQ4));
%%
for i=1:length(ttQ4)  %从0的下一时刻开始
    
    if i>1
        %%
        %Q4区域x方向上的相应计算
        tspan = [ttQ40 ttQ4(i)];
        ComLVxQ40 = [ ComLxQ4(1) ComVxQ4(1)];%表示这一段时间初始横向位移，表示初始横向速度
        [t,y] = ode113(@(t,y) odefunddQ4(t,y,A,L,w,g,d,T,m,MSbp,MLbp,FM,forcefai), tspan, ComLVxQ40);
        ComLxQ4(i)=y(end,1);
        ComVxQ4(i)=y(end,2);
        ComFxQ4(i)=(ComVxQ4(i)-ComVxQ4(i-1))/(ttQ4(i)-ttQ4(i-1))*(m+MSbp+MLbp);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ComAxQ4(i)=(ComVxQ4(i)-ComVxQ4(i-1))/(ttQ4(i)-ttQ4(i-1));
        %Q4区域质心，y方向上的相应计算
        ComLyQ4(i)=L-A+A*cos(w*ttQ4(i));
        ComVyQ4(i)=-A*w*sin(w*ttQ4(i));
        FL=MSbp*g+MLbp*g-FM*cos(w*ttQ4(i)-forcefai);
        ComFyQ4(i)=-A*w^2*cos(w*ttQ4(i))*m+m*g+FL;
        ComAyQ4(i)=-A*w^2*cos(w*ttQ4(i));
        %Q4区域Ltwo,x方向上的相应计算
        LtwoLxQ4(i)=d-ComLxQ4(i);
        LtwoVxQ4(i)=-ComVxQ4(i);
        LtwoFxQ4(i)=-ComFxQ4(i);
        LtwoAxQ4(i)=-ComAxQ4(i);
        %Q4区域Ltwo,y方向上的相应计算
        LtwoLyQ4(i)=ComLyQ4(i);
        LtwoVyQ4(i)=ComVyQ4(i);
        LtwoFyQ4(i)=ComFyQ4(i);
        LtwoAyQ4(i)=ComAyQ4(i);
    end
    %Q4区域Ltwo,xy方向上的相应计算
    LtwoLxyQ4(i)=sqrt(LtwoLxQ4(i)^2+LtwoLyQ4(i)^2);
    LtwoVxyQ4(i)=sqrt(LtwoVxQ4(i)^2+LtwoVyQ4(i)^2);
    LtwoFxyQ4(i)=sqrt(LtwoFxQ4(i)^2+ComFyQ4(i)^2);
    LtwoAxyQ4(i)=sqrt(LtwoAxQ4(i)^2+LtwoAyQ4(i)^2);
    %Q4区域质心，xy方向上的相应计算
    ComLxyQ4(i)=sqrt(ComLxQ4(i)^2+ComLyQ4(i)^2);
    ComVxyQ4(i)=sqrt(ComVxQ4(i)^2+ComVyQ4(i)^2);
    ComFxyQ4(i)=sqrt(ComFxQ4(i)^2+ComFyQ4(i)^2);
    ComAxyQ4(i)=sqrt(ComAxQ4(i)^2+ComAyQ4(i)^2);
end

%%

% figure(1)
% plot(ttQ4/T,ComLxQ4,'-k',ttQ4/T,ComLyQ4,'-.r',ttQ4/T,ComLxyQ4,'--b',ttQ4/T,LtwoLxQ4,'-c',ttQ4/T,LtwoLyQ4,'-.m',ttQ4/T,LtwoLxyQ4,'--g','LineWidth',2)
% title('Com-Ltwo Displacement')
% xlabel('Gait cycle')
% ylabel('Displacement(m)');
% legend('ComH','ComV','Comxy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ4(1)/T,ttQ4(end)/T])
% ylim([0,1.7])
% figure(2)
% plot(ttQ4/T,ComVxQ4,'-k',ttQ4/T,ComVyQ4,'-.r',ttQ4/T,ComVxyQ4,'--b',ttQ4/T,LtwoVxQ4,'-c',ttQ4/T,LtwoVyQ4,'-.m',ttQ4/T,LtwoVxyQ4,'--g','LineWidth',2)
% title('Com-Ltwo Velocity')
% xlabel('Gait cycle')
% ylabel('Velocity(m/s)');
% legend('ComH','ComV','Comxy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ4(1)/T,ttQ4(end)/T])
% figure(3)
% plot(ttQ4/T,ComFxQ4,'-k',ttQ4/T,ComFyQ4,'-.r',ttQ4/T,ComFxyQ4,'--b',ttQ4/T,LtwoFxQ4,'-k',ttQ4/T,LtwoFyQ4,'-.r',ttQ4/T,LtwoFxyQ4,'--b','LineWidth',2)
% title('Com-Ltwo Force')
% xlabel('Gait cycle')
% ylabel('Force(N)');
% legend('ComH','ComV','Comxy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ4(1)/T,ttQ4(end)/T])
% figure(4)
% plot(ttQ4/T,ComAxQ4,'-k',ttQ4/T,ComAyQ4,'-.r',ttQ4/T,ComAxyQ4,'--b',ttQ4/T,LtwoAxQ4,'-c',ttQ4/T,LtwoAyQ4,'-.m',ttQ4/T,LtwoAxyQ4,'--g','LineWidth',2)
% title('Com-Ltwo Acceleration')
% xlabel('Gait cycle')
% ylabel('Acceleration（m/s^2）')
% legend('ComH','ComV','Comxy','LtwoH','LtwoV','Ltwoxy');

save('Q4','ttQ4','Subs','LtwoLxQ4','LtwoVxQ4','LtwoFxQ4','LtwoAxQ4','LtwoLyQ4','LtwoVyQ4','LtwoFyQ4','LtwoAyQ4','LtwoLxyQ4','LtwoVxyQ4','LtwoFxyQ4','LtwoAxyQ4',...
    'ComLxQ4','ComVxQ4','ComFxQ4','ComAxQ4','ComLyQ4','ComVyQ4','ComFyQ4','ComAyQ4','ComLxyQ4','ComVxyQ4','ComFxyQ4','ComAxyQ4');
