%%%考虑双腿支撑期，Y方向余弦运动，0到pi-aerfa*pi区间
function Lyq_CacuCosModQ2()
clc
clear all
%载入人体相关参数,分别为1)腿长m、2)质心振幅m、3)步态周期s、4)人体质量kg、5)重力加速度、6)半步长、7)最高点质心横向速度、8)双腿支撑起比例
%载入Q2区计算结果
load('Q1.mat');
load('par.mat');
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
ttQ200=(1-Q2percent*2)*T/4;
ttQ20=ttQ200-rem(ttQ200,tdeta); 
ttQ2=ttQ20:tdeta:T/4;

%%
%Q2区域Legone,x方向上初始赋值及变量设置
LoneLxQ2=zeros(1,length(ttQ2));
LoneLxQ2(1)=LoneLxQ1(end) ;
LoneVxQ2=zeros(1,length(ttQ2));
LoneVxQ2(1)=LoneVxQ1(end);
LoneFxQ2=zeros(1,length(ttQ2));
LoneFxQ2(1)=LoneFxQ1(end);
LoneAxQ2=zeros(1,length(ttQ2));
LoneAxQ2(1)=LoneAxQ1(end);
%Q2区域Legone,Y方向上初始赋值及变量设置
LoneLyQ2=zeros(1,length(ttQ2));
LoneLyQ2(1)=LoneLyQ1(end);
LoneVyQ2=zeros(1,length(ttQ2));
LoneVyQ2(1)=LoneVyQ1(end) ;
LoneFyQ2=zeros(1,length(ttQ2));
LoneFyQ2(1)=LoneFyQ1(end);
LoneAyQ2=zeros(1,length(ttQ2));
LoneAyQ2(1)=LoneAyQ1(end);
%Q2区域Legone,XY方向上初始赋值及变量设置
LoneLxyQ2=zeros(1,length(ttQ2));
LoneLxyQ2(1)=LoneLxyQ1(end);
LoneVxyQ2=zeros(1,length(ttQ2));
LoneVxyQ2(1)=LoneVxyQ1(end) ;
LoneFxyQ2=zeros(1,length(ttQ2));
LoneFxyQ2(1)=LoneFxyQ1(end);
LoneAxyQ2=zeros(1,length(ttQ2));
LoneAxyQ2(1)=LoneAxyQ1(end);
%Q2区域Legtwo,x方向上初始赋值及变量设置
LtwoLxQ2=zeros(1,length(ttQ2));
LtwoLxQ2(1)=d-LoneLxQ1(end);
LtwoVxQ2=zeros(1,length(ttQ2));
LtwoVxQ2(1)=-LoneVxQ1(end);
LtwoFxQ2=zeros(1,length(ttQ2));
LtwoFxQ2(1)=0;
LtwoAxQ2=zeros(1,length(ttQ2));
LtwoAxQ2(1)=-LoneAxQ1(end);
%Q2区域Legtwo,Y方向上初始赋值及变量设置
LtwoLyQ2=zeros(1,length(ttQ2));
LtwoLyQ2(1)=LoneLyQ1(end);
LtwoVyQ2=zeros(1,length(ttQ2));
LtwoVyQ2(1)=LoneVyQ1(end);
LtwoFyQ2=zeros(1,length(ttQ2));
LtwoFyQ2(1)=0;
LtwoAyQ2=zeros(1,length(ttQ2));
LtwoAyQ2(1)=LoneAyQ1(end);
%Q2区域Legtwo,xY方向上初始赋值及变量设置
LtwoLxyQ2=zeros(1,length(ttQ2));
LtwoVxyQ2=zeros(1,length(ttQ2));
LtwoFxyQ2=zeros(1,length(ttQ2));
LtwoAxyQ2=zeros(1,length(ttQ2));
%Q2区域Legone和Legtwo合力,在X方向上(质心处)
LonetwoFxQ2=zeros(1,length(ttQ2));
LonetwoFxQ2(1)=LoneFxQ1(end);
%Q2区域Legone和Legtwo合力,在Y方向上(质心处)
LonetwoFyQ2=zeros(1,length(ttQ2));
LonetwoFyQ2(1)=LoneFyQ1(end);
%Q2区域质心，X坐标系下,质心相应变量设置
ComLxQ2=zeros(1,length(ttQ2));
ComLxQ2(1)=LoneLxQ2(1);
ComVxQ2=zeros(1,length(ttQ2));
ComVxQ2(1)=LoneVxQ2(1);
ComFxoneQ2=zeros(1,length(ttQ2)); %第一条腿给予质心的x方向力
ComFxtwoQ2=zeros(1,length(ttQ2)); %第二条腿给予质心的x方向力
ComFxQ2=zeros(1,length(ttQ2)); %给予质心的x方向合力
ComAxQ2=zeros(1,length(ttQ2));
%Q2区域质心，X坐标系下,质心相应变量设置
ComLyQ2=zeros(1,length(ttQ2));
ComVyQ2=zeros(1,length(ttQ2));
ComFyoneQ2=zeros(1,length(ttQ2));%第一条腿给予质心的y方向力
ComFytwoQ2=zeros(1,length(ttQ2)); %第二条腿给予质心的y方向力
ComFyQ2=zeros(1,length(ttQ2)); %给予质心的y方向合力
ComAyQ2=zeros(1,length(ttQ2));
%Q2区域质心，XY坐标系下,质心相应变量设置
ComLxyQ2=zeros(1,length(ttQ2));
ComVxyQ2=zeros(1,length(ttQ2));
ComFxyoneQ2=zeros(1,length(ttQ2)); %第一条腿给予质心的xy方向力
ComFxytwoQ2=zeros(1,length(ttQ2)); %第二条腿给予质心的xy方向力
ComFxyQ2=zeros(1,length(ttQ2)); %给予质心的xy方向合力
ComAxyQ2=zeros(1,length(ttQ2));
%%

%%求取变化比例曲线
FLQ2end=MSbp*g+MLbp*g-FM*cos(w*ttQ2(end)-forcefai);
QuaFunParaQ2=FittingCurve(1,ttQ20,LtwoFyQ2(1),ttQ2(end),(Q2percent/(Q2percent+Q3percent))*(-m*A*w^2*cos(w*ttQ2(end))+m*g+FLQ2end));

for i=1:length(ttQ2)  %从0的下一时刻开始
    
    if i>1
        %%
        %Q2区域Legone(除了力),x方向上的相应计算
        tspan = [ttQ20 ttQ2(i)];
        ComLVxQ20 = [ComLxQ2(1) ComVxQ2(1) ];%LoneLxQ20表示这一段时间初始横向位移，LoneLxQ(1)表示初始横向速度
        [t,y] = ode113(@(t,y) odefunddQ2(t,y,A,L,w,g,QuaFunParaQ2(1),QuaFunParaQ2(2),QuaFunParaQ2(3),Q2percent,Q3percent,T,d,m,MSbp,MLbp,FM,forcefai), tspan, ComLVxQ20);
        ComLxQ2(i)=y(end,1);
        ComVxQ2(i)=y(end,2);
        ComFxQ2(i)=(ComVxQ2(i)-ComVxQ2(i-1))/(ttQ2(i)-ttQ2(i-1))*(m+MSbp+MLbp); %(两腿合力)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ComAxQ2(i)=ComFxQ2(i)/(m+MSbp+MLbp); %(两腿合力产生的加速度)
        %Q2区域Legone,x方向上的相应计算
        LoneLxQ2(i)=ComLxQ2(i);
        LoneVxQ2(i)=ComVxQ2(i);
        LoneAxQ2(i)=ComAxQ2(i);
        %Q2区域Legone,y方向上的相应计算
        LoneLyQ2(i)=L-A+A*cos(w*ttQ2(i));
        LoneVyQ2(i)=-A*w*sin(w*ttQ2(i));
        FL=MSbp*g+MLbp*g-FM*cos(w*ttQ2(i)-forcefai); 
        %//修改点
        LtwoFyQ2(i)=QuaFunParaQ2(1)*ttQ2(i)^2+QuaFunParaQ2(2)*ttQ2(i)+QuaFunParaQ2(3);
        LoneFyQ2(i)= -m*A*w^2*cos(w*ttQ2(i))+m*g+FL-LtwoFyQ2(i);
        %LoneFyQ2(i)=-m*A*w^2*cos(w*ttQ2(i))+m*g+FL-(2*ttQ2(i)*(m*A*w^2+m*g+FL)/(aerfa*T))+(1-aerfa)*(m*A*w^2+m*g+FL)/(2*aerfa);
        LoneAyQ2(i)=-A*w^2*cos(w*ttQ2(i));
        %Q2区域Legone（力），x方向上力的计算
        LoneFxQ2(i)=LoneFyQ2(i)*LoneLxQ2(i)/LoneLyQ2(i);
        %Q2区域Legtwo,x方向上的相应计算
        LtwoLxQ2(i)=d-LoneLxQ2(i);
        LtwoVxQ2(i)=-LoneVxQ2(i);
        LtwoFxQ2(i)=LoneFxQ2(i)-ComFxQ2(i); %合力减掉第一条腿的力,(LtwoFxQ2(i)按数值计算)
        LtwoAxQ2(i)=-LoneAxQ2(i);  %加速度为+
        %Q2区域Legtwo,y方向上的相应计算
        LtwoLyQ2(i)= LoneLyQ2(i);
        LtwoVyQ2(i)=LoneVyQ2(i);
        %LtwoFyQ2(i)=;
        %LtwoFyQ2(i)=(2*ttQ2(i)*(m*A*w^2+m*g+FL)/(aerfa*T))-(1-aerfa)*(m*A*w^2+m*g+FL)/(2*aerfa);
        LtwoAyQ2(i)=LoneAyQ2(i);
    end
    %Q2区域Legone,XY方向上初始赋值及变量设置
    LoneLxyQ2(i)=sqrt(LoneLxQ2(i)^2+LoneLyQ2(i)^2);
    LoneVxyQ2(i)=sqrt(LoneVxQ2(i)^2+LoneVyQ2(i)^2);
    LoneFxyQ2(i)=sqrt(LoneFxQ2(i)^2+LoneFyQ2(i)^2);
    LoneAxyQ2(i)=sqrt(LoneAxQ2(i)^2+LoneAyQ2(i)^2);
    %Q2区域Legtwo,xY方向上初始赋值及变量设置
    LtwoLxyQ2(i)=sqrt(LtwoLxQ2(i)^2+LtwoLyQ2(i)^2);
    LtwoVxyQ2(i)=sqrt(LoneVxQ2(i)^2+LtwoVyQ2(i)^2);
    LtwoFxyQ2(i)=sqrt(LtwoFxQ2(i)^2+LtwoFyQ2(i)^2);
    LtwoAxyQ2(i)=sqrt(LtwoAxQ2(i)^2+LtwoAyQ2(i)^2);
    
    %Q2区域质心，X坐标系下,质心相应变量设置
    ComLxQ2(i)=LoneLxQ2(i);
    ComVxQ2(i)=LoneVxQ2(i);
    ComFxoneQ2(i)=LoneFxQ2(i); %第一条腿给予质心的x方向力
    ComFxtwoQ2(i)=LtwoFxQ2(i); %第二条腿给予质心的x方向力
    ComFxQ2(i)=LoneFxQ2(i)-LtwoFxQ2(i);%给予质心的x方向合力
    ComAxQ2(i)=LoneAxQ2(i);
    %Q2区域质心，X坐标系下,质心相应变量设置
    ComLyQ2(i)=LoneLyQ2(i);
    ComVyQ2(i)=LoneVyQ2(i);
    ComFyoneQ2(i)=LoneFyQ2(i);%第一条腿给予质心的y方向力
    ComFytwoQ2(i)=LtwoFyQ2(i); %第二条腿给予质心的y方向力
    ComFyQ2(i)=LoneFyQ2(i)+LtwoFyQ2(i); %给予质心的y方向合力
    ComAyQ2(i)=LoneAyQ2(i);
    %Q2区域质心，XY坐标系下,质心相应变量设置
    ComLxyQ2(i)=sqrt(ComLxQ2(i)^2+ComLyQ2(i)^2);
    ComVxyQ2(i)=sqrt(ComVxQ2(i)^2+ComVyQ2(i)^2);
    ComFxyoneQ2(i)=LoneFxyQ2(i);
    ComFxytwoQ2(i)=LtwoFxyQ2(i);
    ComFxyQ2(i)=sqrt(ComFxQ2(i)^2+ComFyQ2(i)^2);
    ComAxyQ2(i)=LoneAxyQ2(i);
end


%%%步态周期为横坐标
% figure(1)
% plot(ttQ2/T,LoneLxQ2,'-k',ttQ2/T,LoneLyQ2,'-.r',ttQ2/T,LoneLxyQ2,'--b',ttQ2/T,LtwoLxQ2,'-c',ttQ2/T,LtwoLyQ2,'-.m',ttQ2/T,LtwoLxyQ2,'--g','LineWidth',2)
% title('Lone-Com-Ltwo Displacement')
% xlabel('Gait cycle')
% ylabel('Displacement(m)');
% legend('Lone-ComH','Lone-ComV','Lone-Comxy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ2(1)/T,ttQ2(end)/T])
% ylim([0,1.7])
% figure(2)
% plot(ttQ2/T,LoneVxQ2,'-k',ttQ2/T,LoneVyQ2,'-.r',ttQ2/T,LoneVxyQ2,'--b',ttQ2/T,LtwoVxQ2,'-c',ttQ2/T,LtwoVyQ2,'-.m',ttQ2/T,LtwoVxyQ2,'--g','LineWidth',2)
% title('Lone-Com-two Velocity')
% xlabel('Gait cycle')
% ylabel('Velocity(m/s)');
% legend('Lone-ComH','Lone-ComV','Lone-Comxy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ2(1)/T,ttQ2(end)/T])
% figure(3)
% plot(ttQ2/T,LoneFxQ2,'-k',ttQ2/T,LoneFyQ2,'-.r',ttQ2/T,LoneFxyQ2,'--b',ttQ2/T,LtwoFxQ2,'-c',ttQ2/T,LtwoFyQ2,'-.m',ttQ2/T,LtwoFxyQ2,'--g','LineWidth',2)
% title('Lone-two Force')
% xlabel('Gait cycle')
% ylabel('Force(N)');
% legend('LoneH','LoneV','Lonexy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ2(1)/T,ttQ2(end)/T])
% figure(4)
% plot(ttQ2/T,LoneAxQ2,'-k',ttQ2/T,LoneAyQ2,'-.r',ttQ2/T,LoneAxyQ2,'--b',ttQ2/T,LtwoAxQ2,'-c',ttQ2/T,LtwoAyQ2,'-.m',ttQ2/T,LtwoAxyQ2,'--g','LineWidth',2)
% title('Lone-Com-two Acceleration')
% xlabel('Gait cycle')
% ylabel('Acceleration（m/s^2）')
% legend('Lone-ComH','Lone-ComV','Lone-Comxy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ2(1)/T,ttQ2(end)/T])

save('Q2','ttQ2','Subs','LoneLxQ2','LoneVxQ2','LoneFxQ2','LoneAxQ2','LoneLyQ2','LoneVyQ2','LoneFyQ2','LoneAyQ2','LoneLxyQ2','LoneVxyQ2','LoneFxyQ2','LoneAxyQ2',...
                        'LtwoLxQ2','LtwoVxQ2','LtwoFxQ2','LtwoAxQ2','LtwoLyQ2','LtwoVyQ2','LtwoFyQ2','LtwoAyQ2','LtwoLxyQ2','LtwoVxyQ2','LtwoFxyQ2','LtwoAxyQ2',...
                        'ComLxQ2','ComVxQ2','ComFxoneQ2','ComFxtwoQ2','ComFxQ2','ComAxQ2', 'ComLyQ2','ComVyQ2','ComFyoneQ2','ComFytwoQ2','ComFyQ2','ComAyQ2','ComLxyQ2','ComVxyQ2','ComFxyoneQ2','ComFxytwoQ2','ComFxyQ2','ComAxyQ2');
