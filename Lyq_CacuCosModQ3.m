%%%考虑双腿支撑期，Y方向余弦运动，0到pi-aerfa*pi区间
function Lyq_CacuCosModQ3()
clc
clear all
%载入人体相关参数,分别为1)腿长m、2)质心振幅m、3)步态周期s、4)人体质量kg、5)重力加速度、6)半步长、7)最高点质心横向速度、8)双腿支撑起比例
%载入Q3区计算结果
load('par.mat');
load('Q1.mat');
load('Q2.mat');
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
ttQ300=T/4;
ttQ30=ttQ300-rem(ttQ300,tdeta); 
ttQ3=ttQ30:tdeta:((1+Q3percent*2)*T/4);

%%
%Q3区域Legone,x方向上初始赋值及变量设置
LoneLxQ3=zeros(1,length(ttQ3));
LoneLxQ3(1)=LoneLxQ2(end);
LoneVxQ3=zeros(1,length(ttQ3));
LoneVxQ3(1)=LoneVxQ2(end);
LoneFxQ3=zeros(1,length(ttQ3));
LoneFxQ3(1)=LoneFxQ2(end);
LoneAxQ3=zeros(1,length(ttQ3));
LoneAxQ3(1)=LoneAxQ2(end);
%Q3区域Legone,Y方向上初始赋值及变量设置
LoneLyQ3=zeros(1,length(ttQ3));
LoneLyQ3(1)=LoneLyQ2(end);
LoneVyQ3=zeros(1,length(ttQ3));
LoneVyQ3(1)=LoneVyQ2(end) ;
LoneFyQ3=zeros(1,length(ttQ3));
LoneFyQ3(1)=LoneFyQ2(end);
LoneAyQ3=zeros(1,length(ttQ3));
LoneAyQ3(1)=LoneAyQ2(end);
%Q3区域Legone,XY方向上初始赋值及变量设置
LoneLxyQ3=zeros(1,length(ttQ3));
LoneLxyQ3(1)=LoneLxyQ2(end);
LoneVxyQ3=zeros(1,length(ttQ3));
LoneVxyQ3(1)=LoneVxyQ2(end) ;
LoneFxyQ3=zeros(1,length(ttQ3));
LoneFxyQ3(1)=LoneFxyQ2(end);
LoneAxyQ3=zeros(1,length(ttQ3));
LoneAxyQ3(1)=LoneAxyQ2(end);
%Q3区域Legtwo,x方向上初始赋值及变量设置
LtwoLxQ3=zeros(1,length(ttQ3));
LtwoLxQ3(1)=LtwoLxQ2(end);
LtwoVxQ3=zeros(1,length(ttQ3));
LtwoVxQ3(1)=LtwoVxQ2(end);
LtwoFxQ3=zeros(1,length(ttQ3));
LtwoFxQ3(1)=LtwoFxQ2(end);
LtwoAxQ3=zeros(1,length(ttQ3));
LtwoAxQ3(1)=LtwoAxQ2(end);
%Q3区域Legtwo,Y方向上初始赋值及变量设置
LtwoLyQ3=zeros(1,length(ttQ3));
LtwoLyQ3(1)=LtwoLyQ2(end);
LtwoVyQ3=zeros(1,length(ttQ3));
LtwoVyQ3(1)=LtwoVyQ2(end);
LtwoFyQ3=zeros(1,length(ttQ3));
LtwoFyQ3(1)=LtwoFyQ2(end);
LtwoAyQ3=zeros(1,length(ttQ3));
LtwoAyQ3(1)=LtwoAyQ2(end);
%Q3区域Legtwo,xY方向上初始赋值及变量设置
LtwoLxyQ3=zeros(1,length(ttQ3));
LtwoLxyQ3(1)=LtwoLxyQ2(end);
LtwoVxyQ3=zeros(1,length(ttQ3));
LtwoVxyQ3(1)=LtwoVxyQ2(end);
LtwoFxyQ3=zeros(1,length(ttQ3));
LtwoFxyQ3(1)=LtwoFxyQ2(end);
LtwoAxyQ3=zeros(1,length(ttQ3));
LtwoAxyQ3(1)=LtwoAxyQ2(end);

%Q3区域Legone和Legtwo合力,在X方向上(质心处)
LonetwoFxQ3=zeros(1,length(ttQ3));
LonetwoFxQ3(1)=LoneFxQ2(end)-LtwoFxQ2(end);
LonetwoFyQ3=zeros(1,length(ttQ3));
LonetwoFyQ3(1)=LoneFyQ2(end)+LtwoFyQ2(end);

%Q3区域质心，X坐标系下,质心相应变量设置
ComLxQ3=zeros(1,length(ttQ3));
ComLxQ3(1)=ComLxQ3(1);
ComVxQ3=zeros(1,length(ttQ3));
ComVxQ3(1)=ComVxQ3(1);
ComFxoneQ3=zeros(1,length(ttQ3)); %第一条腿给予质心的x方向力
ComFxtwoQ3=zeros(1,length(ttQ3)); %第二条腿给予质心的x方向力
ComFxQ3=zeros(1,length(ttQ3)); %给予质心的x方向合力
ComFxQ3(1)=ComFxQ2(end);
ComAxQ3=zeros(1,length(ttQ3));
ComAxQ3(1)=ComAxQ2(end);
%Q3区域质心，X坐标系下,质心相应变量设置
ComLyQ3=zeros(1,length(ttQ3));
ComVyQ3=zeros(1,length(ttQ3));
ComFyoneQ3=zeros(1,length(ttQ3));%第一条腿给予质心的y方向力
ComFytwoQ3=zeros(1,length(ttQ3)); %第二条腿给予质心的y方向力
ComFyQ3=zeros(1,length(ttQ3)); %给予质心的y方向合力
ComAyQ3=zeros(1,length(ttQ3));
%Q3区域XY坐标系下,质心相应变量设置
ComLxyQ3=zeros(1,length(ttQ3));
ComVxyQ3=zeros(1,length(ttQ3));
ComFxyoneQ3=zeros(1,length(ttQ3));
ComFxytwoQ3=zeros(1,length(ttQ3));
ComFxyQ3=zeros(1,length(ttQ3));
ComAxyQ3=zeros(1,length(ttQ3));
%%
FLQ3end=MSbp*g+MLbp*g-FM*cos(w*ttQ3(end)-forcefai);
QuaFunParaQ3=FittingCurve(2,ttQ30,LtwoFyQ3(1),ttQ3(end),-m*A*w^2*cos(w*ttQ3(end))+m*g+FLQ3end);
%%%还是按照第二条腿的进行计算

for i=1:length(ttQ3)  %从0的下一时刻开始
    
    if i>1
        %%
        %Q3区域Legone(除了力),x方向上的相应计算
        tspan = [ttQ30 ttQ3(i)];
        ComLVxQ30 = [LoneLxQ3(1) LoneVxQ3(1) ];%LoneLxQ30表示这一段时间初始横向位移，LoneLxQ(1)表示初始横向速度
        [t,y] = ode113(@(t,y) odefunddQ3(t,y,A,L,w,g,QuaFunParaQ3(1),QuaFunParaQ3(2),QuaFunParaQ3(3),T,d,m,MSbp,MLbp,FM,forcefai), tspan, ComLVxQ30);
        ComLxQ3(i)=y(end,1);
        ComVxQ3(i)=y(end,2);
        ComFxQ3(i)=(ComVxQ3(i)-ComVxQ3(i-1))/(ttQ3(i)-ttQ3(i-1))*(m+MSbp+MLbp); %(两腿合力)
        ComAxQ3(i)=(ComVxQ3(i)-ComVxQ3(i-1))/(ttQ3(i)-ttQ3(i-1)); %(两腿合力产生的加速度)
        %Q3区域Legone,x方向上的相应计算
        LoneLxQ3(i)=ComLxQ3(i);
        LoneVxQ3(i)=ComVxQ3(i);
        LoneAxQ3(i)=ComAxQ3(i);
        %Q3区域Legtwo,y方向上力
        LtwoFyQ3(i)=QuaFunParaQ3(1)*ttQ3(i)^2+QuaFunParaQ3(2)*ttQ3(i)+QuaFunParaQ3(3);
        
        %Q3区域Legone,y方向上的相应计算
        LoneLyQ3(i)=L-A+A*cos(w*ttQ3(i));
        LoneVyQ3(i)=-A*w*sin(w*ttQ3(i));
        FL=MSbp*g+MLbp*g-FM*cos(w*ttQ3(i)-forcefai);
        LoneFyQ3(i)=-m*A*w^2*cos(w*ttQ3(i))+m*g+FL-LtwoFyQ3(i);
        %LoneFyQ3(i)=-(2*(m*A*w^2+m*g+FL)*ttQ3(i))/(aerfa*T)+(1+aerfa)*(m*A*w^2+m*g+FL)/(2*aerfa);
        LoneAyQ3(i)=-A*w^2*cos(w*ttQ3(i));
        %Q3区域Legone（力），x方向上力的计算
        LoneFxQ3(i)=LoneFyQ3(i)*LoneLxQ3(i)/LoneLyQ3(i);
        %Q3区域Legtwo,x方向上的相应计算
        LtwoLxQ3(i)=d-LoneLxQ3(i);
        LtwoVxQ3(i)=-LoneVxQ3(i);
        LtwoFxQ3(i)=LoneFxQ3(i)-ComFxQ3(i); %合力减掉第一条腿的力
        LtwoAxQ3(i)=-LoneAxQ3(i);
        %Q3区域Legtwo,y方向上的相应计算
        LtwoLyQ3(i)= LoneLyQ3(i);
        LtwoVyQ3(i)=LoneVyQ3(i);
        %LtwoFyQ3(i)=-m*A*w^2*cos(w*ttQ3(i))+m*g+FL+(2*ttQ3(i)*(m*A*w^2+m*g+FL)/(aerfa*T))-(1+aerfa)*(m*A*w^2+m*g+FL)/(2*aerfa);
        LtwoAyQ3(i)=LoneAyQ3(i);    

    end
    %Q3区域Legone,XY方向上初始赋值及变量设置
    LoneLxyQ3(i)=sqrt(LoneLxQ3(i)^2+LoneLyQ3(i)^2);
    LoneVxyQ3(i)=sqrt(LoneVxQ3(i)^2+LoneVyQ3(i)^2);
    LoneFxyQ3(i)=sqrt(LoneFxQ3(i)^2+LoneFyQ3(i)^2);
    LoneAxyQ3(i)=sqrt(LoneAxQ3(i)^2+LoneAyQ3(i)^2);
    %Q3区域Legtwo,xY方向上初始赋值及变量设置
    LtwoLxyQ3(i)=sqrt(LtwoLxQ3(i)^2+LtwoLyQ3(i)^2);
    LtwoVxyQ3(i)=sqrt(LoneVxQ3(i)^2+LtwoVyQ3(i)^2);
    LtwoFxyQ3(i)=sqrt(LtwoFxQ3(i)^2+LtwoFyQ3(i)^2);
    LtwoAxyQ3(i)=sqrt(LtwoAxQ3(i)^2+LtwoAyQ3(i)^2);
    
    %Q3区域质心，X坐标系下,质心相应变量设置
    ComLxQ3(i)=LoneLxQ3(i);
    ComVxQ3(i)=LoneVxQ3(i);
    ComFxoneQ3(i)=LoneFxQ3(i); %第一条腿给予质心的x方向力
    ComFxtwoQ3(i)=LtwoFxQ3(i); %第二条腿给予质心的x方向力
    %ComFxQ3(i)=LoneFxQ3(i)-LtwoFxQ3(i); %给予质心的x方向合力
    %ComAxQ3(i)=LoneAxQ3(i);
    %Q3区域质心，X坐标系下,质心相应变量设置
    ComLyQ3(i)=LoneLyQ3(i);
    ComVyQ3(i)=LoneVyQ3(i);
    ComFyoneQ3(i)=LoneFyQ3(i);%第一条腿给予质心的y方向力
    ComFytwoQ3(i)=LtwoFyQ3(i); %第二条腿给予质心的y方向力
    ComFyQ3(i)=LoneFyQ3(i)+LtwoFyQ3(i); %给予质心的y方向合力
    ComAyQ3(i)=LoneAyQ3(i);
    %Q3区域质心，XY坐标系下,质心相应变量设置
    ComLxyQ3(i)=sqrt(ComLxQ3(i)^2+ComLyQ3(i)^2);
    ComVxyQ3(i)=sqrt(ComVxQ3(i)^2+ComVyQ3(i)^2);
    ComFxyoneQ3(i)=LoneFxyQ3(i);
    ComFxytwoQ3(i)=LtwoFxyQ3(i);
    ComFxyQ3(i)=sqrt(ComFxQ3(i)^2+ComFyQ3(i)^2);
    ComAxyQ3(i)=LoneAxyQ3(i); 
end


%%以步态周期为横坐标 
% figure(1)
% plot(ttQ3/T,LoneLxQ3,'-k',ttQ3/T,LoneLyQ3,'-.r',ttQ3/T,LoneLxyQ3,'--b',ttQ3/T,LtwoLxQ3,'-c',ttQ3/T,LtwoLyQ3,'-.m',ttQ3/T,LtwoLxyQ3,'--g','LineWidth',2)
% title('Lone-Com-Ltwo Displacement')
% xlabel('Gait cycle')
% ylabel('Displacement(m)');
% legend('Lone-ComH','Lone-ComV','Lone-Comxy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ3(1)/T,ttQ3(end)/T])
% ylim([0,1.7])
% figure(2)
% plot(ttQ3/T,LoneVxQ3,'-k',ttQ3/T,LoneVyQ3,'-.r',ttQ3/T,LoneVxyQ3,'--b',ttQ3/T,LtwoVxQ3,'-c',ttQ3/T,LtwoVyQ3,'-.m',ttQ3/T,LtwoVxyQ3,'--g','LineWidth',2)
% title('Lone-Com-two Velocity')
% xlabel('Gait cycle')
% ylabel('Velocity(m/s)');
% legend('Lone-ComH','Lone-ComV','Lone-Comxy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ3(1)/T,ttQ3(end)/T])
% figure(3)
% plot(ttQ3/T,LoneFxQ3,'-k',ttQ3/T,LoneFyQ3,'-.r',ttQ3/T,LoneFxyQ3,'--b',ttQ3/T,LtwoFxQ3,'-c',ttQ3/T,LtwoFyQ3,'-.m',ttQ3/T,LtwoFxyQ3,'--g','LineWidth',2)
% title('Lone-two Force')
% xlabel('Gait cycle')
% ylabel('Force(N)');
% legend('LoneH','LoneV','Lonexy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ3(1)/T,ttQ3(end)/T])
% figure(4)
% plot(ttQ3/T,LoneAxQ3,'-k',ttQ3/T,LoneAyQ3,'-.r',ttQ3/T,LoneAxyQ3,'--b',ttQ3/T,LtwoAxQ3,'-c',ttQ3/T,LtwoAyQ3,'-.m',ttQ3/T,LtwoAxyQ3,'--g','LineWidth',2)
% title('Lone-Com-two Acceleration')
% xlabel('Gait cycle')
% ylabel('Acceleration（m/s^2）')
% legend('Lone-ComH','Lone-ComV','Lone-Comxy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ3(1)/T,ttQ3(end)/T])

save('Q3','ttQ3','Subs','LoneLxQ3','LoneVxQ3','LoneFxQ3','LoneAxQ3','LoneLyQ3','LoneVyQ3','LoneFyQ3','LoneAyQ3','LoneLxyQ3','LoneVxyQ3','LoneFxyQ3','LoneAxyQ3',...
                        'LtwoLxQ3','LtwoVxQ3','LtwoFxQ3','LtwoAxQ3','LtwoLyQ3','LtwoVyQ3','LtwoFyQ3','LtwoAyQ3','LtwoLxyQ3','LtwoVxyQ3','LtwoFxyQ3','LtwoAxyQ3',...
                        'ComLxQ3','ComVxQ3','ComFxoneQ3','ComFxtwoQ3','ComFxQ3','ComAxQ3', 'ComLyQ3','ComVyQ3','ComFyoneQ3','ComFytwoQ3','ComFyQ3','ComAyQ3','ComLxyQ3','ComVxyQ3','ComFxyoneQ3','ComFxytwoQ3','ComFxyQ3','ComAxyQ3');