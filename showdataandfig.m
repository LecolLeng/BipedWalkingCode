%%%计算结果处理和显示

function [sumwork, Pushoffprowork,WorkSection1, WorkSection2,WorkSection3,WorkSection4]=showdataandfig()

load('par.mat');
load('Q1.mat');
load('Q2.mat');
load('Q3.mat');
load('Q4.mat');
L = Subs(1);
A = Subs(2);
T = Subs(3);
w = 2*2*pi/T;
m = Subs(4);
g = Subs(5);
d = Subs(6);
tdeta=Subs(9);

ttQ5=T/2+ttQ1;
ttQ6=T/2+ttQ2;
ttQ7=T/2+ttQ3;
ttQ8=T/2+ttQ4;
%tt=[ttQ1/T,ttQ2(2:end)/T,ttQ3(2:end)/T,ttQ4(2:end)/T,ttQ5(2:end)/T,ttQ6(2:end)/T,ttQ7(2:end)/T,ttQ8(2:end)/T];
tt=[ttQ1,ttQ2(2:end),ttQ3(2:end),ttQ4(2:end),ttQ5(2:end),ttQ6(2:end),ttQ7(2:end),ttQ8(2:end)];

%%两腿长度
LtwoLxQ1=LtwoLxQ2(1)*ones(1,length(ComLxQ1)); %%%%%%%%第1条腿在第一象限腿在x方向上长度为第二象限初始时刻的长度

LoneLxQ4=LoneLxQ3(end)*ones(1,length(ComLxQ4));    %%%%第1条腿在第四象限腿在x方向上长度为第三象限末端时刻的长度
LoneLxQ5=LtwoLxQ1;
LoneLxQ6=LtwoLxQ2;
LoneLxQ7=LtwoLxQ3;
LoneLxQ8=LtwoLxQ4;


LtwoLxQ5=LoneLxQ1;
LtwoLxQ6=LoneLxQ2;
LtwoLxQ7=LoneLxQ3;
LtwoLxQ8=LoneLxQ4;
LoneLx=[ LoneLxQ1  LoneLxQ2(2:end) LoneLxQ3(2:end) LoneLxQ4(2:end) LoneLxQ5(2:end) LoneLxQ6(2:end) LoneLxQ7(2:end) LoneLxQ8(2:end)];
LtwoLx=[ LtwoLxQ1  LtwoLxQ2(2:end) LtwoLxQ3(2:end) LtwoLxQ4(2:end) LtwoLxQ5(2:end) LtwoLxQ6(2:end) LtwoLxQ7(2:end) LtwoLxQ8(2:end)];


LtwoLyQ1=LtwoLyQ2(1)*ones(1,length(ttQ1));   %%%%%%%%第1条腿在第一象限腿在y方向上长度为第二象限初始时刻的长度
LoneLyQ4=LoneLyQ3(end)*ones(1,length(ttQ4));    %%%%第1条腿在第四象限腿在y方向上长度为第三象限末端时刻的长度
LoneLyQ5=LtwoLyQ1;
LoneLyQ6=LtwoLyQ2;
LoneLyQ7=LtwoLyQ3;
LoneLyQ8=LtwoLyQ4;

LtwoLyQ5=LoneLyQ1;
LtwoLyQ6=LoneLyQ2;
LtwoLyQ7=LoneLyQ3;
LtwoLyQ8=LoneLyQ4;
LoneLy=[ LoneLyQ1  LoneLyQ2(2:end)  LoneLyQ3(2:end)  LoneLyQ4(2:end)  LoneLyQ5(2:end)  LoneLyQ6(2:end)  LoneLyQ7(2:end)  LoneLyQ8(2:end) ];
LtwoLy=[ LtwoLyQ1, LtwoLyQ2(2:end)  LtwoLyQ3(2:end)  LtwoLyQ4(2:end)   LtwoLyQ5(2:end)  LtwoLyQ6(2:end)  LtwoLyQ7(2:end)  LtwoLyQ8(2:end)  ];

LtwoLxyQ1=LtwoLxyQ2(1)*ones(1,length(ttQ1));
LoneLxyQ4=LoneLxyQ3(end)*ones(1,length(ttQ4));
LoneLxyQ5=LtwoLxyQ1;
LoneLxyQ6=LtwoLxyQ2;
LoneLxyQ7=LtwoLxyQ3;
LoneLxyQ8=LtwoLxyQ4;


LtwoLxyQ5=LoneLxyQ1;
LtwoLxyQ6=LoneLxyQ2;
LtwoLxyQ7=LoneLxyQ3;
LtwoLxyQ8=LoneLxyQ4;
LoneLxy=[ LoneLxyQ1  LoneLxyQ2(2:end)  LoneLxyQ3(2:end)  LoneLxyQ4(2:end)  LoneLxyQ5(2:end)  LoneLxyQ6(2:end)  LoneLxyQ7(2:end)  LoneLxyQ8(2:end) ];
LtwoLxy=[ LtwoLxyQ1  LtwoLxyQ2(2:end)  LtwoLxyQ3(2:end)  LtwoLxyQ4(2:end)   LtwoLxyQ5(2:end)  LtwoLxyQ6(2:end)  LtwoLxyQ7(2:end)  LtwoLxyQ8(2:end)  ];
figure(1)
plot(tt,LoneLxy,'-r',tt,LtwoLxy,'-b','LineWidth',2);
title('Lone-Ltwo Length')
%xlabel('Gait cycle')
xlabel('time(s)')
ylabel('Length(m)')
legend('LoneLength','LtwoLength');
 xlim([0,tt(end)]);
%%两腿力大小
LoneFxQ4=0*ones(1,length(ttQ4));
LoneFxQ5=0*ones(1,length(ttQ5));
LoneFxQ6=LtwoFxQ2;
LoneFxQ7=LtwoFxQ3;
LoneFxQ8=LtwoFxQ4;

LoneFyQ4=0*ones(1,length(ttQ4));
LoneFyQ5=0*ones(1,length(ttQ5));
LoneFyQ6=LtwoFyQ2;
LoneFyQ7=LtwoFyQ3;
LoneFyQ8=LtwoFyQ4;
LoneFxyQ4=0*ones(1,length(ttQ4));
LoneFxyQ5=0*ones(1,length(ttQ5));
LoneFxyQ6=LtwoFxyQ2;
LoneFxyQ7=LtwoFxyQ3;
LoneFxyQ8=LtwoFxyQ4;

LtwoFxQ1=0*ones(1,length(ttQ1));
LtwoFxQ5=LoneFxQ1;
LtwoFxQ6=LoneFxQ2;
LtwoFxQ7=LoneFxQ3;
LtwoFxQ8=LoneFxQ4;
LtwoFyQ1=0*ones(1,length(ttQ1));
LtwoFyQ5=LoneFyQ1;
LtwoFyQ6=LoneFyQ2;
LtwoFyQ7=LoneFyQ3;
LtwoFyQ8=LoneFyQ4;
LtwoFxyQ1=0*ones(1,length(ttQ1));
LtwoFxyQ5=LoneFxyQ1;
LtwoFxyQ6=LoneFxyQ2;
LtwoFxyQ7=LoneFxyQ3;
LtwoFxyQ8=LoneFxyQ4;
LoneFx=[ LoneFxQ1 LoneFxQ2(2:end)  LoneFxQ3(2:end)  LoneFxQ4(2:end)  LoneFxQ5(2:end)  LoneFxQ6(2:end)  LoneFxQ7(2:end)  LoneFxQ8(2:end) ];
LtwoFx=[ LtwoFxQ1 LtwoFxQ2(2:end)  LtwoFxQ3(2:end)  LtwoFxQ4(2:end)  LtwoFxQ5(2:end)  LtwoFxQ6(2:end)  LtwoFxQ7(2:end)  LtwoFxQ8(2:end) ];

LoneFy=[ LoneFyQ1 LoneFyQ2(2:end)  LoneFyQ3(2:end)  LoneFyQ4(2:end)  LoneFyQ5(2:end)  LoneFyQ6(2:end)  LoneFyQ7(2:end)  LoneFyQ8(2:end) ];
LtwoFy=[ LtwoFyQ1 LtwoFyQ2(2:end)  LtwoFyQ3(2:end)  LtwoFyQ4(2:end)  LtwoFyQ5(2:end)  LtwoFyQ6(2:end)  LtwoFyQ7(2:end)  LtwoFyQ8(2:end) ];

LoneFxy=[ LoneFxyQ1 LoneFxyQ2(2:end)  LoneFxyQ3(2:end)  LoneFxyQ4(2:end)  LoneFxyQ5(2:end)  LoneFxyQ6(2:end)  LoneFxyQ7(2:end)  LoneFxyQ8(2:end) ];
LtwoFxy=[ LtwoFxyQ1 LtwoFxyQ2(2:end)  LtwoFxyQ3(2:end)  LtwoFxyQ4(2:end)  LtwoFxyQ5(2:end)  LtwoFxyQ6(2:end)  LtwoFxyQ7(2:end)  LtwoFxyQ8(2:end) ];
figure(2)
plot(tt,LoneFx,'-.r',tt,LoneFy,'--r',tt,LoneFxy,'-r',tt,LtwoFx,'-.b',tt,LtwoFy,'--b',tt,LtwoFxy,'-b','LineWidth',2);

title('Lone-Ltwo Force')
%xlabel('Gait cycle')
xlabel('time(s)')
ylabel('Force(N)')
legend('LoneFx','LoneFy','LoneFxy','LtwoFx','LtwoFy','LtwoFxy');
 xlim([0,tt(end)]);
%%质心处的力
ComFxQ5=ComFxQ1;
ComFxQ6=ComFxQ2;
ComFxQ7=ComFxQ3;
ComFxQ8=ComFxQ4;
ComFyQ5=ComFyQ1;
ComFyQ6=ComFyQ2;
ComFyQ7=ComFyQ3;
ComFyQ8=ComFyQ4;
ComFxyQ5=ComFxyQ1;
ComFxyQ6=ComFyQ2;
ComFxyQ7=ComFyQ3;
ComFxyQ8=ComFyQ4;
ComFx=[ ComFxQ1  ComFxQ2(2:end)  ComFxQ3(2:end)  ComFxQ4(2:end)  ComFxQ5(2:end)  ComFxQ6(2:end) ComFxQ7(2:end)  ComFxQ8(2:end) ];
ComFy=[ ComFyQ1  ComFyQ2(2:end)  ComFyQ3(2:end)  ComFyQ4(2:end)  ComFyQ5(2:end)  ComFyQ6(2:end)  ComFyQ7(2:end)  ComFyQ8(2:end) ];
ComFxy=[ ComFxyQ1 ComFxyQ2(2:end)  ComFxyQ3(2:end)  ComFxyQ4(2:end)  ComFxyQ5(2:end)  ComFxyQ6(2:end)  ComFxyQ7(2:end)  ComFxyQ8(2:end) ];
figure(3)
plot(tt,ComFx,tt,ComFy,tt,ComFxy,'LineWidth',2);
title('Com Force')
%xlabel('Gait cycle')
xlabel('time(s)')
ylabel('Force(N)')
legend('ComX','ComY','ComXY');
 xlim([0,tt(end)]);
%质心位移
ComLxQ5=d+ComLxQ1;
ComLxQ6=d+ComLxQ2;
ComLxQ7=d+ComLxQ3;
ComLxQ8=d+ComLxQ4;
ComLyQ5=ComLyQ1;
ComLyQ6=ComLyQ2;
ComLyQ7=ComLyQ3;
ComLyQ8=ComLyQ4;
ComLxyQ5=(ComLxQ5.^2+ComLyQ5.^2).^0.5;
ComLxyQ6=(ComLxQ6.^2+ComLyQ6.^2).^0.5;
ComLxyQ7=(ComLxQ7.^2+ComLyQ7.^2).^0.5;
ComLxyQ8=(ComLxQ8.^2+ComLyQ8.^2).^0.5;

ComLx=[ ComLxQ1, ComLxQ2(2:end)  ComLxQ3(2:end)  ComLxQ4(2:end)  ComLxQ5(2:end)  ComLxQ6(2:end)  ComLxQ7(2:end)  ComLxQ8(2:end)  ];
ComLy=[ ComLyQ1, ComLyQ2(2:end)  ComLyQ3(2:end)  ComLyQ4(2:end)  ComLyQ5(2:end)  ComLyQ6(2:end)  ComLyQ7(2:end)  ComLyQ8(2:end) ];
ComLxy=[ ComLxyQ1, ComLxyQ2(2:end)  ComLxyQ3(2:end)  ComLxyQ4(2:end)  ComLxyQ5(2:end)  ComLxyQ6(2:end)  ComLxyQ7(2:end)  ComLxyQ8(2:end) ];
figure(4)
plot(tt,ComLx,tt,ComLy,tt,ComLxy,'LineWidth',2);
title('Com Displacement')
%xlabel('Gait cycle')
xlabel('time(s)')
ylabel('Displacement(m)')
legend('ComX','ComY','ComXY');
 xlim([0,tt(end)]);
%%计算肌肉做的功
LonePx=zeros(1,length(tt));
LonePy=zeros(1,length(tt));
LonePxy=zeros(1,length(tt));
LtwoPx=zeros(1,length(tt));
LtwoPy=zeros(1,length(tt));
LtwoPxy=zeros(1,length(tt));
SumP=zeros(1,length(tt));

LoneWx=zeros(1,length(tt));
LoneWy=zeros(1,length(tt));
LoneWxy=zeros(1,length(tt));
LtwoWx=zeros(1,length(tt));
LtwoWy=zeros(1,length(tt));
LtwoWxy=zeros(1,length(tt));
SumW=zeros(1,length(tt));
for i=1:length(tt);
    if i>1
        
        %         LonePx(i)=((LoneLx(i)-LoneLx(i-1))*LoneFx(i))/tdeta; %x方向上功率的变化
        %         LonePy(i)=((LoneLy(i)-LoneLy(i-1))*LoneFy(i))/tdeta;%y方向上功率的变化
        LonePxy(i)= ((LoneLxy(i)-LoneLxy(i-1))*LoneFxy(i))/tdeta;
        %         LtwoPx(i)=((LtwoLx(i)-LtwoLx(i-1))*LtwoFx(i))/tdeta; %x方向上功率的变化
        %         LtwoPy(i)=((LtwoLy(i)-LtwoLy(i-1))*LtwoFy(i))/tdeta;%y方向上功率的变化
        LtwoPxy(i)=((LtwoLxy(i)-LtwoLxy(i-1))*LtwoFxy(i))/tdeta;
        SumP(i)=LonePxy(i)+LtwoPxy(i);
    end
    
    %     LoneWx(i)=sum(abs(LonePx))*tdeta;
    %     LoneWy(i)=sum(abs(LonePy))*tdeta;
    LoneWxy(i)=sum(abs(LonePxy))*tdeta;
    %     LtwoWx(i)=sum(abs(LtwoPx))*tdeta;
    %     LtwoWy(i)=sum(abs(LtwoPy))*tdeta;
    LtwoWxy(i)=sum(abs(LtwoPxy))*tdeta;
    SumW(i)=LoneWxy(i)+LtwoWxy(i);
end
figure(5)
% plot(tt,LonePx,'-.r',tt,LonePy,'--r',tt,LonePxy,'-r',tt,LtwoPx,'-.b',tt,LtwoPy,'--b',tt,LtwoPxy,'-b',tt,SumP,'-k','LineWidth',2);
% title('Power')
% xlabel('Gait cycle')
% ylabel('Power(W)')
% legend('LonePx','LonePy','LonePxy','LtwoPx','LtwoPy','LtwoPxy','SumP');
plot(tt,LonePxy,'-r',tt,LtwoPxy,'-b','LineWidth',2); %%,tt,SumP,'-k'
title('Power')
%xlabel('Gait cycle')
xlabel('time(s)')
ylabel('Power(W)')
legend('LonePxy','LtwoPxy');  %%,'SumP'
 xlim([0,tt(end)]);
figure(6)
% plot(tt,LoneWx,'-.r',tt,LoneWy,'--r',tt,LoneWxy,'-r',tt,LtwoWx,'-.b',tt,LtwoWy,'--b',tt,LtwoWxy,'-b',tt,SumW,'-k','LineWidth',2);
% title('Work')
% xlabel('Gait cycle')
% ylabel('Wower(J)')
% legend('LoneWx','LoneWy','LoneWxy','LtwoWx','LtwoWy','LtwoWxy','SumW');
plot(tt,LoneWxy,'-r',tt,LtwoWxy,'-b',tt,SumW,'-k','LineWidth',2);
title('Work')
%xlabel('Gait cycle')
xlabel('time(s)')
ylabel('Wower(J)')
legend('LoneWxy','LtwoWxy','SumW');
 xlim([0,tt(end)]);
sumwork=SumW(end);

num=0;
Pushoffprowork=0;
for i=1:length(LtwoPxy)
    
    if LtwoPxy(i)<0 && LtwoPxy(i+1)>0
        num=num+1;
    end
    if num==2
        Pushoffprowork=Pushoffprowork+LtwoPxy(i);
    end 
end
Pushoffprowork=Pushoffprowork*tdeta;

WorkSection1=0;
WorkSection2=0;
for i=1:length(LtwoPxy)/2
    if LtwoPxy(i)<0
        WorkSection1=WorkSection1+LtwoPxy(i)*tdeta;
    end
    if LtwoPxy(i)>0
        WorkSection2=WorkSection2+LtwoPxy(i)*tdeta;
    end
end
WorkSection3=0;
WorkSection4=0;
for i=1:length(LonePxy)/2
    if LonePxy(i)<0
        WorkSection3=WorkSection1+LonePxy(i)*tdeta;
    end
    if LonePxy(i)>0
        WorkSection4=WorkSection2+LonePxy(i)*tdeta;
    end
end












