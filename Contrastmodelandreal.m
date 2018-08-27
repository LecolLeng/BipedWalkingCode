%%%%模拟数据与实际测试数据对比
clc
clear
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
aerfa=Subs(8);
tdeta=Subs(9);
ttQ5=T/2+ttQ1;
ttQ6=T/2+ttQ2;
ttQ7=T/2+ttQ3;
ttQ8=T/2+ttQ4;
tt=[ttQ1/T,ttQ2(2:end)/T,ttQ3(2:end)/T,ttQ4(2:end)/T,ttQ5(2:end)/T,ttQ6(2:end)/T,ttQ7(2:end)/T,ttQ8(2:end)/T];
%tt=[ttQ1,ttQ2(2:end),ttQ3(2:end),ttQ4(2:end),ttQ5(2:end),ttQ6(2:end),ttQ7(2:end),ttQ8(2:end)];


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
% figure(2)
% plot(tt,LoneFx,'-.r',tt,LoneFy,'--r',tt,LoneFxy,'-r',tt,LtwoFx,'-.b',tt,LtwoFy,'--b',tt,LtwoFxy,'-b','LineWidth',2);
% 
% title('Lone-Ltwo Force')
% xlabel('Gait cycle')
% ylabel('Force(N)')
% legend('LoneFx','LoneFy','LoneFxy','LtwoFx','LtwoFy','LtwoFxy');

LtwoFxGaitCir=[LtwoFx(length(ttQ1):end),LtwoFx(1:length(ttQ1)-1)];
LtwoFxGaitCir=[-LtwoFxGaitCir(1:length(ttQ2)+length(ttQ3)+length(ttQ4)),LtwoFxGaitCir(length(ttQ2)+length(ttQ3)+length(ttQ4)+1:end)];
LtwoFyGaitCir=[LtwoFy(length(ttQ1):end),LtwoFy(1:length(ttQ1)-1)];
LtwoFxyGaitCir=[LtwoFxy(length(ttQ1):end),LtwoFxy(1:length(ttQ1)-1)];
real=xlsread('Opensimdata.xlsx');
figure(1)
tt=100*tt;
plot(tt,LtwoFxGaitCir,'-.b',tt,LtwoFyGaitCir,'-b',real(:,2)*100,real(:,3),'-.r',real(:,2)*100,real(:,4),'-r','LineWidth',2);
legend('ModelFx','ModelFy','RealFx','RealFy');
xlabel('Gait cycle(%)')
ylabel('Force(N)')
% x=1:0.1:20;
% z=lagrange(x0,y0,x);
% plot(x,z,':',x0,y0,'ko');

