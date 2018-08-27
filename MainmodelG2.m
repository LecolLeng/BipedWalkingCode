%% Editor:Yuquan Leng
%%%第二版，模型推演，已知人体的身高、质量、平均步速等
clc
clear
Gra_acc=9.8; %重力加速度
g=Gra_acc;
Cyclenum=0;
%% 设定人体的固有参数，身高、体重、平均步速
Human_heigh=1.77; %单位m
Human_weight=70.8;%单位Kg
Human_Velocity=1.7;%步行平均速度，单位m/s
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
    Foot_length=(Human_heigh*100-77.059)/3.536/100; %单位为m
else
    Foot_length=(Human_heigh*100-68.586)/4.0636/100; %单位为m
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
a=0.07;
b=(0.5-a)*(Human_heigh^0.5)/1.87;
tswing=a*Gait_cycle_time+b;
Doublesup_time=Gait_cycle_time-2*tswing;
Double_leg_support_time_rate=Doublesup_time/Gait_cycle_time;

%*******人体质心振幅，单位m*****************参考文献：Hoover J, Meguid S A. Performance assessment of the suspended-load backpack[J]. International Journal of Mechanics and Materials in Design, 2011, 7(2): 111-121.
v=Human_Velocity;
Com_amplitude=(Leg_real_length/2)*(1-(1-(0.936*v/(Leg_real_length*2*1.504*(v/Leg_real_length)^0.57))^2)^0.5);

%% 计算出基于虚拟旋转点模型的腿长（与步长、脚长、双腿支撑期占比、人实际腿长）
%具体公式推导，参考文档(虚拟连杆方式)
Leg_virtual_length=(-Com_amplitude*Foot_length-Foot_length*Com_amplitude*...
    cos(Double_leg_support_time_rate*pi)+Leg_real_length*Gait_step_half_length*...
    (1+Double_leg_support_time_rate))/(Gait_step_half_length*(1+Double_leg_support_time_rate)-Foot_length);
%非虚拟连杆方式
%Leg_virtual_length=0.53*Human_heigh;
%% 初始速度（初始速度小于平均速度），Q2双腿支撑期所占总双腿支撑期的比例，两个作为自变量
%目标函数为做功的总和（纵向是相等的，只需要约束横向的即可）及计算的步长与实际步长的插值，这两者的绝对值应该等于0（理想情况）
%设定目标函数值小于一定的值deita，即可认为需找到了解
%该问题是一个三维曲面求最小点的问题，曲面形式未知，可能存在极值点，因此，需要采用一定的手段解决
%暂时预使用的方法为：设置多个初始点，依次进行寻找，寻找到了解，则跳出循环。如果没有则换一个初始点（只要初始点足够多，则必定能找到解）
%%
A=Com_amplitude;
w=2*2*pi*Gait_frequence; %质心振动频率
T=1/Gait_frequence;%步态周期

%设定自变量范围
Max_Ini_speed=0.95*Human_Velocity;  %初始速度肯定是小于平均速度
Min_Ini_speed=0.85*Human_Velocity;
MaxQ2_Double_support_percent=0.51*Double_leg_support_time_rate;
MinQ2_Double_support_percent=0.49*Double_leg_support_time_rate;
%scale=(Max_Ini_speed-Min_Ini_speed)/(MaxQ2_Double_support_percent-MinQ2_Double_support_percent); %进行范围的比例计算，便于不同方向上步长的对等
scale=4;
%%%%弹性背包的刚度及阻尼的设定
KK=[4000:500:4000 ];
CC=[141:100:141];
sumwork=zeros(length(KK),length(CC));
MaxBrel=zeros(length(KK),length(CC));
MaxFL=zeros(length(KK),length(CC));
MaxHeelStr=zeros(length(KK),length(CC));
MaxHeelPush=zeros(length(KK),length(CC));
Changfai=zeros(length(KK),length(CC));
Pushoffprowork=zeros(length(KK),length(CC));

WorkSection1=zeros(length(KK),length(CC));
WorkSection2=zeros(length(KK),length(CC));
WorkSection3=zeros(length(KK),length(CC));
WorkSection4=zeros(length(KK),length(CC));
ObjectNewAll=zeros(length(KK),length(CC));

for j=1:length(KK)
    for k=1:length(CC)
        %背包相关的参数1）背包固定部分质量2）背包弹性负载质量3）刚度4）阻尼
        Subsbackpack=[0.01,0.01,KK(j),CC(k)];
        MSbp=Subsbackpack(1);
        MLbp=Subsbackpack(2);
        K=Subsbackpack(3);
        C=Subsbackpack(4);
        %力的幅值FM，运动相对的相位差fai，作用在肩部的力相对于运动的相位差forcefai，负载绝对振动的振幅dBabs
        [FM,fai,forcefai,Babs]=elasticbackpack(Subsbackpack,A,w,g);
        %%%%%有待进一步优化%%%%
        ttt=0:0.05:T;
        %%求一下包和负载的相对位移的幅值，作为约束条件
        Brel=Babs*cos(w*ttt-fai)-A*cos(w*ttt);
        %  figure(1)
        %  plot(ttt,Brel)
        MaxBrel(j,k)=max(Brel);
        %figure(2)
        FL=MSbp*g+MLbp*g-FM*cos(w*ttt-forcefai);
        %plot(ttt,FLKKCC)
        MaxFL(j,k)=max(FL);
        if MaxBrel(j,k)<0.2   %如果相对振幅小于0.2m
            %在定义域区域内设置nXn的初始点,内个初始点位于其中的一块内（整个定义域切成nxn块）
            n=3;
            Time_gap=0.001; %每次积分的时间间隔
            Objdeta=30;
            for ni=1:n
                for nj=1:n
                    
                    %在某一个区域块内，随机选取初始速度和Q2支撑期占双腿支撑百分比的初始值
                    Smallarea_speed=(Max_Ini_speed-Min_Ini_speed)/n;
                    Gone_Def_Ini_speed= unifrnd(Min_Ini_speed+Smallarea_speed*(ni-1),Min_Ini_speed+Smallarea_speed*ni);
                    Smallarea_percent=(MaxQ2_Double_support_percent-MinQ2_Double_support_percent)/n;
                    Gone_Def_Q2percent=unifrnd(MinQ2_Double_support_percent+Smallarea_percent*(ni-1),MinQ2_Double_support_percent+Smallarea_percent*ni);
                    
                    %定义人体相关参数,分别为1)腿长m、2)质心振幅m、3)步态周期s、4)人体质量kg、
                    %%5)重力加速度、6)半步长(参考值)、7)最高点质心横向速度、8)双腿支撑起比例 9)时间间隔
                    %%10）Q2双腿支撑起占整个的比例
                    Subs = [Leg_virtual_length,Com_amplitude,Gait_cycle_time,Human_weight,Gra_acc,Gait_step_half_length,Gone_Def_Ini_speed,Double_leg_support_time_rate,Time_gap,Gone_Def_Q2percent];
                    L = Subs(1);
                    A = Subs(2);
                    T = Subs(3);
                    w = 2*2*pi/T;
                    m = Subs(4);
                    g = Subs(5);
                    d = Subs(6);
                    aerfa=Subs(8);
                    tdeta=Subs(9);
                    
                    flage=1;
                    Evo_step_length=0.04; %设定初始的演化步长
                    Numdirection=5;
                    ObjectNew=zeros(1,Numdirection);
                    while (flage) %这个循环中需要找到极值点，假如极值点出的值为0，则认为是最值点，跳出所有循环，输出结果
                        
                        %根据一个点，朝着周边4个方向扩展，直接计算5个点的
%                         Speed_Derive_points=[Gone_Def_Ini_speed+Evo_step_length,Gone_Def_Ini_speed-Evo_step_length,Gone_Def_Ini_speed,Gone_Def_Ini_speed,...
%                             Gone_Def_Ini_speed+Evo_step_length*sqrt(2)/2,Gone_Def_Ini_speed+Evo_step_length*sqrt(2)/2,Gone_Def_Ini_speed-Evo_step_length*sqrt(2)/2,Gone_Def_Ini_speed-Evo_step_length*sqrt(2)/2,Gone_Def_Ini_speed];
                        Speed_Derive_points=[Gone_Def_Ini_speed+Evo_step_length,Gone_Def_Ini_speed-Evo_step_length,Gone_Def_Ini_speed,Gone_Def_Ini_speed,Gone_Def_Ini_speed];
                        Q2percent_Evo_step_length= Evo_step_length/scale;
%                         Q2percent_Derive_points=[Gone_Def_Q2percent,Gone_Def_Q2percent,Gone_Def_Q2percent+Q2percent_Evo_step_length,Gone_Def_Q2percent-Q2percent_Evo_step_length,...
%                             Gone_Def_Q2percent+Q2percent_Evo_step_length*sqrt(2)/2,Gone_Def_Q2percent-Q2percent_Evo_step_length*sqrt(2)/2,Gone_Def_Q2percent+Q2percent_Evo_step_length*sqrt(2)/2,Gone_Def_Q2percent-Q2percent_Evo_step_length*sqrt(2)/2,Gone_Def_Q2percent];
                        Q2percent_Derive_points=[Gone_Def_Q2percent,Gone_Def_Q2percent,Gone_Def_Q2percent+Q2percent_Evo_step_length,Gone_Def_Q2percent-Q2percent_Evo_step_length,Gone_Def_Q2percent];
                        for Directioni=1: Numdirection
                            Subs(7)=Speed_Derive_points(Directioni);
                            Subs(10)=Q2percent_Derive_points(Directioni);
                            save('Par','Subs','FM','fai','forcefai','Babs','MSbp','MLbp');
                            Lyq_CacuCosModQ1();
                            Lyq_CacuCosModQ2();
                            Lyq_CacuCosModQ3();
                            Lyq_CacuCosModQ4();
                            load('Q1.mat');
                            load('Q2.mat');
                            load('Q3.mat');
                            load('Q4.mat');
                            Hor_Psum=ComVxQ1*ComFxQ1'+ComVxQ2*ComFxQ2'+ComVxQ3*ComFxQ3'+ComVxQ4*ComFxQ4';
                            ObjectNew(Directioni)= abs(Gait_step_half_length-ComLxQ4(end))*200+abs(Hor_Psum*Time_gap);%乘以10的目的，是为了调整影响尺度。 %%%横向做功的总和
                            Cyclenum=Cyclenum+1;
                        end
                        % 5个目标数据，从中选取梯度最大的，假如第一个点位于极值点，则所有梯度为正
                        %假如已经是极值点，应难满足边界条件，满足目标函数值接近0条件，
                        %已经是5个点中的极值点，但是不满足目标接近0条件，调整步长，进一步优化。
                        %进一步优化后依然不行，则说明陷入了非0的极值点，选取其他初始点，进行搜索
                        
                        %求取ObjectNew中第5个数值与前4个的差值，差值越大越好，说明梯度大
                        Gradient_obj=ObjectNew(Numdirection)*ones(1,length(ObjectNew)-1)-ObjectNew(1,1: (Numdirection-1));
                        [Max_Gradient,Max_Gradient_num]=max(Gradient_obj); %得到最大梯度值，及编号
                        %如果最大梯度值是负数，说明改点已经是这5个点中的极值点
                        if Max_Gradient<=0
                            %梯度之都小于0，进行相应的判断
                            switch Evo_step_length
                                case 0.04
                                    Evo_step_length=0.02;
                                case 0.02
                                    Evo_step_length=0.01;
                                case 0.01
                                    Evo_step_length=0.005;
                                 case 0.005
                                    Evo_step_length=0.002;
                                  case 0.002
                                      Evo_step_length=0.001;
                                case 0.001
                                    if  ObjectNew(5)<Objdeta%%%%%%%%找到解了，显示并跳出
                                        flage=0;
                                        ni=4;
                                        nj=4;
                                        %找到了合适的参数，来输出相应的计算值
                                        [sumwork(j,k), Pushoffprowork(j,k),WorkSection1(j,k), WorkSection2(j,k),WorkSection3(j,k),WorkSection4(j,k)] =showdataandfig();
                                        MaxHeelStr(j,k)=LtwoFxyQ3(end);
                                        MaxHeelPush(j,k)=LoneFxyQ1(end);
                                        Changfai(j,k)=fai; 
                                        ObjectNewAll(j,k)=ObjectNew(5); %%存储极点值
                                        break;%%%%%%%%找到解了，
                                    else
                                        %输出这一区块找到了极小值，没有找到最小值
                                        sprintf('%s，%d 区块找到了极小值，没有找到最小值，将继续寻找',ni,nj);
                                    end
                            end
                        else
                            %有梯度值大于0的，将最大梯度点的数值作为Newpoint，最速梯度下降法
                            Gone_Def_Ini_speed=Speed_Derive_points(Max_Gradient_num);
                            Gone_Def_Q2percent=Q2percent_Derive_points(Max_Gradient_num);
                        end
                    end
                end
            end
        end
    end
end


figure(8)
for k=1:length(CC)
    hold on
    plot(KK,sumwork(:,k),'LineWidth',2)
end
title('sumwork')
xlabel('Stiffness(N/m)')
ylabel('Work(J)')
legend('c=50','c=100','c=150','c=200','c=250');


figure(9)
for k=1:length(CC)
    hold on
    plot(KK,MaxBrel(:,k),'LineWidth',2)
end
title('MaxBrel')
xlabel('Stiffness(N/m)')
ylabel('Displacement(m)')
legend('c=50','c=100','c=150','c=200','c=250');


figure(10)
for k=1:length(CC)
    hold on
    plot(KK,MaxFL(:,k),'LineWidth',2)
end
title('MaxFL')
xlabel('Stiffness(N/m)')
ylabel('Force(N)')
legend('c=50','c=100','c=150','c=200','c=250');

figure(11)
for k=1:length(CC)
    hold on
    plot(KK,MaxHeelStr(:,k),'LineWidth',2)
end
title('MaxHeelStr')
xlabel('Stiffness(N/m)')
ylabel('Force(N)')
legend('c=50','c=100','c=150','c=200','c=250');

figure(12)
for k=1:length(CC)
    hold on
    plot(KK,MaxHeelPush(:,k),'LineWidth',2)
end
title('MaxHeelPush')
xlabel('Stiffness(N/m)')
ylabel('Force(N)')
legend('c=50','c=100','c=150','c=200','c=250');

figure(13)
for k=1:length(CC)
    hold on
    plot(KK,Changfai(:,k),'LineWidth',2)
end
title('Changfai')
xlabel('Stiffness(N/m)')
ylabel('Degree(°)')
legend('c=50','c=100','c=150','c=200','c=250');

figure(14)
for k=1:length(CC)
    hold on
    plot(KK,Pushoffprowork(:,k),'LineWidth',2)
    
end
title('Pushoffprowork')
xlabel('Stiffness(N/m)')
ylabel('Work(J)')
legend('c=50','c=100','c=150','c=200','c=250');


figure(15)
for k=1:length(CC)
    hold on
    plot(KK,WorkSection1(:,k),'LineWidth',2)
    
end
title('WorkSection1')
xlabel('Stiffness(N/m)')
ylabel('Work(J)')
legend('c=50','c=100','c=150','c=200','c=250');


figure(16)
for k=1:length(CC)
    hold on
    plot(KK,WorkSection2(:,k),'LineWidth',2)
    
end
title('WorkSection2')
xlabel('Stiffness(N/m)')
ylabel('Work(J)')
legend('c=50','c=100','c=150','c=200','c=250');

figure(16)
for k=1:length(CC)
    hold on
    plot(KK,WorkSection3(:,k),'LineWidth',2) 
end
title('WorkSection3')
xlabel('Stiffness(N/m)')
ylabel('Work(J)')
legend('c=50','c=100','c=150','c=200','c=250');

figure(17)
for k=1:length(CC)
    hold on
    plot(KK,WorkSection4(:,k),'LineWidth',2)
    
end
title('WorkSection4')
xlabel('Stiffness(N/m)')
ylabel('Work(J)')
legend('c=50','c=100','c=150','c=200','c=250');

a=0;