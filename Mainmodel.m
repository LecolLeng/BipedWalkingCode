%%%
clc
clear
Step=0.650:0.001:0.650;
Stepdeta=zeros(1,length(Step));
Axend=zeros(1,length(Step));
for i=1:length(Step)    %%为了匹配合适的步长，使设定的步长和积分出的步长基本相等
    %定义人体相关参数,分别为1)腿长m、2)质心振幅m、3)步态周期s、4)人体质量kg、
    %%5)重力加速度、6)半步长(参考值)、7)最高点质心横向速度、8)双腿支撑起比例 9)时间间隔
    Subs = [1.5,0.02,1.25,70,9.8,Step(i),0.95,0.34,0.001];
    L = Subs(1);
    A = Subs(2);
    T = Subs(3);
    w = 2*2*pi/T;
    m = Subs(4);
    g = Subs(5);
    d = Subs(6);
    aerfa=Subs(8);
    tdeta=Subs(9);
    KK=[30000:400:30000 ];
    CC=[100:40:100];
    
    sumwork=zeros(length(KK),length(CC));
    MaxBrel=zeros(length(KK),length(CC));
    MaxFL=zeros(length(KK),length(CC));
    MaxHeelStr=zeros(length(KK),length(CC));
    MaxHeelPush=zeros(length(KK),length(CC));
    Changfai=zeros(length(KK),length(CC));
    Pushoffprowork=zeros(length(KK),length(CC));
    for j=1:length(KK)
        for k=1:length(CC)
            %背包相关的参数1）背包固定部分质量2）背包弹性负载质量3）刚度4）阻尼
            Subsbackpack=[0.02,20,KK(j),CC(k)];
            MSbp=Subsbackpack(1);
            MLbp=Subsbackpack(2);
            K=Subsbackpack(3);
            C=Subsbackpack(4);
            %需要得到什么??   %力的幅值FM，运动相对的相位差fai，作用在肩部的力相对于运动的相位差forcefai，负载绝对振动的振幅dBabs
            [FM,fai,forcefai,Babs]=elasticbackpack(Subsbackpack,A,w,g);
            ttt=0:0.05:T;
            %     %求一下包和负载的相对位移的幅值，作为约束条件
            Brel=Babs*cos(w*ttt-fai)-A*cos(w*ttt);
%             figure(1)
%             plot(ttt,Brel)
            MaxBrel(j,k)=max(Brel);
            %figure(2)
            FL=MSbp*g+MLbp*g-FM*cos(w*ttt-forcefai);
            %plot(ttt,FLKKCC)
            MaxFL(j,k)=max(FL);
            if MaxBrel(j,k)<0.1
                save('Par','Subs','FM','fai','forcefai','Babs','MSbp','MLbp');
                Lyq_CacuCosModQ1();
                Lyq_CacuCosModQ2();
                Lyq_CacuCosModQ3();
                Lyq_CacuCosModQ4();
                load('Q1.mat');
                load('Q2.mat');
                load('Q3.mat');
                load('Q4.mat');
                Stepdeta(i)=abs(Step(i)-ComLxQ4(end));
                Axend(i)=ComAxQ4(end);
                [sumwork(j,k), Pushoffprowork(j,k)] =showdataandfig();
                MaxHeelStr(j,k)=LtwoFxyQ3(end);
                MaxHeelPush(j,k)=LoneFxyQ1(end);
                Changfai(j,k)=fai;
                
            end
        end
    end
    
end


%用于检测步长是否合理
% figure(6)
% plot(Step,Stepdeta);
% figure(7)
% plot(Step,Axend);
figure(8)
for k=1:length(CC)
    hold on
    plot(KK,sumwork(:,k))
end
title('sumwork')
figure(9)
for k=1:length(CC)
    hold on
    plot(KK,MaxBrel(:,k))
end
title('MaxBrel')
figure(10)
for k=1:length(CC)
    hold on
    plot(KK,MaxFL(:,k))
end
title('MaxFL')
figure(11)
for k=1:length(CC)
    hold on
    plot(KK,MaxHeelStr(:,k))
end
title('MaxHeelStr')
figure(12)
for k=1:length(CC)
    hold on
    plot(KK,MaxHeelPush(:,k))
end
title('MaxHeelPush')
figure(13)
for k=1:length(CC)
    hold on
    plot(KK,Changfai(:,k))
end
title('Changfai')
figure(14)
for k=1:length(CC)
    hold on
    plot(KK,Pushoffprowork(:,k),'LineWidth',2)
    
end
title('Pushoffprowork')
xlabel('刚度(N/m)')
ylabel('Wower(J)')
legend('c=80','c=120','c=160','c=200','c=240','c=280');

%计算四个过程功率的总和（正负都考虑）
ComVxQ1*ComFxQ1'+ComVxQ2*ComFxQ2'+ComVxQ3*ComFxQ3'+ComVxQ4*ComFxQ4'
