%%%
clc
clear
Step=0.650:0.001:0.650;
Stepdeta=zeros(1,length(Step));
Axend=zeros(1,length(Step));
for i=1:length(Step)    %%Ϊ��ƥ����ʵĲ�����ʹ�趨�Ĳ����ͻ��ֳ��Ĳ����������
    %����������ز���,�ֱ�Ϊ1)�ȳ�m��2)�������m��3)��̬����s��4)��������kg��
    %%5)�������ٶȡ�6)�벽��(�ο�ֵ)��7)��ߵ����ĺ����ٶȡ�8)˫��֧������� 9)ʱ����
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
            %������صĲ���1�������̶���������2���������Ը�������3���ն�4������
            Subsbackpack=[0.02,20,KK(j),CC(k)];
            MSbp=Subsbackpack(1);
            MLbp=Subsbackpack(2);
            K=Subsbackpack(3);
            C=Subsbackpack(4);
            %��Ҫ�õ�ʲô??   %���ķ�ֵFM���˶���Ե���λ��fai�������ڼ粿����������˶�����λ��forcefai�����ؾ����񶯵����dBabs
            [FM,fai,forcefai,Babs]=elasticbackpack(Subsbackpack,A,w,g);
            ttt=0:0.05:T;
            %     %��һ�°��͸��ص����λ�Ƶķ�ֵ����ΪԼ������
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


%���ڼ�ⲽ���Ƿ����
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
xlabel('�ն�(N/m)')
ylabel('Wower(J)')
legend('c=80','c=120','c=160','c=200','c=240','c=280');

%�����ĸ����̹��ʵ��ܺͣ����������ǣ�
ComVxQ1*ComFxQ1'+ComVxQ2*ComFxQ2'+ComVxQ3*ComFxQ3'+ComVxQ4*ComFxQ4'
