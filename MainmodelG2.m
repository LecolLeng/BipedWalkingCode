%% Editor:Yuquan Leng
%%%�ڶ��棬ģ�����ݣ���֪�������ߡ�������ƽ�����ٵ�
clc
clear
Gra_acc=9.8; %�������ٶ�
g=Gra_acc;
Cyclenum=0;
%% �趨����Ĺ��в�������ߡ����ء�ƽ������
Human_heigh=1.77; %��λm
Human_weight=70.8;%��λKg
Human_Velocity=1.7;%����ƽ���ٶȣ���λm/s
Human_sex=0;%�� 0 Ů1
%% �������岽��ʱ���в������Ƶ���
Vr=Human_Velocity/Human_heigh;%����ٶ�
Gait_frequence=64.8*(Vr^0.57)/60;%���岽Ƶ �ο����ף�Grieve D W, Gear R J. The relationships between length of stride, step frequency, time of swing and speed of walking for children and adults[J]. Ergonomics, 1966, 9(5): 379-399.
Gait_cycle_time=1/Gait_frequence;%��̬���ڣ�
Gait_step_length=Human_Velocity*Gait_cycle_time; %���� Stride length���ο����ף�
Gait_step_half_length=Gait_step_length/2; %�������ο����ף�
%***********�ȳ�(m)*********�ο����ף�Roebuck, J.A., Kroemer, K.H.E., Thomson, W.G.: Engineering Anthropometry Methods. Wiley-interscience, New York (1975)
Leg_real_length=0.53*Human_heigh;
%***********�ų�(cm)*********�ο����ף�Agnihotri A K, Purwar B, Googoolye K, et al. Estimation of stature by foot length[J]. Journal of forensic and legal medicine, 2007, 14(5): 279-283.
if Human_sex>0
    Foot_length=(Human_heigh*100-77.059)/3.536/100; %��λΪm
else
    Foot_length=(Human_heigh*100-68.586)/4.0636/100; %��λΪm
end
%*****˫��֧����ռ������̬���ڵı���*********�ο����ף�Grieve D W, Gear R J. The relationships between length of stride, step frequency, time of swing and speed of walking for children and adults[J]. Ergonomics, 1966, 9(5): 379-399.
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

%*******���������������λm*****************�ο����ף�Hoover J, Meguid S A. Performance assessment of the suspended-load backpack[J]. International Journal of Mechanics and Materials in Design, 2011, 7(2): 111-121.
v=Human_Velocity;
Com_amplitude=(Leg_real_length/2)*(1-(1-(0.936*v/(Leg_real_length*2*1.504*(v/Leg_real_length)^0.57))^2)^0.5);

%% ���������������ת��ģ�͵��ȳ����벽�����ų���˫��֧����ռ�ȡ���ʵ���ȳ���
%���幫ʽ�Ƶ����ο��ĵ�(�������˷�ʽ)
Leg_virtual_length=(-Com_amplitude*Foot_length-Foot_length*Com_amplitude*...
    cos(Double_leg_support_time_rate*pi)+Leg_real_length*Gait_step_half_length*...
    (1+Double_leg_support_time_rate))/(Gait_step_half_length*(1+Double_leg_support_time_rate)-Foot_length);
%���������˷�ʽ
%Leg_virtual_length=0.53*Human_heigh;
%% ��ʼ�ٶȣ���ʼ�ٶ�С��ƽ���ٶȣ���Q2˫��֧������ռ��˫��֧���ڵı�����������Ϊ�Ա���
%Ŀ�꺯��Ϊ�������ܺͣ���������ȵģ�ֻ��ҪԼ������ļ��ɣ�������Ĳ�����ʵ�ʲ����Ĳ�ֵ�������ߵľ���ֵӦ�õ���0�����������
%�趨Ŀ�꺯��ֵС��һ����ֵdeita��������Ϊ���ҵ��˽�
%��������һ����ά��������С������⣬������ʽδ֪�����ܴ��ڼ�ֵ�㣬��ˣ���Ҫ����һ�����ֶν��
%��ʱԤʹ�õķ���Ϊ�����ö����ʼ�㣬���ν���Ѱ�ң�Ѱ�ҵ��˽⣬������ѭ�������û����һ����ʼ�㣨ֻҪ��ʼ���㹻�࣬��ض����ҵ��⣩
%%
A=Com_amplitude;
w=2*2*pi*Gait_frequence; %������Ƶ��
T=1/Gait_frequence;%��̬����

%�趨�Ա�����Χ
Max_Ini_speed=0.95*Human_Velocity;  %��ʼ�ٶȿ϶���С��ƽ���ٶ�
Min_Ini_speed=0.85*Human_Velocity;
MaxQ2_Double_support_percent=0.51*Double_leg_support_time_rate;
MinQ2_Double_support_percent=0.49*Double_leg_support_time_rate;
%scale=(Max_Ini_speed-Min_Ini_speed)/(MaxQ2_Double_support_percent-MinQ2_Double_support_percent); %���з�Χ�ı������㣬���ڲ�ͬ�����ϲ����ĶԵ�
scale=4;
%%%%���Ա����ĸնȼ�������趨
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
        %������صĲ���1�������̶���������2���������Ը�������3���ն�4������
        Subsbackpack=[0.01,0.01,KK(j),CC(k)];
        MSbp=Subsbackpack(1);
        MLbp=Subsbackpack(2);
        K=Subsbackpack(3);
        C=Subsbackpack(4);
        %���ķ�ֵFM���˶���Ե���λ��fai�������ڼ粿����������˶�����λ��forcefai�����ؾ����񶯵����dBabs
        [FM,fai,forcefai,Babs]=elasticbackpack(Subsbackpack,A,w,g);
        %%%%%�д���һ���Ż�%%%%
        ttt=0:0.05:T;
        %%��һ�°��͸��ص����λ�Ƶķ�ֵ����ΪԼ������
        Brel=Babs*cos(w*ttt-fai)-A*cos(w*ttt);
        %  figure(1)
        %  plot(ttt,Brel)
        MaxBrel(j,k)=max(Brel);
        %figure(2)
        FL=MSbp*g+MLbp*g-FM*cos(w*ttt-forcefai);
        %plot(ttt,FLKKCC)
        MaxFL(j,k)=max(FL);
        if MaxBrel(j,k)<0.2   %���������С��0.2m
            %�ڶ���������������nXn�ĳ�ʼ��,�ڸ���ʼ��λ�����е�һ���ڣ������������г�nxn�飩
            n=3;
            Time_gap=0.001; %ÿ�λ��ֵ�ʱ����
            Objdeta=30;
            for ni=1:n
                for nj=1:n
                    
                    %��ĳһ��������ڣ����ѡȡ��ʼ�ٶȺ�Q2֧����ռ˫��֧�Űٷֱȵĳ�ʼֵ
                    Smallarea_speed=(Max_Ini_speed-Min_Ini_speed)/n;
                    Gone_Def_Ini_speed= unifrnd(Min_Ini_speed+Smallarea_speed*(ni-1),Min_Ini_speed+Smallarea_speed*ni);
                    Smallarea_percent=(MaxQ2_Double_support_percent-MinQ2_Double_support_percent)/n;
                    Gone_Def_Q2percent=unifrnd(MinQ2_Double_support_percent+Smallarea_percent*(ni-1),MinQ2_Double_support_percent+Smallarea_percent*ni);
                    
                    %����������ز���,�ֱ�Ϊ1)�ȳ�m��2)�������m��3)��̬����s��4)��������kg��
                    %%5)�������ٶȡ�6)�벽��(�ο�ֵ)��7)��ߵ����ĺ����ٶȡ�8)˫��֧������� 9)ʱ����
                    %%10��Q2˫��֧����ռ�����ı���
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
                    Evo_step_length=0.04; %�趨��ʼ���ݻ�����
                    Numdirection=5;
                    ObjectNew=zeros(1,Numdirection);
                    while (flage) %���ѭ������Ҫ�ҵ���ֵ�㣬���缫ֵ�����ֵΪ0������Ϊ����ֵ�㣬��������ѭ����������
                        
                        %����һ���㣬�����ܱ�4��������չ��ֱ�Ӽ���5�����
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
                            ObjectNew(Directioni)= abs(Gait_step_half_length-ComLxQ4(end))*200+abs(Hor_Psum*Time_gap);%����10��Ŀ�ģ���Ϊ�˵���Ӱ��߶ȡ� %%%�����������ܺ�
                            Cyclenum=Cyclenum+1;
                        end
                        % 5��Ŀ�����ݣ�����ѡȡ�ݶ����ģ������һ����λ�ڼ�ֵ�㣬�������ݶ�Ϊ��
                        %�����Ѿ��Ǽ�ֵ�㣬Ӧ������߽�����������Ŀ�꺯��ֵ�ӽ�0������
                        %�Ѿ���5�����еļ�ֵ�㣬���ǲ�����Ŀ��ӽ�0������������������һ���Ż���
                        %��һ���Ż�����Ȼ���У���˵�������˷�0�ļ�ֵ�㣬ѡȡ������ʼ�㣬��������
                        
                        %��ȡObjectNew�е�5����ֵ��ǰ4���Ĳ�ֵ����ֵԽ��Խ�ã�˵���ݶȴ�
                        Gradient_obj=ObjectNew(Numdirection)*ones(1,length(ObjectNew)-1)-ObjectNew(1,1: (Numdirection-1));
                        [Max_Gradient,Max_Gradient_num]=max(Gradient_obj); %�õ�����ݶ�ֵ�������
                        %�������ݶ�ֵ�Ǹ�����˵���ĵ��Ѿ�����5�����еļ�ֵ��
                        if Max_Gradient<=0
                            %�ݶ�֮��С��0��������Ӧ���ж�
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
                                    if  ObjectNew(5)<Objdeta%%%%%%%%�ҵ����ˣ���ʾ������
                                        flage=0;
                                        ni=4;
                                        nj=4;
                                        %�ҵ��˺��ʵĲ������������Ӧ�ļ���ֵ
                                        [sumwork(j,k), Pushoffprowork(j,k),WorkSection1(j,k), WorkSection2(j,k),WorkSection3(j,k),WorkSection4(j,k)] =showdataandfig();
                                        MaxHeelStr(j,k)=LtwoFxyQ3(end);
                                        MaxHeelPush(j,k)=LoneFxyQ1(end);
                                        Changfai(j,k)=fai; 
                                        ObjectNewAll(j,k)=ObjectNew(5); %%�洢����ֵ
                                        break;%%%%%%%%�ҵ����ˣ�
                                    else
                                        %�����һ�����ҵ��˼�Сֵ��û���ҵ���Сֵ
                                        sprintf('%s��%d �����ҵ��˼�Сֵ��û���ҵ���Сֵ��������Ѱ��',ni,nj);
                                    end
                            end
                        else
                            %���ݶ�ֵ����0�ģ�������ݶȵ����ֵ��ΪNewpoint�������ݶ��½���
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
ylabel('Degree(��)')
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