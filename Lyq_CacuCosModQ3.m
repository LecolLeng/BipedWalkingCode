%%%����˫��֧���ڣ�Y���������˶���0��pi-aerfa*pi����
function Lyq_CacuCosModQ3()
clc
clear all
%����������ز���,�ֱ�Ϊ1)�ȳ�m��2)�������m��3)��̬����s��4)��������kg��5)�������ٶȡ�6)�벽����7)��ߵ����ĺ����ٶȡ�8)˫��֧�������
%����Q3��������
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
%Q3����Legone,x�����ϳ�ʼ��ֵ����������
LoneLxQ3=zeros(1,length(ttQ3));
LoneLxQ3(1)=LoneLxQ2(end);
LoneVxQ3=zeros(1,length(ttQ3));
LoneVxQ3(1)=LoneVxQ2(end);
LoneFxQ3=zeros(1,length(ttQ3));
LoneFxQ3(1)=LoneFxQ2(end);
LoneAxQ3=zeros(1,length(ttQ3));
LoneAxQ3(1)=LoneAxQ2(end);
%Q3����Legone,Y�����ϳ�ʼ��ֵ����������
LoneLyQ3=zeros(1,length(ttQ3));
LoneLyQ3(1)=LoneLyQ2(end);
LoneVyQ3=zeros(1,length(ttQ3));
LoneVyQ3(1)=LoneVyQ2(end) ;
LoneFyQ3=zeros(1,length(ttQ3));
LoneFyQ3(1)=LoneFyQ2(end);
LoneAyQ3=zeros(1,length(ttQ3));
LoneAyQ3(1)=LoneAyQ2(end);
%Q3����Legone,XY�����ϳ�ʼ��ֵ����������
LoneLxyQ3=zeros(1,length(ttQ3));
LoneLxyQ3(1)=LoneLxyQ2(end);
LoneVxyQ3=zeros(1,length(ttQ3));
LoneVxyQ3(1)=LoneVxyQ2(end) ;
LoneFxyQ3=zeros(1,length(ttQ3));
LoneFxyQ3(1)=LoneFxyQ2(end);
LoneAxyQ3=zeros(1,length(ttQ3));
LoneAxyQ3(1)=LoneAxyQ2(end);
%Q3����Legtwo,x�����ϳ�ʼ��ֵ����������
LtwoLxQ3=zeros(1,length(ttQ3));
LtwoLxQ3(1)=LtwoLxQ2(end);
LtwoVxQ3=zeros(1,length(ttQ3));
LtwoVxQ3(1)=LtwoVxQ2(end);
LtwoFxQ3=zeros(1,length(ttQ3));
LtwoFxQ3(1)=LtwoFxQ2(end);
LtwoAxQ3=zeros(1,length(ttQ3));
LtwoAxQ3(1)=LtwoAxQ2(end);
%Q3����Legtwo,Y�����ϳ�ʼ��ֵ����������
LtwoLyQ3=zeros(1,length(ttQ3));
LtwoLyQ3(1)=LtwoLyQ2(end);
LtwoVyQ3=zeros(1,length(ttQ3));
LtwoVyQ3(1)=LtwoVyQ2(end);
LtwoFyQ3=zeros(1,length(ttQ3));
LtwoFyQ3(1)=LtwoFyQ2(end);
LtwoAyQ3=zeros(1,length(ttQ3));
LtwoAyQ3(1)=LtwoAyQ2(end);
%Q3����Legtwo,xY�����ϳ�ʼ��ֵ����������
LtwoLxyQ3=zeros(1,length(ttQ3));
LtwoLxyQ3(1)=LtwoLxyQ2(end);
LtwoVxyQ3=zeros(1,length(ttQ3));
LtwoVxyQ3(1)=LtwoVxyQ2(end);
LtwoFxyQ3=zeros(1,length(ttQ3));
LtwoFxyQ3(1)=LtwoFxyQ2(end);
LtwoAxyQ3=zeros(1,length(ttQ3));
LtwoAxyQ3(1)=LtwoAxyQ2(end);

%Q3����Legone��Legtwo����,��X������(���Ĵ�)
LonetwoFxQ3=zeros(1,length(ttQ3));
LonetwoFxQ3(1)=LoneFxQ2(end)-LtwoFxQ2(end);
LonetwoFyQ3=zeros(1,length(ttQ3));
LonetwoFyQ3(1)=LoneFyQ2(end)+LtwoFyQ2(end);

%Q3�������ģ�X����ϵ��,������Ӧ��������
ComLxQ3=zeros(1,length(ttQ3));
ComLxQ3(1)=ComLxQ3(1);
ComVxQ3=zeros(1,length(ttQ3));
ComVxQ3(1)=ComVxQ3(1);
ComFxoneQ3=zeros(1,length(ttQ3)); %��һ���ȸ������ĵ�x������
ComFxtwoQ3=zeros(1,length(ttQ3)); %�ڶ����ȸ������ĵ�x������
ComFxQ3=zeros(1,length(ttQ3)); %�������ĵ�x�������
ComFxQ3(1)=ComFxQ2(end);
ComAxQ3=zeros(1,length(ttQ3));
ComAxQ3(1)=ComAxQ2(end);
%Q3�������ģ�X����ϵ��,������Ӧ��������
ComLyQ3=zeros(1,length(ttQ3));
ComVyQ3=zeros(1,length(ttQ3));
ComFyoneQ3=zeros(1,length(ttQ3));%��һ���ȸ������ĵ�y������
ComFytwoQ3=zeros(1,length(ttQ3)); %�ڶ����ȸ������ĵ�y������
ComFyQ3=zeros(1,length(ttQ3)); %�������ĵ�y�������
ComAyQ3=zeros(1,length(ttQ3));
%Q3����XY����ϵ��,������Ӧ��������
ComLxyQ3=zeros(1,length(ttQ3));
ComVxyQ3=zeros(1,length(ttQ3));
ComFxyoneQ3=zeros(1,length(ttQ3));
ComFxytwoQ3=zeros(1,length(ttQ3));
ComFxyQ3=zeros(1,length(ttQ3));
ComAxyQ3=zeros(1,length(ttQ3));
%%
FLQ3end=MSbp*g+MLbp*g-FM*cos(w*ttQ3(end)-forcefai);
QuaFunParaQ3=FittingCurve(2,ttQ30,LtwoFyQ3(1),ttQ3(end),-m*A*w^2*cos(w*ttQ3(end))+m*g+FLQ3end);
%%%���ǰ��յڶ����ȵĽ��м���

for i=1:length(ttQ3)  %��0����һʱ�̿�ʼ
    
    if i>1
        %%
        %Q3����Legone(������),x�����ϵ���Ӧ����
        tspan = [ttQ30 ttQ3(i)];
        ComLVxQ30 = [LoneLxQ3(1) LoneVxQ3(1) ];%LoneLxQ30��ʾ��һ��ʱ���ʼ����λ�ƣ�LoneLxQ(1)��ʾ��ʼ�����ٶ�
        [t,y] = ode113(@(t,y) odefunddQ3(t,y,A,L,w,g,QuaFunParaQ3(1),QuaFunParaQ3(2),QuaFunParaQ3(3),T,d,m,MSbp,MLbp,FM,forcefai), tspan, ComLVxQ30);
        ComLxQ3(i)=y(end,1);
        ComVxQ3(i)=y(end,2);
        ComFxQ3(i)=(ComVxQ3(i)-ComVxQ3(i-1))/(ttQ3(i)-ttQ3(i-1))*(m+MSbp+MLbp); %(���Ⱥ���)
        ComAxQ3(i)=(ComVxQ3(i)-ComVxQ3(i-1))/(ttQ3(i)-ttQ3(i-1)); %(���Ⱥ��������ļ��ٶ�)
        %Q3����Legone,x�����ϵ���Ӧ����
        LoneLxQ3(i)=ComLxQ3(i);
        LoneVxQ3(i)=ComVxQ3(i);
        LoneAxQ3(i)=ComAxQ3(i);
        %Q3����Legtwo,y��������
        LtwoFyQ3(i)=QuaFunParaQ3(1)*ttQ3(i)^2+QuaFunParaQ3(2)*ttQ3(i)+QuaFunParaQ3(3);
        
        %Q3����Legone,y�����ϵ���Ӧ����
        LoneLyQ3(i)=L-A+A*cos(w*ttQ3(i));
        LoneVyQ3(i)=-A*w*sin(w*ttQ3(i));
        FL=MSbp*g+MLbp*g-FM*cos(w*ttQ3(i)-forcefai);
        LoneFyQ3(i)=-m*A*w^2*cos(w*ttQ3(i))+m*g+FL-LtwoFyQ3(i);
        %LoneFyQ3(i)=-(2*(m*A*w^2+m*g+FL)*ttQ3(i))/(aerfa*T)+(1+aerfa)*(m*A*w^2+m*g+FL)/(2*aerfa);
        LoneAyQ3(i)=-A*w^2*cos(w*ttQ3(i));
        %Q3����Legone��������x���������ļ���
        LoneFxQ3(i)=LoneFyQ3(i)*LoneLxQ3(i)/LoneLyQ3(i);
        %Q3����Legtwo,x�����ϵ���Ӧ����
        LtwoLxQ3(i)=d-LoneLxQ3(i);
        LtwoVxQ3(i)=-LoneVxQ3(i);
        LtwoFxQ3(i)=LoneFxQ3(i)-ComFxQ3(i); %����������һ���ȵ���
        LtwoAxQ3(i)=-LoneAxQ3(i);
        %Q3����Legtwo,y�����ϵ���Ӧ����
        LtwoLyQ3(i)= LoneLyQ3(i);
        LtwoVyQ3(i)=LoneVyQ3(i);
        %LtwoFyQ3(i)=-m*A*w^2*cos(w*ttQ3(i))+m*g+FL+(2*ttQ3(i)*(m*A*w^2+m*g+FL)/(aerfa*T))-(1+aerfa)*(m*A*w^2+m*g+FL)/(2*aerfa);
        LtwoAyQ3(i)=LoneAyQ3(i);    

    end
    %Q3����Legone,XY�����ϳ�ʼ��ֵ����������
    LoneLxyQ3(i)=sqrt(LoneLxQ3(i)^2+LoneLyQ3(i)^2);
    LoneVxyQ3(i)=sqrt(LoneVxQ3(i)^2+LoneVyQ3(i)^2);
    LoneFxyQ3(i)=sqrt(LoneFxQ3(i)^2+LoneFyQ3(i)^2);
    LoneAxyQ3(i)=sqrt(LoneAxQ3(i)^2+LoneAyQ3(i)^2);
    %Q3����Legtwo,xY�����ϳ�ʼ��ֵ����������
    LtwoLxyQ3(i)=sqrt(LtwoLxQ3(i)^2+LtwoLyQ3(i)^2);
    LtwoVxyQ3(i)=sqrt(LoneVxQ3(i)^2+LtwoVyQ3(i)^2);
    LtwoFxyQ3(i)=sqrt(LtwoFxQ3(i)^2+LtwoFyQ3(i)^2);
    LtwoAxyQ3(i)=sqrt(LtwoAxQ3(i)^2+LtwoAyQ3(i)^2);
    
    %Q3�������ģ�X����ϵ��,������Ӧ��������
    ComLxQ3(i)=LoneLxQ3(i);
    ComVxQ3(i)=LoneVxQ3(i);
    ComFxoneQ3(i)=LoneFxQ3(i); %��һ���ȸ������ĵ�x������
    ComFxtwoQ3(i)=LtwoFxQ3(i); %�ڶ����ȸ������ĵ�x������
    %ComFxQ3(i)=LoneFxQ3(i)-LtwoFxQ3(i); %�������ĵ�x�������
    %ComAxQ3(i)=LoneAxQ3(i);
    %Q3�������ģ�X����ϵ��,������Ӧ��������
    ComLyQ3(i)=LoneLyQ3(i);
    ComVyQ3(i)=LoneVyQ3(i);
    ComFyoneQ3(i)=LoneFyQ3(i);%��һ���ȸ������ĵ�y������
    ComFytwoQ3(i)=LtwoFyQ3(i); %�ڶ����ȸ������ĵ�y������
    ComFyQ3(i)=LoneFyQ3(i)+LtwoFyQ3(i); %�������ĵ�y�������
    ComAyQ3(i)=LoneAyQ3(i);
    %Q3�������ģ�XY����ϵ��,������Ӧ��������
    ComLxyQ3(i)=sqrt(ComLxQ3(i)^2+ComLyQ3(i)^2);
    ComVxyQ3(i)=sqrt(ComVxQ3(i)^2+ComVyQ3(i)^2);
    ComFxyoneQ3(i)=LoneFxyQ3(i);
    ComFxytwoQ3(i)=LtwoFxyQ3(i);
    ComFxyQ3(i)=sqrt(ComFxQ3(i)^2+ComFyQ3(i)^2);
    ComAxyQ3(i)=LoneAxyQ3(i); 
end


%%�Բ�̬����Ϊ������ 
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
% ylabel('Acceleration��m/s^2��')
% legend('Lone-ComH','Lone-ComV','Lone-Comxy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ3(1)/T,ttQ3(end)/T])

save('Q3','ttQ3','Subs','LoneLxQ3','LoneVxQ3','LoneFxQ3','LoneAxQ3','LoneLyQ3','LoneVyQ3','LoneFyQ3','LoneAyQ3','LoneLxyQ3','LoneVxyQ3','LoneFxyQ3','LoneAxyQ3',...
                        'LtwoLxQ3','LtwoVxQ3','LtwoFxQ3','LtwoAxQ3','LtwoLyQ3','LtwoVyQ3','LtwoFyQ3','LtwoAyQ3','LtwoLxyQ3','LtwoVxyQ3','LtwoFxyQ3','LtwoAxyQ3',...
                        'ComLxQ3','ComVxQ3','ComFxoneQ3','ComFxtwoQ3','ComFxQ3','ComAxQ3', 'ComLyQ3','ComVyQ3','ComFyoneQ3','ComFytwoQ3','ComFyQ3','ComAyQ3','ComLxyQ3','ComVxyQ3','ComFxyoneQ3','ComFxytwoQ3','ComFxyQ3','ComAxyQ3');