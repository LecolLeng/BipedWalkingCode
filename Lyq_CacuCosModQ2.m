%%%����˫��֧���ڣ�Y���������˶���0��pi-aerfa*pi����
function Lyq_CacuCosModQ2()
clc
clear all
%����������ز���,�ֱ�Ϊ1)�ȳ�m��2)�������m��3)��̬����s��4)��������kg��5)�������ٶȡ�6)�벽����7)��ߵ����ĺ����ٶȡ�8)˫��֧�������
%����Q2��������
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
%Q2����Legone,x�����ϳ�ʼ��ֵ����������
LoneLxQ2=zeros(1,length(ttQ2));
LoneLxQ2(1)=LoneLxQ1(end) ;
LoneVxQ2=zeros(1,length(ttQ2));
LoneVxQ2(1)=LoneVxQ1(end);
LoneFxQ2=zeros(1,length(ttQ2));
LoneFxQ2(1)=LoneFxQ1(end);
LoneAxQ2=zeros(1,length(ttQ2));
LoneAxQ2(1)=LoneAxQ1(end);
%Q2����Legone,Y�����ϳ�ʼ��ֵ����������
LoneLyQ2=zeros(1,length(ttQ2));
LoneLyQ2(1)=LoneLyQ1(end);
LoneVyQ2=zeros(1,length(ttQ2));
LoneVyQ2(1)=LoneVyQ1(end) ;
LoneFyQ2=zeros(1,length(ttQ2));
LoneFyQ2(1)=LoneFyQ1(end);
LoneAyQ2=zeros(1,length(ttQ2));
LoneAyQ2(1)=LoneAyQ1(end);
%Q2����Legone,XY�����ϳ�ʼ��ֵ����������
LoneLxyQ2=zeros(1,length(ttQ2));
LoneLxyQ2(1)=LoneLxyQ1(end);
LoneVxyQ2=zeros(1,length(ttQ2));
LoneVxyQ2(1)=LoneVxyQ1(end) ;
LoneFxyQ2=zeros(1,length(ttQ2));
LoneFxyQ2(1)=LoneFxyQ1(end);
LoneAxyQ2=zeros(1,length(ttQ2));
LoneAxyQ2(1)=LoneAxyQ1(end);
%Q2����Legtwo,x�����ϳ�ʼ��ֵ����������
LtwoLxQ2=zeros(1,length(ttQ2));
LtwoLxQ2(1)=d-LoneLxQ1(end);
LtwoVxQ2=zeros(1,length(ttQ2));
LtwoVxQ2(1)=-LoneVxQ1(end);
LtwoFxQ2=zeros(1,length(ttQ2));
LtwoFxQ2(1)=0;
LtwoAxQ2=zeros(1,length(ttQ2));
LtwoAxQ2(1)=-LoneAxQ1(end);
%Q2����Legtwo,Y�����ϳ�ʼ��ֵ����������
LtwoLyQ2=zeros(1,length(ttQ2));
LtwoLyQ2(1)=LoneLyQ1(end);
LtwoVyQ2=zeros(1,length(ttQ2));
LtwoVyQ2(1)=LoneVyQ1(end);
LtwoFyQ2=zeros(1,length(ttQ2));
LtwoFyQ2(1)=0;
LtwoAyQ2=zeros(1,length(ttQ2));
LtwoAyQ2(1)=LoneAyQ1(end);
%Q2����Legtwo,xY�����ϳ�ʼ��ֵ����������
LtwoLxyQ2=zeros(1,length(ttQ2));
LtwoVxyQ2=zeros(1,length(ttQ2));
LtwoFxyQ2=zeros(1,length(ttQ2));
LtwoAxyQ2=zeros(1,length(ttQ2));
%Q2����Legone��Legtwo����,��X������(���Ĵ�)
LonetwoFxQ2=zeros(1,length(ttQ2));
LonetwoFxQ2(1)=LoneFxQ1(end);
%Q2����Legone��Legtwo����,��Y������(���Ĵ�)
LonetwoFyQ2=zeros(1,length(ttQ2));
LonetwoFyQ2(1)=LoneFyQ1(end);
%Q2�������ģ�X����ϵ��,������Ӧ��������
ComLxQ2=zeros(1,length(ttQ2));
ComLxQ2(1)=LoneLxQ2(1);
ComVxQ2=zeros(1,length(ttQ2));
ComVxQ2(1)=LoneVxQ2(1);
ComFxoneQ2=zeros(1,length(ttQ2)); %��һ���ȸ������ĵ�x������
ComFxtwoQ2=zeros(1,length(ttQ2)); %�ڶ����ȸ������ĵ�x������
ComFxQ2=zeros(1,length(ttQ2)); %�������ĵ�x�������
ComAxQ2=zeros(1,length(ttQ2));
%Q2�������ģ�X����ϵ��,������Ӧ��������
ComLyQ2=zeros(1,length(ttQ2));
ComVyQ2=zeros(1,length(ttQ2));
ComFyoneQ2=zeros(1,length(ttQ2));%��һ���ȸ������ĵ�y������
ComFytwoQ2=zeros(1,length(ttQ2)); %�ڶ����ȸ������ĵ�y������
ComFyQ2=zeros(1,length(ttQ2)); %�������ĵ�y�������
ComAyQ2=zeros(1,length(ttQ2));
%Q2�������ģ�XY����ϵ��,������Ӧ��������
ComLxyQ2=zeros(1,length(ttQ2));
ComVxyQ2=zeros(1,length(ttQ2));
ComFxyoneQ2=zeros(1,length(ttQ2)); %��һ���ȸ������ĵ�xy������
ComFxytwoQ2=zeros(1,length(ttQ2)); %�ڶ����ȸ������ĵ�xy������
ComFxyQ2=zeros(1,length(ttQ2)); %�������ĵ�xy�������
ComAxyQ2=zeros(1,length(ttQ2));
%%

%%��ȡ�仯��������
FLQ2end=MSbp*g+MLbp*g-FM*cos(w*ttQ2(end)-forcefai);
QuaFunParaQ2=FittingCurve(1,ttQ20,LtwoFyQ2(1),ttQ2(end),(Q2percent/(Q2percent+Q3percent))*(-m*A*w^2*cos(w*ttQ2(end))+m*g+FLQ2end));

for i=1:length(ttQ2)  %��0����һʱ�̿�ʼ
    
    if i>1
        %%
        %Q2����Legone(������),x�����ϵ���Ӧ����
        tspan = [ttQ20 ttQ2(i)];
        ComLVxQ20 = [ComLxQ2(1) ComVxQ2(1) ];%LoneLxQ20��ʾ��һ��ʱ���ʼ����λ�ƣ�LoneLxQ(1)��ʾ��ʼ�����ٶ�
        [t,y] = ode113(@(t,y) odefunddQ2(t,y,A,L,w,g,QuaFunParaQ2(1),QuaFunParaQ2(2),QuaFunParaQ2(3),Q2percent,Q3percent,T,d,m,MSbp,MLbp,FM,forcefai), tspan, ComLVxQ20);
        ComLxQ2(i)=y(end,1);
        ComVxQ2(i)=y(end,2);
        ComFxQ2(i)=(ComVxQ2(i)-ComVxQ2(i-1))/(ttQ2(i)-ttQ2(i-1))*(m+MSbp+MLbp); %(���Ⱥ���)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ComAxQ2(i)=ComFxQ2(i)/(m+MSbp+MLbp); %(���Ⱥ��������ļ��ٶ�)
        %Q2����Legone,x�����ϵ���Ӧ����
        LoneLxQ2(i)=ComLxQ2(i);
        LoneVxQ2(i)=ComVxQ2(i);
        LoneAxQ2(i)=ComAxQ2(i);
        %Q2����Legone,y�����ϵ���Ӧ����
        LoneLyQ2(i)=L-A+A*cos(w*ttQ2(i));
        LoneVyQ2(i)=-A*w*sin(w*ttQ2(i));
        FL=MSbp*g+MLbp*g-FM*cos(w*ttQ2(i)-forcefai); 
        %//�޸ĵ�
        LtwoFyQ2(i)=QuaFunParaQ2(1)*ttQ2(i)^2+QuaFunParaQ2(2)*ttQ2(i)+QuaFunParaQ2(3);
        LoneFyQ2(i)= -m*A*w^2*cos(w*ttQ2(i))+m*g+FL-LtwoFyQ2(i);
        %LoneFyQ2(i)=-m*A*w^2*cos(w*ttQ2(i))+m*g+FL-(2*ttQ2(i)*(m*A*w^2+m*g+FL)/(aerfa*T))+(1-aerfa)*(m*A*w^2+m*g+FL)/(2*aerfa);
        LoneAyQ2(i)=-A*w^2*cos(w*ttQ2(i));
        %Q2����Legone��������x���������ļ���
        LoneFxQ2(i)=LoneFyQ2(i)*LoneLxQ2(i)/LoneLyQ2(i);
        %Q2����Legtwo,x�����ϵ���Ӧ����
        LtwoLxQ2(i)=d-LoneLxQ2(i);
        LtwoVxQ2(i)=-LoneVxQ2(i);
        LtwoFxQ2(i)=LoneFxQ2(i)-ComFxQ2(i); %����������һ���ȵ���,(LtwoFxQ2(i)����ֵ����)
        LtwoAxQ2(i)=-LoneAxQ2(i);  %���ٶ�Ϊ+
        %Q2����Legtwo,y�����ϵ���Ӧ����
        LtwoLyQ2(i)= LoneLyQ2(i);
        LtwoVyQ2(i)=LoneVyQ2(i);
        %LtwoFyQ2(i)=;
        %LtwoFyQ2(i)=(2*ttQ2(i)*(m*A*w^2+m*g+FL)/(aerfa*T))-(1-aerfa)*(m*A*w^2+m*g+FL)/(2*aerfa);
        LtwoAyQ2(i)=LoneAyQ2(i);
    end
    %Q2����Legone,XY�����ϳ�ʼ��ֵ����������
    LoneLxyQ2(i)=sqrt(LoneLxQ2(i)^2+LoneLyQ2(i)^2);
    LoneVxyQ2(i)=sqrt(LoneVxQ2(i)^2+LoneVyQ2(i)^2);
    LoneFxyQ2(i)=sqrt(LoneFxQ2(i)^2+LoneFyQ2(i)^2);
    LoneAxyQ2(i)=sqrt(LoneAxQ2(i)^2+LoneAyQ2(i)^2);
    %Q2����Legtwo,xY�����ϳ�ʼ��ֵ����������
    LtwoLxyQ2(i)=sqrt(LtwoLxQ2(i)^2+LtwoLyQ2(i)^2);
    LtwoVxyQ2(i)=sqrt(LoneVxQ2(i)^2+LtwoVyQ2(i)^2);
    LtwoFxyQ2(i)=sqrt(LtwoFxQ2(i)^2+LtwoFyQ2(i)^2);
    LtwoAxyQ2(i)=sqrt(LtwoAxQ2(i)^2+LtwoAyQ2(i)^2);
    
    %Q2�������ģ�X����ϵ��,������Ӧ��������
    ComLxQ2(i)=LoneLxQ2(i);
    ComVxQ2(i)=LoneVxQ2(i);
    ComFxoneQ2(i)=LoneFxQ2(i); %��һ���ȸ������ĵ�x������
    ComFxtwoQ2(i)=LtwoFxQ2(i); %�ڶ����ȸ������ĵ�x������
    ComFxQ2(i)=LoneFxQ2(i)-LtwoFxQ2(i);%�������ĵ�x�������
    ComAxQ2(i)=LoneAxQ2(i);
    %Q2�������ģ�X����ϵ��,������Ӧ��������
    ComLyQ2(i)=LoneLyQ2(i);
    ComVyQ2(i)=LoneVyQ2(i);
    ComFyoneQ2(i)=LoneFyQ2(i);%��һ���ȸ������ĵ�y������
    ComFytwoQ2(i)=LtwoFyQ2(i); %�ڶ����ȸ������ĵ�y������
    ComFyQ2(i)=LoneFyQ2(i)+LtwoFyQ2(i); %�������ĵ�y�������
    ComAyQ2(i)=LoneAyQ2(i);
    %Q2�������ģ�XY����ϵ��,������Ӧ��������
    ComLxyQ2(i)=sqrt(ComLxQ2(i)^2+ComLyQ2(i)^2);
    ComVxyQ2(i)=sqrt(ComVxQ2(i)^2+ComVyQ2(i)^2);
    ComFxyoneQ2(i)=LoneFxyQ2(i);
    ComFxytwoQ2(i)=LtwoFxyQ2(i);
    ComFxyQ2(i)=sqrt(ComFxQ2(i)^2+ComFyQ2(i)^2);
    ComAxyQ2(i)=LoneAxyQ2(i);
end


%%%��̬����Ϊ������
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
% ylabel('Acceleration��m/s^2��')
% legend('Lone-ComH','Lone-ComV','Lone-Comxy','LtwoH','LtwoV','Ltwoxy');
% xlim([ttQ2(1)/T,ttQ2(end)/T])

save('Q2','ttQ2','Subs','LoneLxQ2','LoneVxQ2','LoneFxQ2','LoneAxQ2','LoneLyQ2','LoneVyQ2','LoneFyQ2','LoneAyQ2','LoneLxyQ2','LoneVxyQ2','LoneFxyQ2','LoneAxyQ2',...
                        'LtwoLxQ2','LtwoVxQ2','LtwoFxQ2','LtwoAxQ2','LtwoLyQ2','LtwoVyQ2','LtwoFyQ2','LtwoAyQ2','LtwoLxyQ2','LtwoVxyQ2','LtwoFxyQ2','LtwoAxyQ2',...
                        'ComLxQ2','ComVxQ2','ComFxoneQ2','ComFxtwoQ2','ComFxQ2','ComAxQ2', 'ComLyQ2','ComVyQ2','ComFyoneQ2','ComFytwoQ2','ComFyQ2','ComAyQ2','ComLxyQ2','ComVxyQ2','ComFxyoneQ2','ComFxytwoQ2','ComFxyQ2','ComAxyQ2');
