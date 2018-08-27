%�����仯��ĳһ�����£����������ڼ�������ķ�ֵ�仯���Լ���λ�仯
clc
clear
w=2*3.14*2;
H=0.03;
g=9.8;
M=3;
FMshoulder=zeros(1,20); %����ܵ��������Ա���
fai=zeros(1,20);
zhenfuB=zeros(1,20);
c=50;
m=20:39;

k=2000;
for i=1:20
   fai(i)=faical(m(i),k,c,w);
   zhenfuB(i)=zhenfu(H,m(i),k,c,w);
   FMshoulder(i)=FMPack(M,m(i),zhenfuB(i),w,H,fai(i));
   if fai(i)<0
       fai(i)=3.14+fai(i);
   end
end
figure(1);
plot(m,fai*180/3.14);
hold on

k=1500;
for i=1:20 
   fai(i)=faical(m(i),k,c,w);
   zhenfuB(i)=zhenfu(H,m(i),k,c,w);
   FMshoulder(i)=FMPack(M,m(i),zhenfuB(i),w,H,fai(i));
   if fai(i)<0
       fai(i)=3.14+fai(i);
   end
end
hold on
figure(1);
plot(m,fai*180/3.14)
xlabel('����(Kg)')
ylabel('�Ƕ�(N)')
legend('k=2000','k=1500');


k=2000;
for i=1:20
   fai(i)=faical(m(i),k,c,w);
   zhenfuB(i)=zhenfu(H,m(i),k,c,w);
   FMshoulder(i)=FMPack(M,m(i), zhenfuB(i),w,H,fai(i));
end
figure(2);
plot(m, FMshoulder)
hold on

k=1500;
for i=1:20
   fai(i)=faical(m(i),k,c,w);
   zhenfuB(i)=zhenfu(H,m(i),k,c,w);
   FMshoulder(i)=FMPack(M,m(i), zhenfuB(i),w,H,fai(i));
end
figure(2);
plot(m, FMshoulder)
xlabel('����(Kg)')
ylabel('����ֵ(N)')
legend('k=2000','k=1500');



