%质量变化，某一阻尼下，绘制作用在肩膀上力的幅值变化，以及相位变化
clc
clear
w=2*3.14*2;
H=0.03;
g=9.8;
M=3;
FMshoulder=zeros(1,20); %肩膀受的力，弹性背包
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
xlabel('质量(Kg)')
ylabel('角度(N)')
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
xlabel('质量(Kg)')
ylabel('力幅值(N)')
legend('k=2000','k=1500');



