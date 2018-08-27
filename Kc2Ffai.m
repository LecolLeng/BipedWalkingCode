%刚度变化，不同的阻尼，绘制作用在肩膀上力的幅值变化，以及相位变化
clc;
clear;
m=13.5;
w=2.4*3.14*2;
H=0.03;
g=9.8;
M=3;
FMshoulder=zeros(1,8000); %肩膀受的力，弹性背包
fai=zeros(1,8000);
zhenfuB=zeros(1,8000);
k=1:8000;


c=5;
for i=1:8000
   fai(i)=faical(m,k(i),c,w);
   zhenfuB(i)=zhenfu(H,m,k(i),c,w);
   FMshoulder(i)=FMPack(M,m, zhenfuB(i),w,H,fai(i));
end
figure(1);
plot(k, FMshoulder)



c=100;
for i=1:8000
   fai(i)=faical(m,k(i),c,w);
   zhenfuB(i)=zhenfu(H,m,k(i),c,w);
   FMshoulder(i)=FMPack(M,m, zhenfuB(i),w,H,fai(i));
end

hold on 
plot(k, FMshoulder)



c=200;
for i=1:8000
   fai(i)=faical(m,k(i),c,w);
   zhenfuB(i)=zhenfu(H,m,k(i),c,w);
   FMshoulder(i)=FMPack(M,m, zhenfuB(i),w,H,fai(i));
end

hold on 
plot(k, FMshoulder)



c=300;
for i=1:8000
   fai(i)=faical(m,k(i),c,w);
   zhenfuB(i)=zhenfu(H,m,k(i),c,w);
   FMshoulder(i)=FMPack(M,m, zhenfuB(i),w,H,fai(i));
end

hold on 
plot(k, FMshoulder)

c=400;
for i=1:8000
   fai(i)=faical(m,k(i),c,w);
   zhenfuB(i)=zhenfu(H,m,k(i),c,w);
   FMshoulder(i)=FMPack(M,m, zhenfuB(i),w,H,fai(i));
end

hold on 
plot(k, FMshoulder)
xlabel('刚度(N/m)')
ylabel('力幅值(N)')
legend('c=50','c=100','c=200','c=300','c=400');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c=30;
for i=1:8000
   fai(i)=faical(m,k(i),c,w);
   zhenfuB(i)=zhenfu(H,m,k(i),c,w);
   FMshoulder(i)=FMPack(M,m, zhenfuB(i),w,H,fai(i));
   if fai(i)<0
       fai(i)=3.14+fai(i);
   end
end

figure(2);
plot(k,fai*180/3.14);


c=100;
for i=1:8000
   fai(i)=faical(m,k(i),c,w);
   zhenfuB(i)=zhenfu(H,m,k(i),c,w);
   FMshoulder(i)=FMPack(M,m, zhenfuB(i),w,H,fai(i));
   if fai(i)<0
       fai(i)=3.14+fai(i);
   end
end
hold on 
plot(k,fai*180/3.14);


c=200;
for i=1:8000
   fai(i)=faical(m,k(i),c,w);
   zhenfuB(i)=zhenfu(H,m,k(i),c,w);
   FMshoulder(i)=FMPack(M,m, zhenfuB(i),w,H,fai(i));
   if fai(i)<0
       fai(i)=3.14+fai(i);
   end
end
plot(k,fai*180/3.14);


c=300;
for i=1:8000
   fai(i)=faical(m,k(i),c,w);
   zhenfuB(i)=zhenfu(H,m,k(i),c,w);
   FMshoulder(i)=FMPack(M,m, zhenfuB(i),w,H,fai(i));
   if fai(i)<0
       fai(i)=3.14+fai(i);
   end
end

plot(k,fai*180/3.14);


c=400;
for i=1:8000
   fai(i)=faical(m,k(i),c,w);
   zhenfuB(i)=zhenfu(H,m,k(i),c,w);
   FMshoulder(i)=FMPack(M,m, zhenfuB(i),w,H,fai(i));
   if fai(i)<0
       fai(i)=3.14+fai(i);
   end
end

plot(k,fai*180/3.14);
legend('c=5','c=100','c=200','c=300','c=400');
xlabel('刚度(N/m)')
ylabel('相位差(°)')

