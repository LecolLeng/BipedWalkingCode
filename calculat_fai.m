clc;
clear;
m=20;
c=100;
w=2*3.14*2;
k=1800;
H=0.03;
g=9.8;
M=3;
i=1;
ffshoulder=zeros(1,101); %肩膀受的力，弹性背包
fGshoulder=zeros(1,101); %固定背包
t=0:0.01:1;
fai=faical(m,k,c,w);
angulefai=fai*180/3.14;
zhenfuB=zhenfu(H,m,k,c,w);

for i=1:101
fGshoulder(i)=(m+M)*g-(m+M)*H*w*w*sin(w*t(i));
end

for i=1:101
ffshoulder(i)=fshoulder(m,M,zhenfuB,H,k,c,w,t(i),g,fai);
end

plot(t,ffshoulder)
hold on 
plot(t,fGshoulder)

k=2300;
for i=1:101
ffshoulder(i)=fshoulder(m,M,zhenfuB,H,k,c,w,t(i),g,fai);
end
plot(t,ffshoulder)


w=3*3.14*2;
fai=faical(m,k,c,w);
angulefai2=atan(fai)*180/3.14;


w=4*3.14*2;
fai=faical(m,k,c,w);
angulefai3=atan(fai)*180/3.14;




