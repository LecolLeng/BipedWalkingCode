%%%%%用于示意运动及力的相关曲线
clc
clear
t=0:0.01:4*pi;
fai=0.7*pi;
plot(t,cos(t),'--b')
hold on
plot(t,cos(t-fai),'--r')
hold on
plot(t,-cos(t),'-b')
hold on
plot(t,-cos(t-fai),'-r')
legend('Commotion','Loadmotion','ComForce','LoadForce');
hold on
plot([0.76*pi, 0.76*pi],[1,-1])
hold on
plot([1.24*pi, 1.24*pi],[1,-1])

hold on
plot([2.76*pi, 2.76*pi],[1,-1])
hold on
plot([3.24*pi, 3.24*pi],[1,-1])

