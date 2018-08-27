function [WorkSection1, WorkSection2,WorkSection3,WorkSection4]=Sectionwork()

load('par.mat');
load('Q1.mat');
load('Q2.mat');
load('Q3.mat');
load('Q4.mat');
L = Subs(1);
A = Subs(2);
T = Subs(3);
w = 2*2*pi/T;
m = Subs(4);
g = Subs(5);
d = Subs(6);
tdeta=Subs(9);

timelengthQ1Q3=length(ComLxQ1)+length(ComLxQ2)+length(ComLxQ3);
timelengthQ2Q4=length(ComLxQ2)+length(ComLxQ3)+length(ComLxQ4);
WorkSection1=0;
WorkSection2=0;
WorkSection3=0;
WorkSection4=0;

LoneFxy=[LoneFxy]


for i=1:timelengthQ1Q3
    
    PLeg=
    
    
    if LtwoPxy(i)<0 && LtwoPxy(i+1)>0
        num=num+1;
    end
    if num==2
        Pushoffprowork=Pushoffprowork+LtwoPxy(i);
    end 
end

