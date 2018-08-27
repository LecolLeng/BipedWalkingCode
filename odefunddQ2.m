function dydt = odefunddQ2(t,y,A,L,w,g,QuaFunParaQ2_1,QuaFunParaQ2_2,QuaFunParaQ2_3,Q2percent,Q3percent,T,d,m,MSbp,MLbp,FM,forcefai)
dydt = zeros(2,1);
dydt(1) = y(2);
coswt=1-0.5*(w*t)^2+(w*t)^4/(4*3*2)-(w*t)^6/(6*5*4*3*2)+(w*t)^8/(8*7*6*5*4*3*2);
FL=MSbp*g+MLbp*g-FM*cos(w*t-forcefai);
% FL_Q2E=MSbp*g+MLbp*g-FM*cos(w*T/4-forcefai);
% t1Q2=(1-2*Q2percent)*T/4;%Q2初始时刻点的时间
% t2Q2=T/4;%Q2结束时刻点的时间
% Fl2yQ2E=Q2percent*(A*m*w^2+m*g+FL_Q2E)/(Q2percent+Q3percent);
% 
% al2q2=Fl2yQ2E/(t1Q2-t2Q2)^2; 
% bl2q2=-2*al2q2*t1Q2;
% cl2q2=al2q2*t1Q2^2;
al2q2=QuaFunParaQ2_1;
bl2q2=QuaFunParaQ2_2;
cl2q2=QuaFunParaQ2_3;

fenzi1=-A*m*w^2*coswt+m*g+FL-al2q2*t^2-bl2q2*t-cl2q2;
fenmu1=L-A+A*coswt;
fenzi2=al2q2*t^2+bl2q2*t+cl2q2;
fenmu2=fenmu1;
dydt(2) =(fenzi1*y(1)/(fenmu1)-(d-y(1))*fenzi2/(fenmu2))/(m+MSbp+MLbp);

