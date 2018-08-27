function dydt = odefunddQ3(t,y,A,L,w,g,QuaFunPara1,QuaFunPara2,QuaFunPara3,T,d,m,MSbp,MLbp,FM,forcefai)
dydt = zeros(2,1);
dydt(1) = y(2);
coswt=1-0.5*(w*(t-T/2))^2+(w*(t-T/2))^4/(4*3*2)-(w*(t-T/2))^6/(6*5*4*3*2)+(w*(t-T/2))^8/(8*7*6*5*4*3*2);
FL=MSbp*g+MLbp*g-FM*cos(w*t-forcefai);

al2Q3=QuaFunPara1;
bl2Q3=QuaFunPara2;
cl2Q3=QuaFunPara3;

fenzi1=-m*A*w^2*coswt+m*g+FL-al2Q3*t^2-bl2Q3*t-cl2Q3;
fenmu1=L-A+A*coswt;
fenzi2=al2Q3*t^2+bl2Q3*t+cl2Q3;
fenmu2=L-A+A*coswt;
dydt(2) =(fenzi1*y(1)/fenmu1-(d-y(1))*fenzi2/fenmu2)/(m+MSbp+MLbp);

