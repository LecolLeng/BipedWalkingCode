function dydt = odefunddQ4(t,y,A,L,w,g,d,T,m,MSbp,MLbp,FM,forcefai)
dydt = zeros(2,1);
dydt(1) = y(2);
coswt=1-0.5*(w*(t-T/2))^2+(w*(t-T/2))^4/(4*3*2)-(w*(t-T/2))^6/(6*5*4*3*2)+(w*(t-T/2))^8/(8*7*6*5*4*3*2);
FL=MSbp*g+MLbp*g-FM*cos(w*t-forcefai);
dydt(2) =-(d-y(1))*(-m*A*w^2*coswt+m*g+FL)/((L-A+A*coswt)*(m+MSbp+MLbp));

