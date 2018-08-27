function dydt = odefunddQ1(t,y,A,L,w,m,g,MSbp,MLbp,FM,forcefai)
dydt = zeros(2,1);
dydt(1) = y(2);
coswt=1-0.5*(w*t)^2+(w*t)^4/(4*3*2)-(w*t)^6/(6*5*4*3*2)+(w*t)^8/(8*7*6*5*4*3*2);

FL=MSbp*g+MLbp*g-FM*cos(w*t-forcefai);
dydt(2) =(m*g-m*A*w^2*coswt+FL)*y(1)/((L-A+A*coswt)*(m+MSbp+MLbp));