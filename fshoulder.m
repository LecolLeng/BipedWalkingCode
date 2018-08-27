

function ffshoulder=fshoulder(m,M,B,H,k,c,w,t,g,fai)

ffshoulder=(m+M)*g+k*B*sin(w*t-fai)+c*B*w*cos(w*t-fai)-k*H*sin(w*t)-c*H*w*cos(w*t);