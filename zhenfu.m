function zhenfuB=zhenfu(H,m,k,c,w)

zhenfuB=(H*sqrt(k*k+c*c*w*w))/(sqrt((k-m*w*w)^2+c*c*w*w));