function fai=faical(m,k,c,w)

fai=atan((m*c*w*w*w)/(c*c*w*w+k*k-k*m*w*w));
if fai<0
    fai=pi+fai;
end