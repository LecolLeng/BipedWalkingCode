%作用在肩膀背包力振动时候的幅值
function [FMPack, forcefai]=FMPack(MSbp,MLbp,Babs,w,H,fai)  %fai是关于k，c的函数  
FMPack=sqrt((MLbp*Babs*w^2)^2+(MSbp*H*w^2)^2+2*MLbp*Babs*MSbp*H*w^4*cos(-fai));
maxfai=atan((MLbp*Babs-MSbp*H)/(MLbp*Babs+MSbp*H)*tan(-fai/2));
forcefai=fai/2-maxfai;