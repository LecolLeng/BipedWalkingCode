%�����ڼ�򱳰�����ʱ��ķ�ֵ
function [FMPack, forcefai]=FMPack(MSbp,MLbp,Babs,w,H,fai)  %fai�ǹ���k��c�ĺ���  
FMPack=sqrt((MLbp*Babs*w^2)^2+(MSbp*H*w^2)^2+2*MLbp*Babs*MSbp*H*w^4*cos(-fai));
maxfai=atan((MLbp*Babs-MSbp*H)/(MLbp*Babs+MSbp*H)*tan(-fai/2));
forcefai=fai/2-maxfai;