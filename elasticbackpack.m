function [FM,fai,forcefai,Babs]=elasticbackpack(Subsbackpack,A,w,g)


    MSbp=Subsbackpack(1);
    MLbp=Subsbackpack(2);
    K=Subsbackpack(3);
    C=Subsbackpack(4);

   fai=faical(MLbp,K,C,w);
   Babs=zhenfu(A,MLbp,K,C,w);
   [FM,forcefai]=FMPack(MSbp,MLbp,Babs,w,A,fai);