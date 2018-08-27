function QuaFunPara = FittingCurve(phase,Point1X,Point1Y,Point2X,Point2Y)
%README:
%------->phase1:the open direction is towards up/phase2:the open direction is towards down
% Fx = a*x^2+b*x+c;
if phase == 1
    A = [Point1X^2,Point1X,1;Point2X^2,Point2X,1;2*Point1X,1,0];
    B = [Point1Y;Point2Y;0];%%%参数波动上下限
    QuaFunPara =A^-1*B; %二次曲线的参数%           Fx = QuaFunPara(1)*x.^2+QuaFunPara(2)*x+QuaFunPara(3)
end
if phase ==2
    
    A = [Point1X^2,Point1X,1;Point2X^2,Point2X,1;2*Point2X,1,0];
    B = [Point1Y;Point2Y;0];%%%参数波动上下限
    QuaFunPara =A^-1*B; %二次曲线的参数%           Fx = QuaFunPara(1)*x.^2+QuaFunPara(2)*x+QuaFunPara(3)
end
end

