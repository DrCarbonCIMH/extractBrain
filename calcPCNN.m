function [Fn,Ln,Tn,Yn] = calcPCNN(F,L,T,Y,M,S)
    K = convn(Y,M,'same');
    Fn = exp(-log(2)/0.3).*F + 0.01*K + S;
    Ln = exp(-log(2)/1).*L + 0.2.*K;
    Un = Fn.*(1+0.2*Ln);
    Tn = exp(-log(2)/10).*T + 20.*Y;
    Yn = double(Un>Tn);
end

