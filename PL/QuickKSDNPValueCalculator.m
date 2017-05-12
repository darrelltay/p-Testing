function [PKS,PDN] = QuickKSDNPValueCalculator(KS,DN,N)
options = optimoptions('lsqnonlin','Display','off','TolFun',1e-12);

ParaKS = [0.2744,0.1760,0.4924];
PKS = 100./(1+(exp(ParaKS(1,1))*KS.* (N.^(ParaKS(1,3)))).^(-1/ParaKS(1,2)));

EDN = sqrt( 1./(2-exp(-sqrt(5)./N)) .* (1-(1+N.^2).*exp(-N)) ./ (2+(2*N+N.^2).*exp(-N)) );

ParaDN = [0.4298,0.3023,0.4946];
Miu = DN-EDN;
if Miu<0
    PDNFun =@(PDN)   PDN.^(ParaDN(1,1))+ log(abs(Miu).*(N.^ParaDN(1,3))).*((50-PDN).^(ParaDN(1,2))) ;
    PDN = lsqnonlin(PDNFun,25,0,50,options);
elseif Miu==0
    PDN = 50;
else
    PDNFun =@(PDN) (100-PDN).^(ParaDN(1,1))+ log(abs(Miu).*(N.^ParaDN(1,3))).*((PDN-50).^(ParaDN(1,2)));
    PDN = lsqnonlin(PDNFun,75,50,100,options);
end
end