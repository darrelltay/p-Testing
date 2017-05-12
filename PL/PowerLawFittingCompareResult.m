function [WithoutAdjustment,AdjustedForAlpha,AdjustedForAlphaAndKSDist] = PowerLawFittingCompareResult(Data)
warning('off','all');
options = optimoptions('lsqnonlin','Display','off');
Data = sort(Data,2);
PL_Alpha = zeros(size(Data,1),1);
PL_Alpha_Xmin = zeros(size(Data,1),1);
PL_Alpha_Dist = zeros(size(Data,1),1);
PL_Alpha_DistNoise = zeros(size(Data,1),1);
PL_Alpha_NumSample = zeros(size(Data,1),1);
PL_Alpha_PreKS_Dist = inf*ones(size(Data,1),1);

PL_AlphaAdj_Alpha = zeros(size(Data,1),1);
PL_AlphaAdj_Xmin = zeros(size(Data,1),1);
PL_AlphaAdjKS_Dist = zeros(size(Data,1),1);
PL_AlphaAdj_DistNoise = zeros(size(Data,1),1);
PL_AlphaAdj_NumSample = zeros(size(Data,1),1);
PL_AlphaAdj_PreKS_Dist = inf*ones(size(Data,1),1);

PL_AlphaAdjKSAdj_Alpha = zeros(size(Data,1),1);
PL_AlphaAdjKSAdj_Xmin = zeros(size(Data,1),1);
PL_AlphaAdjKSAdj_Dist = zeros(size(Data,1),1);
PL_AlphaAdjKSAdj_DistNoise = zeros(size(Data,1),1);
PL_AlphaAdjKSAdj_NumSample = zeros(size(Data,1),1);
PL_AlphaAdj_PreKS_DistAdj = inf*ones(size(Data,1),1);
    
XminArray = unique(Data);
XminArray = reshape(XminArray,[1,size(XminArray,1)*size(XminArray,2)]);
Ind = XminArray>=min(Data(:,end+1-min(size(Data,2),50)));  XminArray(:,Ind) = [];
if size(XminArray,2)>100
    XminArray = linspace(min(XminArray),max(XminArray),100);
end
for Xmin_i = 1:size(XminArray,2)-1
    Xmin = XminArray(1,Xmin_i);
    Ind =sum(log(Data).*(Data>=Xmin),2);
    LogXPerXminMean = Ind./sum(Data>=Xmin,2)-log(Xmin);
%% Parameter Estimation
    %%% Alpha
    Alpha = 1+1./LogXPerXminMean; XMax = max(Data')';
    %%% Alpha Adjusted
    AlphaFunction =@(A) abs( (LogXPerXminMean+1./(1-A))+ (log(XMax./Xmin)+1./(A-1)-LogXPerXminMean).*(XMax./Xmin).^(1-A) );
    AlphaAdj = lsqnonlin(AlphaFunction,Alpha,2,inf,options);

%% KS-Distance Measure    
    Ind = Data>=Xmin;
    CDF_Th = cumsum(Ind')'; CDF_Th = CDF_Th./repmat(sum(Ind,2),[1,size(CDF_Th,2)]);
    %%% Without any adjustment
    CumF_Alpha = Ind.*(1-(Data./Xmin).^(1-repmat(Alpha,[1,size(Data,2)])));
    UpperBound = 1-(XMax./Xmin).^(1-Alpha);    CumF_Alpha = CumF_Alpha./repmat(UpperBound,[1,size(Data,2)]);
    Dist = max(abs(CDF_Th-CumF_Alpha)')';
    Dist(isnan(Dist))=inf;
    Weight = diff(CumF_Alpha')';
    DistributionNoise = (Weight.^2).*(repmat(1./sum(Ind,2),[1,size(Weight,2)])./Weight-1).^2;
    DistributionNoise(isnan(DistributionNoise))=0;
    DistributionNoise = sqrt(sum(DistributionNoise,2)./sum(Weight.^2,2));
    PL_Alpha    = (Dist<PL_Alpha_PreKS_Dist).*Alpha + (Dist>=PL_Alpha_PreKS_Dist).*PL_Alpha;
    PL_Alpha_Xmin    = (Dist<PL_Alpha_PreKS_Dist).*Xmin + (Dist>=PL_Alpha_PreKS_Dist).*PL_Alpha_Xmin;
    PL_Alpha_Dist = (Dist<PL_Alpha_PreKS_Dist).*Dist + (Dist>=PL_Alpha_PreKS_Dist).*PL_Alpha_Dist;
    PL_Alpha_DistNoise = (Dist<PL_Alpha_PreKS_Dist).*DistributionNoise + (Dist>=PL_Alpha_PreKS_Dist).*PL_Alpha_DistNoise;
    PL_Alpha_NumSample = (Dist<PL_Alpha_PreKS_Dist).*sum(Ind,2) + (Dist>=PL_Alpha_PreKS_Dist).*PL_Alpha_NumSample;
    PL_Alpha_PreKS_Dist = min(PL_Alpha_PreKS_Dist,Dist);
    %%% with parameter estimation adjustment
    CumF_AlphaAdj = Ind.*(1-(Data./Xmin).^(1-repmat(AlphaAdj,[1,size(Data,2)])));
    Dist = max(abs(CDF_Th-CumF_AlphaAdj)')';
    Dist(isnan(Dist))=inf;
    Weight = diff(CumF_AlphaAdj')';
    DistributionNoise = (Weight.^2).*(repmat(1./sum(Ind,2),[1,size(Weight,2)])./Weight-1).^2;
    DistributionNoise(isnan(DistributionNoise))=0;
    DistributionNoise = sqrt(sum(DistributionNoise,2)./sum(Weight.^2,2));
    PL_AlphaAdj_Alpha    = (Dist<PL_AlphaAdj_PreKS_Dist).*AlphaAdj + (Dist>=PL_AlphaAdj_PreKS_Dist).*PL_AlphaAdj_Alpha;
    PL_AlphaAdj_Xmin    = (Dist<PL_AlphaAdj_PreKS_Dist).*Xmin + (Dist>=PL_AlphaAdj_PreKS_Dist).*PL_AlphaAdj_Xmin;
    PL_AlphaAdjKS_Dist    = (Dist<PL_AlphaAdj_PreKS_Dist).*Dist + (Dist>=PL_AlphaAdj_PreKS_Dist).*PL_AlphaAdjKS_Dist;
    PL_AlphaAdj_DistNoise    = (Dist<PL_AlphaAdj_PreKS_Dist).*DistributionNoise + (Dist>=PL_AlphaAdj_PreKS_Dist).*PL_AlphaAdj_DistNoise;
    PL_AlphaAdj_NumSample    = (Dist<PL_AlphaAdj_PreKS_Dist).*sum(Ind,2) + (Dist>=PL_AlphaAdj_PreKS_Dist).*PL_AlphaAdj_NumSample;
    PL_AlphaAdj_PreKS_Dist = min(PL_AlphaAdj_PreKS_Dist,Dist);
    %%% with parameter estimation adjustment and rescale KS Dist
    UpperBound = 1-(XMax./Xmin).^(1-AlphaAdj);
    CumF_AlphaAdjKSAdj = CumF_AlphaAdj./repmat(UpperBound,[1,size(Data,2)]);
    Dist = max(abs(CDF_Th-CumF_AlphaAdjKSAdj)')';
    Dist(isnan(Dist))=inf;
    Weight = diff(CumF_AlphaAdjKSAdj')';
    DistributionNoise = (Weight.^2).*(repmat(1./sum(Ind,2),[1,size(Weight,2)])./Weight-1).^2;
    DistributionNoise(isnan(DistributionNoise))=0;
    DistributionNoise = sqrt(sum(DistributionNoise,2)./sum(Weight.^2,2));
    PL_AlphaAdjKSAdj_Alpha    = (Dist<PL_AlphaAdj_PreKS_DistAdj).*AlphaAdj + (Dist>=PL_AlphaAdj_PreKS_DistAdj).*PL_AlphaAdjKSAdj_Alpha;
    PL_AlphaAdjKSAdj_Xmin    = (Dist<PL_AlphaAdj_PreKS_DistAdj).*Xmin + (Dist>=PL_AlphaAdj_PreKS_DistAdj).*PL_AlphaAdjKSAdj_Xmin;
    PL_AlphaAdjKSAdj_Dist    = (Dist<PL_AlphaAdj_PreKS_DistAdj).*Dist + (Dist>=PL_AlphaAdj_PreKS_DistAdj).*PL_AlphaAdjKSAdj_Dist;
    PL_AlphaAdjKSAdj_DistNoise    = (Dist<PL_AlphaAdj_PreKS_DistAdj).*DistributionNoise + (Dist>=PL_AlphaAdj_PreKS_DistAdj).*PL_AlphaAdjKSAdj_DistNoise;
    PL_AlphaAdjKSAdj_NumSample    = (Dist<PL_AlphaAdj_PreKS_DistAdj).*sum(Ind,2) + (Dist>=PL_AlphaAdj_PreKS_DistAdj).*PL_AlphaAdjKSAdj_NumSample;
    PL_AlphaAdj_PreKS_DistAdj = min(PL_AlphaAdj_PreKS_DistAdj,Dist);
end
WithoutAdjustment.PL_Alpha = PL_Alpha;
WithoutAdjustment.PL_Xmin = PL_Alpha_Xmin;
WithoutAdjustment.PL_KS_Dist = PL_Alpha_Dist;
WithoutAdjustment.PL_DistNoise = PL_Alpha_DistNoise;
WithoutAdjustment.PL_NumSample = PL_Alpha_NumSample;

AdjustedForAlpha.PL_Alpha = PL_AlphaAdj_Alpha;
AdjustedForAlpha.PL_Xmin = PL_AlphaAdj_Xmin;
AdjustedForAlpha.PL_KS_Dist = PL_AlphaAdjKS_Dist;
AdjustedForAlpha.PL_DistNoise = PL_AlphaAdj_DistNoise;
AdjustedForAlpha.PL_NumSample = PL_AlphaAdj_NumSample;
    
AdjustedForAlphaAndKSDist.PL_Alpha =PL_AlphaAdjKSAdj_Alpha;
AdjustedForAlphaAndKSDist.PL_Xmin = PL_AlphaAdjKSAdj_Xmin;
AdjustedForAlphaAndKSDist.PL_KS_Dist = PL_AlphaAdjKSAdj_Dist;
AdjustedForAlphaAndKSDist.PL_DistNoise = PL_AlphaAdjKSAdj_DistNoise;
AdjustedForAlphaAndKSDist.PL_NumSample = PL_AlphaAdjKSAdj_NumSample;
end
