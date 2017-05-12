function [WithoutAdjustment,AdjustedForBeta,AdjustedForBetaAndKSDist] = ExponentialFittingCompareResult(Data)
warning('off','all');
options = optimoptions('lsqnonlin','Display','off');
Data = sort(Data,2);
EXP_Beta = zeros(size(Data,1),1);
EXP_Beta_Xmin = zeros(size(Data,1),1);
EXP_Beta_Dist = zeros(size(Data,1),1);
EXP_Beta_DistNoise = zeros(size(Data,1),1);
EXP_Beta_NumSample = zeros(size(Data,1),1);
EXP_Beta_PreKS_Dist = inf*ones(size(Data,1),1);

EXP_BetaAdj_Beta = zeros(size(Data,1),1);
EXP_BetaAdj_Xmin = zeros(size(Data,1),1);
EXP_BetaAdjKS_Dist = zeros(size(Data,1),1);
EXP_BetaAdj_DistNoise = zeros(size(Data,1),1);
EXP_BetaAdj_NumSample = zeros(size(Data,1),1);
EXP_BetaAdj_PreKS_Dist = inf*ones(size(Data,1),1);
 
EXP_BetaAdjKSAdj_Beta = zeros(size(Data,1),1);
EXP_BetaAdjKSAdj_Xmin = zeros(size(Data,1),1);
EXP_BetaAdjKSAdj_Dist = zeros(size(Data,1),1);
EXP_BetaAdjKSAdj_DistNoise = zeros(size(Data,1),1);
EXP_BetaAdjKSAdj_NumSample = zeros(size(Data,1),1);
EXP_BetaAdj_PreKS_DistAdj = inf*ones(size(Data,1),1);
    
XminArray = unique(Data);
XminArray = reshape(XminArray,[1,size(XminArray,1)*size(XminArray,2)]);
Ind = XminArray>=min(Data(:,end+1-min(size(Data,2),50)));  XminArray(:,Ind) = [];
if size(XminArray,2)>100
    XminArray = linspace(min(XminArray),max(XminArray),100);
end
for Xmin_i = 1:size(XminArray,2)-1
    Xmin = XminArray(1,Xmin_i);
    Ind =sum(Data.*(Data>=Xmin),2);
    XMean = Ind./sum(Data>=Xmin,2);
%% Parameter Estimation
    %%% Beta
    Beta = 1./(XMean-Xmin); XMax = max(Data')';
    %%% Beta Adjusted
    BetaFunction =@(B) abs(B.*(XMean-Xmin)+(B.*(XMax-XMean)+1).*exp(-B.*(XMax-Xmin))-1);
    BetaAdj = lsqnonlin(BetaFunction,Beta,0,inf,options);

%% KS-Distance Measure    
    Ind = Data>=Xmin;
    CDF_Th = cumsum(Ind')'; CDF_Th = CDF_Th./repmat(sum(Ind,2),[1,size(CDF_Th,2)]);
    %%% Without any adjustment
    CumF_Beta = Ind.*(1-exp(-repmat(Beta,[1,size(Data,2)]).*(Data-Xmin)));
%     UpperBound = 1-exp(-Beta.*(XMax-Xmin));    CumF_Beta = CumF_Beta./repmat(UpperBound,[1,size(Data,2)]);
    Dist = max(abs(CDF_Th-CumF_Beta)')';
    Dist(isnan(Dist))=inf;
    Weight = diff(CumF_Beta')';
    DistributionNoise = (Weight.^2).*(repmat(1./sum(Ind,2),[1,size(Weight,2)])./Weight-1).^2;
    DistributionNoise(isnan(DistributionNoise))=0;
    DistributionNoise = sqrt(sum(DistributionNoise,2)./sum(Weight.^2,2));
    EXP_Beta    = (Dist<EXP_Beta_PreKS_Dist).*Beta + (Dist>=EXP_Beta_PreKS_Dist).*EXP_Beta;
    EXP_Beta_Xmin    = (Dist<EXP_Beta_PreKS_Dist).*Xmin + (Dist>=EXP_Beta_PreKS_Dist).*EXP_Beta_Xmin;
    EXP_Beta_Dist = (Dist<EXP_Beta_PreKS_Dist).*Dist + (Dist>=EXP_Beta_PreKS_Dist).*EXP_Beta_Dist;
    EXP_Beta_DistNoise = (Dist<EXP_Beta_PreKS_Dist).*DistributionNoise + (Dist>=EXP_Beta_PreKS_Dist).*EXP_Beta_DistNoise;
    EXP_Beta_NumSample    = (Dist<EXP_BetaAdj_PreKS_DistAdj).*sum(Ind,2) + (Dist>=EXP_BetaAdj_PreKS_DistAdj).*EXP_Beta_NumSample;
    EXP_Beta_PreKS_Dist = min(EXP_Beta_PreKS_Dist,Dist);
    %%% with parameter estimation adjustment
    CumF_BetaAdj = Ind.*(1-exp(-repmat(BetaAdj,[1,size(Data,2)]).*(Data-Xmin)));
    Dist = max(abs(CDF_Th-CumF_BetaAdj)')';
    Dist(isnan(Dist))=inf;
    Weight = diff(CumF_BetaAdj')';
    DistributionNoise = (Weight.^2).*(repmat(1./sum(Ind,2),[1,size(Weight,2)])./Weight-1).^2;
    DistributionNoise(isnan(DistributionNoise))=0;
    DistributionNoise = sqrt(sum(DistributionNoise,2)./sum(Weight.^2,2));
    EXP_BetaAdj_Beta    = (Dist<EXP_BetaAdj_PreKS_Dist).*BetaAdj + (Dist>=EXP_BetaAdj_PreKS_Dist).*EXP_BetaAdj_Beta;
    EXP_BetaAdj_Xmin    = (Dist<EXP_BetaAdj_PreKS_Dist).*Xmin + (Dist>=EXP_BetaAdj_PreKS_Dist).*EXP_BetaAdj_Xmin;
    EXP_BetaAdjKS_Dist    = (Dist<EXP_BetaAdj_PreKS_Dist).*Dist + (Dist>=EXP_BetaAdj_PreKS_Dist).*EXP_BetaAdjKS_Dist;
    EXP_BetaAdj_DistNoise    = (Dist<EXP_BetaAdj_PreKS_Dist).*DistributionNoise + (Dist>=EXP_BetaAdj_PreKS_Dist).*EXP_BetaAdj_DistNoise;
    EXP_BetaAdj_NumSample    = (Dist<EXP_BetaAdj_PreKS_DistAdj).*sum(Ind,2) + (Dist>=EXP_BetaAdj_PreKS_DistAdj).*EXP_BetaAdj_NumSample;
    EXP_BetaAdj_PreKS_Dist = min(EXP_BetaAdj_PreKS_Dist,Dist);
    %%% with parameter estimation adjustment and rescale KS Dist
    UpperBound = 1-exp(-BetaAdj.*(XMax-Xmin));
    CumF_BetaAdjKSAdj = CumF_BetaAdj./repmat(UpperBound,[1,size(Data,2)]);
    Dist = max(abs(CDF_Th-CumF_BetaAdjKSAdj)')';
    Dist(isnan(Dist))=inf;
    Weight = diff(CumF_BetaAdjKSAdj')';
    DistributionNoise = (Weight.^2).*(repmat(1./sum(Ind,2),[1,size(Weight,2)])./Weight-1).^2;
    DistributionNoise(isnan(DistributionNoise))=0;
    DistributionNoise = sqrt(sum(DistributionNoise,2)./sum(Weight.^2,2));
    EXP_BetaAdjKSAdj_Beta    = (Dist<EXP_BetaAdj_PreKS_DistAdj).*BetaAdj + (Dist>=EXP_BetaAdj_PreKS_DistAdj).*EXP_BetaAdjKSAdj_Beta;
    EXP_BetaAdjKSAdj_Xmin    = (Dist<EXP_BetaAdj_PreKS_DistAdj).*Xmin + (Dist>=EXP_BetaAdj_PreKS_DistAdj).*EXP_BetaAdjKSAdj_Xmin;
    EXP_BetaAdjKSAdj_Dist    = (Dist<EXP_BetaAdj_PreKS_DistAdj).*Dist + (Dist>=EXP_BetaAdj_PreKS_DistAdj).*EXP_BetaAdjKSAdj_Dist;
    EXP_BetaAdjKSAdj_DistNoise    = (Dist<EXP_BetaAdj_PreKS_DistAdj).*DistributionNoise + (Dist>=EXP_BetaAdj_PreKS_DistAdj).*EXP_BetaAdjKSAdj_DistNoise;
    EXP_BetaAdjKSAdj_NumSample    = (Dist<EXP_BetaAdj_PreKS_DistAdj).*sum(Ind,2) + (Dist>=EXP_BetaAdj_PreKS_DistAdj).*EXP_BetaAdjKSAdj_NumSample;
    EXP_BetaAdj_PreKS_DistAdj = min(EXP_BetaAdj_PreKS_DistAdj,Dist);
end
WithoutAdjustment.EXP_Beta = EXP_Beta;
WithoutAdjustment.EXP_Xmin = EXP_Beta_Xmin;
WithoutAdjustment.EXP_KS_Dist = EXP_Beta_Dist;
WithoutAdjustment.EXP_DistNoise = EXP_Beta_DistNoise;
WithoutAdjustment.EXP_NumSample = EXP_Beta_NumSample;

AdjustedForBeta.EXP_Beta = EXP_BetaAdj_Beta;
AdjustedForBeta.EXP_Xmin = EXP_BetaAdj_Xmin;
AdjustedForBeta.EXP_KS_Dist = EXP_BetaAdjKS_Dist;
AdjustedForBeta.EXP_DistNoise = EXP_BetaAdj_DistNoise;
AdjustedForBeta.EXP_NumSample = EXP_BetaAdj_NumSample;

AdjustedForBetaAndKSDist.EXP_Beta =EXP_BetaAdjKSAdj_Beta;
AdjustedForBetaAndKSDist.EXP_Xmin = EXP_BetaAdjKSAdj_Xmin;
AdjustedForBetaAndKSDist.EXP_KS_Dist = EXP_BetaAdjKSAdj_Dist;
AdjustedForBetaAndKSDist.EXP_DistNoise = EXP_BetaAdjKSAdj_DistNoise;
AdjustedForBetaAndKSDist.EXP_NumSample = EXP_BetaAdjKSAdj_NumSample;
end
