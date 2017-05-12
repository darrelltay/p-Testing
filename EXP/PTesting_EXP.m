function [Result,PreResult] = PTesting_EXP(Sample,SaveDir)
%% This code need to work with QuickKSDNPValueCalculator,ExponentialFittingCompareResult
%%%%%% Sample is empirical data in Nx1 array
Sample = sort(Sample);
Sample = reshape(Sample,[max(size(Sample)),1]);
%% Estimate the parameter X_min and alpah, also determine KS-distance and Distribution noise. 
[WithoutAdjustment,AdjustedForBeta,AdjustedForBetaAndKSDist] = ExponentialFittingCompareResult(Sample');
BetaArraySample = [WithoutAdjustment.EXP_Beta, AdjustedForBeta.EXP_Beta, AdjustedForBetaAndKSDist.EXP_Beta];
XminArraySample = [WithoutAdjustment.EXP_Xmin, AdjustedForBeta.EXP_Xmin, AdjustedForBetaAndKSDist.EXP_Xmin];
KSDistArraySample = [WithoutAdjustment.EXP_KS_Dist, AdjustedForBeta.EXP_KS_Dist, AdjustedForBetaAndKSDist.EXP_KS_Dist];
DistNoiseArraySample = [WithoutAdjustment.EXP_DistNoise, AdjustedForBeta.EXP_DistNoise, AdjustedForBetaAndKSDist.EXP_DistNoise];
NumSampleArraySample = [WithoutAdjustment.EXP_NumSample, AdjustedForBeta.EXP_NumSample, AdjustedForBetaAndKSDist.EXP_NumSample];

Unadjusted.Parameter.Beta = BetaArraySample(1,1);
Unadjusted.Parameter.Xmin = XminArraySample(1,1);
Unadjusted.Parameter.KS = KSDistArraySample(1,1);
Unadjusted.Parameter.DN = DistNoiseArraySample(1,1);
Unadjusted.Parameter.SS = NumSampleArraySample(1,1);

Adjusted.Parameter.Beta = BetaArraySample(1,3);
Adjusted.Parameter.Xmin = XminArraySample(1,3);
Adjusted.Parameter.KS = KSDistArraySample(1,3);
Adjusted.Parameter.DN = DistNoiseArraySample(1,3);
Adjusted.Parameter.SS = NumSampleArraySample(1,3);

%% CSN p-testing without adjustment. Resimulate NumSimulation (1000) times of sample using unadjusted parameter. 
NumSimulation = 1000;

XminPWo = XminArraySample(1,1);
BetaPWo = BetaArraySample(1,1);
TruncatePWo = Sample(Sample<XminPWo);
SampleSize = size(Sample,1)-size(TruncatePWo,1);

SimulatedSample  = XminPWo -log(1-(1-0)*rand(SampleSize,NumSimulation))/BetaPWo;
SimulatedSample = sort([SimulatedSample;repmat(TruncatePWo,[1,NumSimulation])]);

[WithoutAdjustment,AdjustedForBeta,AdjustedForBetaAndKSDist] = ExponentialFittingCompareResult(SimulatedSample');
BetaArrayPWo = [WithoutAdjustment.EXP_Beta, AdjustedForBeta.EXP_Beta, AdjustedForBetaAndKSDist.EXP_Beta];
XminArrayPWo = [WithoutAdjustment.EXP_Xmin, AdjustedForBeta.EXP_Xmin, AdjustedForBetaAndKSDist.EXP_Xmin];
KSDistArrayPWo = [WithoutAdjustment.EXP_KS_Dist, AdjustedForBeta.EXP_KS_Dist, AdjustedForBetaAndKSDist.EXP_KS_Dist];
DistNoiseArrayPWo = [WithoutAdjustment.EXP_DistNoise, AdjustedForBeta.EXP_DistNoise, AdjustedForBetaAndKSDist.EXP_DistNoise];
NumSampleArrayPWo = [WithoutAdjustment.EXP_NumSample, AdjustedForBeta.EXP_NumSample, AdjustedForBetaAndKSDist.EXP_NumSample];

%% CSN p-testing without adjustment. Resimulate NumSimulation (1000) times of sample using unadjusted parameter. 
XminPW = XminArraySample(1,3);
BetaPW = BetaArraySample(1,3);
TruncatePW = Sample(Sample<XminPW);
SampleSize = size(Sample,1)-size(TruncatePW,1);

SimulatedSample  = XminPW -log(1-(1-0)*rand(SampleSize,NumSimulation))/BetaPW;
SimulatedSample = sort([SimulatedSample;repmat(TruncatePW,[1,NumSimulation])]);

[WithoutAdjustment,AdjustedForBeta,AdjustedForBetaAndKSDist] = ExponentialFittingCompareResult(SimulatedSample');
BetaArrayPW = [WithoutAdjustment.EXP_Beta, AdjustedForBeta.EXP_Beta, AdjustedForBetaAndKSDist.EXP_Beta];
XminArrayPW = [WithoutAdjustment.EXP_Xmin, AdjustedForBeta.EXP_Xmin, AdjustedForBetaAndKSDist.EXP_Xmin];
KSDistArrayPW = [WithoutAdjustment.EXP_KS_Dist, AdjustedForBeta.EXP_KS_Dist, AdjustedForBetaAndKSDist.EXP_KS_Dist];
DistNoiseArrayPW = [WithoutAdjustment.EXP_DistNoise, AdjustedForBeta.EXP_DistNoise, AdjustedForBetaAndKSDist.EXP_DistNoise];
NumSampleArrayPW = [WithoutAdjustment.EXP_NumSample, AdjustedForBeta.EXP_NumSample, AdjustedForBetaAndKSDist.EXP_NumSample];

%% Calculate p-value
KS_SSDependencePower = 0.4924;
DN_SSDependencePower = 0.4946;
%%% CSN unadjusted
KS_Sample_Unadj = KSDistArraySample(1,1);
DN_Sample_Unadj = DistNoiseArraySample(1,1)-1/sqrt(2);
SS_Sample_Unadj = NumSampleArraySample(1,1);
KS_CSN_Unadj = KSDistArrayPWo(:,1);
DN_CSN_Unadj = DistNoiseArrayPWo(:,1)-1/sqrt(2);
SS_CSN_Unadj = NumSampleArrayPWo(:,1);

KSPValue_CSN_Unadj = 100-100*sum(KS_CSN_Unadj.*(SS_CSN_Unadj.^KS_SSDependencePower)<KS_Sample_Unadj.*(SS_Sample_Unadj.^KS_SSDependencePower))/NumSimulation;
DNPValue_CSN_Unadj = 100-100*sum(DN_CSN_Unadj.*(SS_CSN_Unadj.^DN_SSDependencePower)<DN_Sample_Unadj.*(SS_Sample_Unadj.^DN_SSDependencePower))/NumSimulation;

Unadjusted.PValue.CSN_KS = KSPValue_CSN_Unadj;
Unadjusted.PValue.CSN_DN = DNPValue_CSN_Unadj;
%%% Inversion formulae unadjusted
[PKS,PDN] = QuickKSDNPValueCalculator(KS_Sample_Unadj,DN_Sample_Unadj+1/sqrt(2),SS_Sample_Unadj);
Unadjusted.PValue.IF_KS = 100-PKS;
Unadjusted.PValue.IF_DN = 100-PDN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CSN adjusted
KS_Sample_adj = KSDistArraySample(1,3);
DN_Sample_adj = DistNoiseArraySample(1,3)-1/sqrt(2);
SS_Sample_adj = NumSampleArraySample(1,3);
KS_CSN_adj = KSDistArrayPW(:,3);
DN_CSN_adj = DistNoiseArrayPW(:,3)-1/sqrt(2);
SS_CSN_adj = NumSampleArrayPW(:,3);

KSPValue_CSN_adj = 100-100*sum(KS_CSN_adj.*(SS_CSN_adj.^KS_SSDependencePower)<KS_Sample_adj.*(SS_Sample_adj.^KS_SSDependencePower))/NumSimulation;
DNPValue_CSN_adj = 100-100*sum(DN_CSN_adj.*(SS_CSN_adj.^DN_SSDependencePower)<DN_Sample_adj.*(SS_Sample_adj.^DN_SSDependencePower))/NumSimulation;

Adjusted.PValue.CSN_KS = KSPValue_CSN_adj;
Adjusted.PValue.CSN_DN = DNPValue_CSN_adj;
%%% Inversion formulae adjusted
[PKS,PDN] = QuickKSDNPValueCalculator(KS_Sample_adj,DN_Sample_adj+1/sqrt(2),SS_Sample_adj);
Adjusted.PValue.IF_KS = 100-PKS;
Adjusted.PValue.IF_DN = 100-PDN;

Result.Sample = Sample;
Result.Adjusted = Adjusted;
Result.Unadjusted = Unadjusted;
%% Organized Data with estimated parameters details
PreResult.Sample = Sample;
PreResult.BetaArraySample = BetaArraySample;
PreResult.XminArraySample = XminArraySample;
PreResult.KSDistArraySample = KSDistArraySample;
PreResult.DistNoiseArraySample = DistNoiseArraySample;
PreResult.NumSampleArraySample = NumSampleArraySample;

PreResult.BetaArrayPWo = BetaArrayPWo;
PreResult.XminArrayPWo = XminArrayPWo;
PreResult.KSDistArrayPWo = KSDistArrayPWo;
PreResult.DistNoiseArrayPWo = DistNoiseArrayPWo;
PreResult.NumSampleArrayPWo = NumSampleArrayPWo;

PreResult.BetaArrayPW = BetaArrayPW;
PreResult.XminArrayPW = XminArrayPW;
PreResult.KSDistArrayPW = KSDistArrayPW;
PreResult.DistNoiseArrayPW = DistNoiseArrayPW;
PreResult.NumSampleArrayPW = NumSampleArrayPW;

save([SaveDir '/Result_CSNTest.mat'],'Result','PreResult');